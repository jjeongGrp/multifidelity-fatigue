%% XCT Pore Subdomain Extraction and Visualization
% -------------------------------------------------------------------------
% This script extracts a cubic pore subdomain from a previously saved
% whole-domain XCT pore dataset.
%
% The input file, cubedomain.mat, is expected to contain a structure named
% `cube_domain`. In this script, `cube_domain` represents the whole processed
% XCT pore domain, not necessarily a cubic subdomain. The pore mask is stored
% as:
%
%       cube_domain.pore_mask
%
% with array ordering:
%
%       [Y, X, Z]
%
% where:
%       Y = image row direction
%       X = image column direction
%       Z = cropped slice direction
%
% The script performs the following operations:
%
%   1. Loads the whole processed XCT pore mask from cubedomain.mat.
%   2. Reads voxel-size and slice metadata from cube_domain when available.
%   3. Defines a cubic subdomain using:
%          - cube_vox  : number of voxels per cube edge
%          - center_mm : cube center in source-domain physical coordinates, mm
%   4. Converts the physical cube center to voxel indices using the
%      voxel-center coordinate convention.
%   5. Extracts the requested cube from the whole-domain pore mask.
%   6. Applies zero-padding if the requested cube extends outside the
%      available XCT domain.
%   7. Saves the extracted subdomain to a new MAT file.
%   8. Generates a 3D pore visualization of the extracted cube.
%
% Coordinate convention:
%
%   The physical coordinates are based on voxel centers:
%
%       X_mm = (X_idx - 0.5) * voxel_mm
%       Y_mm = (Y_idx - 0.5) * voxel_mm
%       Z_mm = (Z_idx - 0.5) * voxel_mm
%
%   Therefore, the voxel-center index corresponding to a physical coordinate is:
%
%       X_idx = X_mm / voxel_mm + 0.5
%       Y_idx = Y_mm / voxel_mm + 0.5
%       Z_idx = Z_mm / voxel_mm + 0.5
%
%   In the code, these index values are rounded to the nearest voxel center.
%
%   The Z-coordinate is measured in the cropped whole-domain coordinate
%   system. Thus, Z = 0 corresponds to the cropped plane associated with
%   slice_first_abs in the original XCT stack.
%
% Input:
%
%   cubedomain.mat
%       MAT-file containing:
%
%           cube_domain.pore_mask
%               Logical 3D pore mask for the whole processed XCT domain.
%               true  = pore
%               false = solid/non-pore
%
%           cube_domain.voxel_um or cube_domain.dx_mm
%               Voxel size.
%
%           cube_domain.slice_first_abs
%               Absolute slice number in the original TIFF stack
%               corresponding to cropped Z = 0.
%
%           cube_domain.slice_last_abs
%               Last absolute slice number used in the processed domain.
%
% Output:
%
%   CubeXXXmm_Nvox_Ngrid_PoreDomain.mat
%       MAT-file containing the structure `subdomain`, where:
%
%           subdomain.pore_mask
%               Extracted cubic pore mask, logical array.
%               true  = pore
%               false = solid/non-pore
%
%           subdomain.phase
%               Phase map for the extracted cube.
%               1 = solid
%               2 = pore
%
%           subdomain.cube_mm
%               Physical edge length of the cube in mm.
%
%           subdomain.cube_vox
%               Number of voxels per cube edge.
%
%           subdomain.grid_n
%               Number of grid points per cube edge, equal to cube_vox + 1.
%
%           subdomain.center_mm_cropped
%               Requested cube center in source-domain cropped physical
%               coordinates, mm.
%
%           subdomain.center_idx_cropped
%               Requested cube center in source-domain cropped voxel-index
%               coordinates.
%
%           subdomain.req_window_cropped
%               Requested extraction window in source-domain cropped
%               voxel-index space.
%
%           subdomain.clamp_window_cropped
%               Actual extraction window after clamping to available data.
%
%           subdomain.zero_padded
%               Boolean flag indicating whether zero-padding was applied.
%
% Notes:
%
%   - The extracted cube and the source pore mask both use array ordering:
%
%         [Y, X, Z]
%
%   - The extracted cube is visualized in a local coordinate system:
%
%         X = 0, Y = 0, Z = 0
%
%     at the lower corner of the extracted subdomain.
%
%   - The visualization coordinates do not include the source-domain offset.
%     The source-domain location of the cube is stored in the requested and
%     clamped window metadata.
%
%   - Pores are visualized as a gray isosurface with optional blue sampled
%     pore points.
%
%   - If the requested subdomain extends outside the available whole-domain
%     XCT data, missing regions are filled with false values, corresponding
%     to solid/non-pore material.
%
%   - The current extraction formula assumes cube_vox is even.
%
% Authors:       Juyoung Jeong, Veera Sundararaghavan
% Affiliation:   Department of Aerospace Engineering, University of Michigan,
%                Ann Arbor, MI 48109, USA.
% Repository:    https://github.com/jjeongGrp/multifidelity-fatigue
%
% Acknowledgment:
%   This work was supported by the Defense Advanced Research Projects Agency
%   (DARPA) SURGE program under Cooperative Agreement No. HR0011-25-2-0009,
%   'Predictive Real-time Intelligence for Metallic Endurance (PRIME)'.
%
% License:       MIT License. See LICENSE file in the repository.
% -------------------------------------------------------------------------
clear; close all;

%% ===== LOAD WHOLE DOMAIN FROM cubedomain.mat =====
S = load('cubedomain.mat');

if ~isfield(S, 'cube_domain')
    error('cubedomain.mat does not contain variable cube_domain.');
end

whole_domain = S.cube_domain;

if ~isfield(whole_domain, 'pore_mask')
    error('cube_domain does not contain field pore_mask.');
end

% Extract the logical pore mask from the loaded structure.
% Array order is [Y, X, Z]:
%   Y = image rows
%   X = image columns
%   Z = cropped slice index
%
% true  = pore
% false = solid
pore_volume = whole_domain.pore_mask;

[img_height, img_width, num_slices] = size(pore_volume);

% ===== READ VOXEL SIZE AND SOURCE-DOMAIN METADATA =====
% Read voxel size and slice metadata from the saved structure when possible.
% If the metadata fields are missing, default values are used.
if isfield(whole_domain, 'voxel_um')
    voxel_um = whole_domain.voxel_um;
elseif isfield(whole_domain, 'dx_mm')
    voxel_um = whole_domain.dx_mm * 1000;
else
    voxel_um = 5.5;
    warning('voxel_um not found in cubedomain.mat. Using default voxel_um = %.3f.', voxel_um);
end

voxel_mm = voxel_um * 1e-3;

if isfield(whole_domain, 'slice_first_abs')
    slice_first = whole_domain.slice_first_abs;
else
    slice_first = 1;
    warning('slice_first_abs not found. Assuming slice_first = 1.');
end

if isfield(whole_domain, 'slice_last_abs')
    slice_last = whole_domain.slice_last_abs;
else
    slice_last = slice_first + num_slices - 1;
end

% ===== DEFINE CUBIC SUBDOMAIN PARAMETERS =====
% Number of voxels along each cube edge.
% This extraction formula assumes cube_vox is even so that half = cube_vox/2
% is an integer and the requested index window has exactly cube_vox voxels.
cube_vox = 546;
if mod(cube_vox, 2) ~= 0
    error('cube_vox must be even for this centered extraction method.');
end
grid_n   = cube_vox + 1;
cube_mm  = cube_vox * voxel_mm;
dx_mm    = cube_mm / cube_vox;

% Physical center of the desired cube in mm.
% These coordinates are defined in the source whole-domain cropped coordinate
% system, not in the local coordinate system of the extracted cube.
%
% In particular, Z = 0 corresponds to the cropped plane associated with
% slice_first_abs from the original XCT stack.
center_mm = [2.750, 2.750, 2.000];

% Confirm that the cube voxel spacing is identical to the XCT voxel spacing.
if abs(dx_mm - voxel_mm) > 1e-9
    error('Mismatch: cube dx_mm=%.12f but voxel_mm=%.12f.', dx_mm, voxel_mm);
end

% ===== CONVERT PHYSICAL CUBE CENTER TO VOXEL INDICES =====
% Convert the requested physical cube center into voxel-center indices in
% the cropped whole-domain coordinate system.

cx = round(center_mm(1) / voxel_mm + 0.5);
cy = round(center_mm(2) / voxel_mm + 0.5);
cz = round(center_mm(3) / voxel_mm + 0.5);

half = cube_vox / 2;

% ===== COMPUTE REQUESTED EXTRACTION WINDOW =====
% Define the requested cubic window in cropped voxel-index space.
% The requested window may extend outside the available XCT domain.
req_x_start = cx - half + 1;
req_x_end   = cx + half;

req_y_start = cy - half + 1;
req_y_end   = cy + half;

req_z_start = cz - half + 1;
req_z_end   = cz + half;

% ===== CLAMP EXTRACTION WINDOW TO AVAILABLE DATA =====
% Clamp the requested window so that all indexing operations remain inside
% the available whole-domain pore mask.
clamp_x_start = max(1, req_x_start);
clamp_x_end   = min(img_width, req_x_end);

clamp_y_start = max(1, req_y_start);
clamp_y_end   = min(img_height, req_y_end);

clamp_z_start = max(1, req_z_start);
clamp_z_end   = min(num_slices, req_z_end);

% ===== COMPUTE BUFFER INDICES FOR ZERO-PADDED OUTPUT CUBE =====
% If the requested cube extends outside the XCT domain, only the available
% portion is copied into the output cube. The remaining voxels stay false.
% In the pore mask, false means non-pore. 
%
% Therefore, padded regions are treated as solid/non-pore material.
buf_x_start = clamp_x_start - req_x_start + 1;
buf_x_end   = clamp_x_end   - req_x_start + 1;

buf_y_start = clamp_y_start - req_y_start + 1;
buf_y_end   = clamp_y_end   - req_y_start + 1;

buf_z_start = clamp_z_start - req_z_start + 1;
buf_z_end   = clamp_z_end   - req_z_start + 1;

is_padded = (req_x_start < 1) || (req_x_end > img_width)  || ...
            (req_y_start < 1) || (req_y_end > img_height) || ...
            (req_z_start < 1) || (req_z_end > num_slices);


% ===== EXTRACT PORE SUBDOMAIN FROM WHOLE XCT DOMAIN =====
% Extract the clamped pore data from the whole-domain mask and insert it into
% the zero-initialized cubic output array.
%
% Both arrays use [Y, X, Z] ordering:
%
%   pore_volume : whole processed XCT pore mask
%   cube_pores  : extracted cubic pore subdomain

cube_pores = false(cube_vox, cube_vox, cube_vox);

cube_pores(buf_y_start:buf_y_end, ...
           buf_x_start:buf_x_end, ...
           buf_z_start:buf_z_end) = ...
    pore_volume(clamp_y_start:clamp_y_end, ...
                clamp_x_start:clamp_x_end, ...
                clamp_z_start:clamp_z_end);

cube_porosity = 100 * nnz(cube_pores) / numel(cube_pores);


% ===== DEFINE LOCAL SUBDOMAIN GRID AND VOXEL-CENTER COORDINATES =====
% Grid points define the cube boundaries and have length cube_vox + 1.
% Voxel-center coordinates define the centers of the cube voxels and have
% length cube_vox.
%
% The local coordinate system starts at the lower corner of the extracted cube.
xg_mm = linspace(0, cube_mm, grid_n);
yg_mm = linspace(0, cube_mm, grid_n);
zg_mm = linspace(0, cube_mm, grid_n);

xc_mm = ((1:cube_vox) - 0.5) * dx_mm;
yc_mm = ((1:cube_vox) - 0.5) * dx_mm;
zc_mm = ((1:cube_vox) - 0.5) * dx_mm;

% ===== BUILD OUTPUT SUBDOMAIN STRUCTURE =====
% Store the extracted pore mask, phase map, grid information, voxel size,
% extraction metadata, and source-domain metadata in a single structure.
subdomain = struct();

subdomain.pore_mask = cube_pores;                    % logical mask: true = pore, false = solid
subdomain.phase     = uint8(cube_pores) + uint8(1);  % phase map: 1 = solid, 2 = pore

subdomain.domain_type = 'subdomain_extracted_from_cubedomain';
subdomain.array_order = '[Y, X, Z]';
subdomain.phase_labels = struct('solid', uint8(1), 'pore', uint8(2));

subdomain.cube_mm  = cube_mm;
subdomain.cube_vox = cube_vox;
subdomain.grid_n   = grid_n;

subdomain.voxel_um = voxel_um;
subdomain.voxel_mm = voxel_mm;
subdomain.dx_mm    = dx_mm;

subdomain.center_mm_cropped  = center_mm;
subdomain.center_idx_cropped = [cx, cy, cz];

subdomain.req_window_cropped = [req_x_start, req_x_end; ...
                                req_y_start, req_y_end; ...
                                req_z_start, req_z_end];

subdomain.clamp_window_cropped = [clamp_x_start, clamp_x_end; ...
                                  clamp_y_start, clamp_y_end; ...
                                  clamp_z_start, clamp_z_end];

subdomain.buffer_window = [buf_x_start, buf_x_end; ...
                           buf_y_start, buf_y_end; ...
                           buf_z_start, buf_z_end];

subdomain.zero_padded = is_padded;

subdomain.x_grid_mm = xg_mm;
subdomain.y_grid_mm = yg_mm;
subdomain.z_grid_mm = zg_mm;

subdomain.x_center_mm = xc_mm;
subdomain.y_center_mm = yc_mm;
subdomain.z_center_mm = zc_mm;

subdomain.pore_voxels = nnz(cube_pores);
subdomain.total_voxels = numel(cube_pores);
subdomain.porosity_percent = cube_porosity;

subdomain.source_file = 'cubedomain.mat';
subdomain.source_slice_first_abs = slice_first;
subdomain.source_slice_last_abs  = slice_last;

% ===== SAVE EXTRACTED SUBDOMAIN =====
% Save the extracted pore cube and all relevant metadata to a new MAT-file.
cube_mm_str = strrep(sprintf('%.3f', cube_mm), '.', 'p');
out_file = sprintf('Cube%smm_%dvox_%dgrid_PoreDomain.mat', ...
    cube_mm_str, cube_vox, grid_n);

fprintf('\nSaving extracted subdomain to MAT file...\n');
save(out_file, 'subdomain', '-v7.3');
fprintf('Subdomain saved: %s\n', out_file);





%% ===== 3D VISUALIZATION OF EXTRACTED PORE SUBDOMAIN =====
% Visualize the extracted cubic pore domain using:
%   - gray pore isosurface
%   - optional blue pore point cloud
%   - cube-local physical coordinates in micrometers
mm2um  = 1000;
cube_um = cube_mm * mm2um;
dx_um   = dx_mm   * mm2um;

% Downsampling factor used only for visualization.
% The saved subdomain remains full resolution.
ds_cube = 1;  % 1 = full, 2 = faster, 3 or 4 = much faster
if ds_cube == 1
    cube_viz = cube_pores;
else
    if exist('imresize3', 'file')
        cube_viz = imresize3(cube_pores, 1 / ds_cube, 'nearest') > 0;
    else
        cube_viz = cube_pores(1:ds_cube:end, ...
                              1:ds_cube:end, ...
                              1:ds_cube:end);
    end
end

[nYv, nXv, nZv] = size(cube_viz);

% Generate pore isosurface
% Smooth the logical pore mask slightly before extracting the isosurface.
% This improves the rendered surface appearance but does not modify the
% saved pore mask.
cube_smooth = smooth3(single(cube_viz), 'box', 3);
[faces_cube, verts_cube] = isosurface(cube_smooth, 0.5);

if isempty(faces_cube) || isempty(verts_cube)

    warning('No isosurface generated for cube. Pores may be too sparse after downsampling.');

else

    % Clean duplicate mesh vertices
    % Remove duplicate isosurface vertices before mesh reduction.
    [verts_cube, ~, ic] = unique(round(verts_cube * 1000) / 1000, ...
        'rows', 'stable');

    faces_cube = ic(faces_cube);

    % Reduce mesh complexity
	% Reduce the number of triangular faces to improve rendering speed and
    % decrease figure-file size.



    f0 = size(faces_cube, 1);

    if f0 > 200000
        red = 0.05;
    elseif f0 > 80000
        red = 0.08;
    else
        red = 0.15;
    end

    [faces_cube, verts_cube] = reducepatch(faces_cube, verts_cube, red);


    % Map isosurface vertices to cube-local physical coordinates
    % isosurface returns vertices in visualization-grid index coordinates.
    % These are converted back to original cube voxel coordinates and then to
    % physical coordinates in micrometers.
    %
    % The plotted coordinate system is local to the extracted cube:
    %   X = 0 at the lower-X cube face
    %   Y = 0 at the lower-Y cube face
    %   Z = 0 at the lower-Z cube face
    verts_cube_orig = (verts_cube - 1) * ds_cube + 1;

    verts_cube_um = zeros(size(verts_cube_orig));
    verts_cube_um(:, 1) = (verts_cube_orig(:, 1) - 0.5) * dx_um;   % X um
    verts_cube_um(:, 2) = (verts_cube_orig(:, 2) - 0.5) * dx_um;   % Y um
    verts_cube_um(:, 3) = (verts_cube_orig(:, 3) - 0.5) * dx_um;   % Z um

    % Create figure
    fig_cube = figure; clf;

    ax = axes(fig_cube);
    hold(ax, 'on');
    set(ax, 'SortMethod', 'depth');

    %% Isosurface patch
    h_cube = patch('Vertices', verts_cube_um, ...
        'Faces', faces_cube, ...
        'FaceColor', [0.5 0.5 0.5], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.6, ...
        'AmbientStrength', 0.5, ...
        'DiffuseStrength', 0.6, ...
        'SpecularStrength', 0.3, ...
        'FaceLighting', 'flat', ...
        'BackFaceLighting', 'unlit', ...
        'SpecularExponent', 5);

    set(h_cube, 'DiffuseStrength', 0.6, 'SpecularColorReflectance', 0);

    try
        isonormals(cube_smooth, h_cube);
    catch
        warning('isonormals failed. Continuing without smoothed normals.');
    end

    % Add sampled pore point cloud overlay
	% Plot a randomly sampled subset of pore voxels to help show the internal
    % pore distribution. This is for visualization only.
    max_pts = 40000;
    lin_p = find(cube_viz);

    if ~isempty(lin_p)

        if numel(lin_p) > max_pts
            lin_p = lin_p(randperm(numel(lin_p), max_pts));
        end

        [py, px, pz] = ind2sub(size(cube_viz), lin_p);

        px_um = ((px - 1) * ds_cube + 1 - 0.5) * dx_um;
        py_um = ((py - 1) * ds_cube + 1 - 0.5) * dx_um;
        pz_um = ((pz - 1) * ds_cube + 1 - 0.5) * dx_um;

        h_pts = scatter3(px_um, py_um, pz_um, ...
            5, [0.1 0.3 0.7], 'filled', ...
            'MarkerFaceAlpha', 0.2);

    end

    % Draw cube bounding box
	% Add a black wireframe box showing the physical extent of the extracted cube.
    bx = [0 cube_um];
    by = [0 cube_um];
    bz = [0 cube_um];

    box_ex = [bx(1) bx(2); bx(1) bx(2); bx(1) bx(2); bx(1) bx(2); ...
              bx(1) bx(1); bx(2) bx(2); bx(1) bx(1); bx(2) bx(2); ...
              bx(1) bx(1); bx(2) bx(2); bx(1) bx(1); bx(2) bx(2)];

    box_ey = [by(1) by(1); by(2) by(2); by(1) by(1); by(2) by(2); ...
              by(1) by(2); by(1) by(2); by(1) by(2); by(1) by(2); ...
              by(1) by(1); by(1) by(1); by(2) by(2); by(2) by(2)];

    box_ez = [bz(1) bz(1); bz(1) bz(1); bz(2) bz(2); bz(2) bz(2); ...
              bz(1) bz(1); bz(1) bz(1); bz(2) bz(2); bz(2) bz(2); ...
              bz(1) bz(2); bz(1) bz(2); bz(1) bz(2); bz(1) bz(2)];

    for e = 1:size(box_ex, 1)
        plot3(box_ex(e, :), box_ey(e, :), box_ez(e, :), ...
            'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end

    % Axes formatting
    axis(ax, 'vis3d');
    daspect([1 1 1]);
    xlim([0 cube_um]);
    ylim([0 cube_um]);
    zlim([0 cube_um]);
    num_ticks = 6;
    xticks(linspace(0, cube_um, num_ticks));
    yticks(linspace(0, cube_um, num_ticks));
    zticks(linspace(0, cube_um, num_ticks));
    hx = xlabel('X_1 (\mum)', 'FontSize', 12, 'FontWeight', 'bold');
    hy = ylabel('X_2 (\mum)', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('X_3 (\mum)', 'FontSize', 12, 'FontWeight', 'bold');
    set(hx, 'Rotation', 21);
    set(hy, 'Rotation', -31);
    delete(findall(ax, 'Type', 'light'));

    camlight(ax, 'headlight');
    camlight(ax, 'right');
    camlight(ax, 'left');
    lighting(ax, 'gouraud');
    set(ax, 'FontSize', 12, 'GridAlpha', 0.3, 'LineWidth', 1.2);
    view(ax, -37.5, 28);
    camzoom(0.55);

end