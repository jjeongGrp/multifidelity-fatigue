%% PRISMS-Plasticity Subdomain Mesh and Orientation Export
% -------------------------------------------------------------------------
% This script generates PRISMS-Plasticity input files for a selected
% microstructural subdomain extracted from a voxelized polycrystalline
% microstructure.
%
% The script loads a previously generated SS316L polycrystal dataset,
% extracts a user-defined rectangular subdomain, removes pore elements from
% the finite element mesh, remaps global grain IDs to compact local grain IDs,
% and writes:
%
%   - a Gmsh-format hexahedral mesh file, micro.msh
%   - a Rodrigues-vector orientation file, ori.txt
%   - a grain ID mapping file, grain_id_map_subdomain.txt
%
% The exported mesh is intended for crystal-plasticity simulations in
% PRISMS-Plasticity. Each solid voxel in the selected subdomain becomes one
% 8-node hexahedral element. Pore voxels are skipped and therefore are not
% included as mesh elements.
%
% Main workflow:
%   1. Load the full microstructure and orientation data.
%   2. Define the physical subdomain bounds in X, Y, and Z.
%   3. Select voxels fully contained within the requested subdomain.
%   4. Optionally add fully solid padding layers at the top Z face.
%   5. Build a compact local grain-ID map for grains present in the subdomain.
%   6. Generate nodal coordinates with the subdomain origin shifted to zero.
%   7. Write the Gmsh 2.2 mesh file using 8-node hexahedral elements.
%   8. Skip pore elements during mesh export.
%   9. Write grain orientations for only the grains present in the subdomain.
%  10. Write a local-to-global grain ID mapping file.
%
% Inputs:
%   - xct_poly_params_SS316L_full.mat
%       MAT-file containing the full microstructure and orientation data,
%       expected to include:
%         * phase      : 3D voxel-wise phase/grain ID map
%         * xvec       : X node coordinate vector, mm
%         * yvec       : Y node coordinate vector, mm
%         * zvec       : Z node coordinate vector, mm
%         * N_grains   : number of solid grain phases in the full domain
%         * ori_vec    : full-domain grain orientations as Rodrigues vectors
%
% User-defined settings:
%   - x_range, y_range, z_range
%       Physical bounds of the exported subdomain, in mm.
%
%   - pad_z_top
%       Number of additional fully solid layers to append at the top Z face.
%       These layers are filled using the nearest solid grain ID from each
%       vertical column to avoid pore elements in the padding region.
%
% Outputs:
%   - micro.msh
%       Gmsh version 2.2 ASCII mesh file containing 8-node hexahedral elements.
%       Each exported solid voxel becomes one element. Pore voxels are omitted.
%
%   - ori.txt
%       PRISMS-Plasticity orientation file containing local grain IDs and
%       Rodrigues vectors:
%
%         local_grain_id, r_x, r_y, r_z
%
%   - grain_id_map_subdomain.txt
%       Text file mapping compact local grain IDs used in the exported mesh to
%       the original global grain IDs from the full microstructure.
%
% Phase convention:
%   - Solid grains in the full domain are labeled 1:N_grains.
%   - The global pore phase is assumed to be N_grains + 1.
%   - In the exported subdomain, solid grains are remapped to compact local IDs:
%
%         1:N_grains_sub
%
%     and the local pore ID is:
%
%         pore_id_sub = N_grains_sub + 1
%
%     However, pore elements are skipped and are not written to micro.msh.
%
% Coordinate convention:
%   - Input coordinate vectors xvec, yvec, and zvec are node coordinates in mm.
%   - The selected voxel/elements are defined using voxel left-node coordinates.
%   - Exported mesh coordinates are shifted so that the minimum corner of the
%     selected subdomain is at the origin.
%   - Output mesh coordinates remain in mm.
%
% Mesh convention:
%   - Element type:
%       Gmsh element type 5 = 8-node hexahedron.
%
%   - Element tags:
%       The physical and elementary tags are both set to the local grain ID.
%
% Notes:
%   - The selected ranges should align with the underlying voxel grid spacing
%     if exact subdomain dimensions are desired.
%   - Voxels are selected only when their left-node coordinate lies within:
%
%         [range_min, range_max)
%
%     so the voxel is fully contained in the requested subdomain.
%   - If pad_z_top > 0, the added top layers are forced to be solid to simplify
%     downstream boundary conditions or loading setup.
%
% Authors:       Juyoung Jeong, Veera Sundararaghavan
% Affiliation:   Department of Aerospace Engineering, University of Michigan,
%                Ann Arbor, MI 48109, USA.
% Repository:    https://github.com/jjeongGrp/multifidelity-fatigue
%
% Acknowledgment:
%   This work was supported by the Defense Advanced Research Projects Agency
%   (DARPA) SURGE program under Cooperative Agreement No. HR0011-25-2-0009,
%   "Predictive Real-time Intelligence for Metallic Endurance (PRIME)".
%
% License:       MIT License. See LICENSE file in the repository.
% -------------------------------------------------------------------------
clear; close all;

%% ---------------------- LOAD FULL MICROSTRUCTURE DATA --------------------
% Load the full-domain microstructure and orientation data generated by the
% preprocessing/microstructure-generation workflow.
%
% Expected variables include:
%   phase
%     3D voxel-wise phase map with dimensions [Ny, Nx, Nz].
%     Solid grains are labeled 1:N_grains and pores are labeled N_grains+1.
%
%   xvec, yvec, zvec
%     Node coordinate vectors in mm. Each vector has one more entry than the
%     corresponding voxel count.
%
%   N_grains
%     Number of solid grain phases in the full domain.
%
%   ori_vec
%     Full-domain Rodrigues orientation vectors arranged as:
%       [r1_x; r1_y; r1_z; r2_x; r2_y; r2_z; ...]
load('xct_poly_params_SS316L_full.mat');

%% ---------------------- EXPORT SETTINGS ---------------------------------
% Output file names for PRISMS-Plasticity input.
%
% micro.msh:
%   Gmsh 2.2 ASCII mesh file containing 8-node hexahedral elements.
%
% ori.txt:
%   Orientation file containing local grain IDs and Rodrigues vectors.
msh_filename = 'micro.msh';
ori_filename = 'ori.txt';

disp(['Writing GMSH mesh file: ', msh_filename])

% --- Optional top Z padding -----------------------------------------------
% Number of additional fully solid voxel layers appended to the top Z face
% of the selected subdomain.
%
% Purpose:
%   A fully solid top pad can simplify boundary condition application or avoid
%   loading directly on a porous top surface in downstream PRISMS-Plasticity
%   simulations.
%
% Behavior:
%   For each X-Y column, the added pad layers inherit the nearest solid grain
%   ID found by scanning downward from the top of the selected core subdomain.
%   This guarantees that the added pad region contains no pore elements.
%
% Set pad_z_top = 0 to disable this feature.
pad_z_top = 0;    % number of fully solid padding layers at top Z face

% Working copy of the full-domain phase map.
% Array convention:
%   phase3d(y_index, x_index, z_index)
phase3d = phase;


%% ---------------------- SELECT SUBDOMAIN -------------------------------
% Define the physical bounds of the PRISMS-Plasticity subdomain.
%
% Units:
%   All coordinate ranges are in millimeters.
%
% Convention:
%   x_range = [xmin, xmax]
%   y_range = [ymin, ymax]
%   z_range = [zmin, zmax]
%
% Voxels are selected using their left/lower node coordinates. A voxel is
% included if its left-node coordinate lies in:
%
%   [range_min, range_max)
%
% This selection rule avoids including voxels that extend beyond the requested
% upper bound, assuming the requested range aligns with the voxel grid.
x_range = [0.0000, 0.4950]; % mm
y_range = [1.2100, 1.7050]; % mm
z_range = [2.0515, 2.5465]; % mm

%% ---------------------- GRID AND VOXEL COORDINATES ----------------------
% phase3d is voxel-based and stored with array dimensions:
%
%   phase3d(y_index, x_index, z_index)
%
% Therefore:
%   Ny_vox = number of voxels in Y
%   Nx_vox = number of voxels in X
%   Nz_vox = number of voxels in Z
Ny_vox = size(phase3d, 1);
Nx_vox = size(phase3d, 2);
Nz_vox = size(phase3d, 3);

% xvec, yvec, and zvec are node coordinate vectors.
% For N voxels in one direction, the corresponding node vector has N+1 entries.
%
% The left-node coordinate of voxel i is xvec(i), so the final node is removed
% to obtain one coordinate per voxel.
xvec_vox = xvec(1:end-1);
yvec_vox = yvec(1:end-1);
zvec_vox = zvec(1:end-1);

% Grid spacing in each coordinate direction.
dx = xvec(2) - xvec(1);
dy = yvec(2) - yvec(1);
dz = zvec(2) - zvec(1);

% Voxel-center coordinate vectors.
% These are not required for the Gmsh node coordinates, but they are useful for
% diagnostics or possible future selection/visualization operations.
xvec_ctr = xvec_vox + 0.5*dx;
yvec_ctr = yvec_vox + 0.5*dy;
zvec_ctr = zvec_vox + 0.5*dz;


%% ---------------------- FIND SUBDOMAIN VOXEL INDICES --------------------
% Floating-point tolerance used when comparing physical coordinates to the
% requested subdomain bounds.
tol = 1e-9;

% Select voxels whose left-node coordinate lies within the requested range.
%
% The upper-bound condition uses < range_max rather than <= range_max so that
% voxels are selected by their lower/left node and remain fully inside the
% requested physical domain.
ix = find(xvec_vox >= x_range(1) - tol & xvec_vox < x_range(2) - tol);
iy = find(yvec_vox >= y_range(1) - tol & yvec_vox < y_range(2) - tol);
iz = find(zvec_vox >= z_range(1) - tol & zvec_vox < z_range(2) - tol);

% Safety clamp to valid voxel index bounds.
ix = ix(ix >= 1 & ix <= Nx_vox);
iy = iy(iy >= 1 & iy <= Ny_vox);
iz = iz(iz >= 1 & iz <= Nz_vox);

% Abort if any direction produced an empty subdomain.
if isempty(ix) || isempty(iy) || isempty(iz)
    error(['Empty subdomain.\n' ...
        'xvec_vox: [%.6f, %.6f], requested: [%.6f, %.6f]\n' ...
        'yvec_vox: [%.6f, %.6f], requested: [%.6f, %.6f]\n' ...
        'zvec_vox: [%.6f, %.6f], requested: [%.6f, %.6f]'], ...
        xvec_vox(1), xvec_vox(end), x_range(1), x_range(2), ...
        yvec_vox(1), yvec_vox(end), y_range(1), y_range(2), ...
        zvec_vox(1), zvec_vox(end), z_range(1), z_range(2));
end

% Report requested physical size.
fprintf('Requested domain size: %.5f x %.5f x %.5f mm^3\n', ...
    diff(x_range), diff(y_range), diff(z_range));

%% ---------------------- EXTRACT SUBDOMAIN PHASE MAP ---------------------
% Extract voxel left-node coordinate vectors for the selected subdomain.
xvec_sub  = xvec_vox(ix);
yvec_sub  = yvec_vox(iy);
zvec_sub  = zvec_vox(iz);

% Extract the phase map for the selected subdomain.
%
% Remember phase array ordering:
%   phase_sub(y_index, x_index, z_index)
phase_sub = phase3d(iy, ix, iz);

% Subdomain voxel counts in each direction.
Ny_sub = numel(iy);
Nx_sub = numel(ix);
Nz_sub = numel(iz);

fprintf('Core subdomain voxels: Nx=%d x Ny=%d x Nz=%d\n', ...
    Nx_sub, Ny_sub, Nz_sub);

%% ---------------------- OPTIONAL TOP Z SOLID PAD -------------------------
% Build optional fully solid padding layers on the top Z face of the selected
% subdomain.
%
% Purpose:
%   If the top surface of the selected subdomain contains pores, applying
%   boundary conditions directly on that surface can be inconvenient for
%   downstream PRISMS-Plasticity simulations. The optional top pad creates a
%   pore-free solid layer above the selected core subdomain.
%
% Method:
%   For each X-Y column, scan downward from the top of the selected core
%   subdomain until a solid grain voxel is found. Use that grain ID to fill all
%   top-pad voxels in the same X-Y column.
%
% Result:
%   - phase_sub_padded contains the original selected subdomain plus pad_z_top
%     additional solid layers at the top Z face.
%   - No pore voxels are introduced into the added pad region.
%   - If pad_z_top = 0, phase_sub_padded is identical to phase_sub.
%
% Phase convention:
%   Solid grains are labeled 1:N_grains.
%   The global pore phase is labeled N_grains + 1.
pore_phase_val = N_grains + 1;   % global pore phase ID

% For each X-Y column, store the nearest solid grain ID found by scanning
% downward from the top of the selected core subdomain.
%
% Array convention:
%   dimension 1 = Y
%   dimension 2 = X
%   dimension 3 = Z
solid_id_top = zeros(Ny_sub, Nx_sub, 'uint32');

for ky = 1:Ny_sub
    for kx = 1:Nx_sub
        found = false;

        % Search from top to bottom for the nearest solid voxel in this column.        
        for kz = Nz_sub:-1:1
            pid = phase_sub(ky, kx, kz);
            if pid >= 1 && pid <= N_grains
                solid_id_top(ky, kx) = pid;
                found = true;
                break;
            end
        end

        % If no solid voxel exists in this X-Y column, mark it for fallback fill.
        if ~found
            solid_id_top(ky, kx) = 0;
        end
    end
end

% Determine a fallback solid grain ID for columns that are entirely pore.
%
% Preferred fallback:
%   the most common solid grain ID in the top slice.
%
% Secondary fallback:
%   the most common solid grain ID anywhere in the selected core subdomain.
top_slice_solid = phase_sub(:, :, end);
top_slice_solid(top_slice_solid == pore_phase_val) = 0;
solid_vals = top_slice_solid(top_slice_solid > 0);

if ~isempty(solid_vals)
    fallback_id = uint32(mode(double(solid_vals(:))));
else
    core_solid  = phase_sub(phase_sub >= 1 & phase_sub <= N_grains);
    fallback_id = uint32(mode(double(core_solid(:))));
end

% Build the top pad block.
%
% pad_block has dimensions:
%   [Ny_sub, Nx_sub, pad_z_top]
%
% Each added layer repeats the nearest-solid/fallback grain ID for that X-Y
% column.
pad_block        = repmat(solid_id_top, [1, 1, pad_z_top]);
phase_sub_padded = cat(3, phase_sub, pad_block);
Nz_sub_padded    = Nz_sub + pad_z_top;


% Extend the subdomain Z left-node vector for the padded layers.
%
% zvec_sub contains left-node coordinates for the original selected voxels.
% zvec_sub_padded contains left-node coordinates for both the selected core
% voxels and the appended top-pad voxels.
if pad_z_top > 0
    dz_sub = zvec_sub(2) - zvec_sub(1);
    zvec_pad_extra  = zvec_sub(end) + dz_sub*(1:pad_z_top)';
    zvec_sub_padded = [zvec_sub; zvec_pad_extra];
else
    zvec_sub_padded = zvec_sub;
end


%% ---------------------- LOCAL GRAIN-ID REMAPPING -------------------------
% Build a compact local grain-ID system for the exported PRISMS-Plasticity
% subdomain.
%
% Full-domain phase convention:
%   - Solid grains are labeled 1:N_grains.
%   - The global pore phase is labeled N_grains + 1.
%
% Export/subdomain convention:
%   - Only grains present in the selected, padded subdomain are retained.
%   - Present solid grains are remapped to compact local IDs:
%
%         1:N_grains_sub
%
%   - A local pore phase ID is defined as:
%
%         pore_id_sub = N_grains_sub + 1
%
%     but pore elements are skipped when writing the Gmsh mesh.
%
% Why remap?
%   PRISMS-Plasticity orientation files are easier to manage when grain IDs are
%   compact and correspond directly to rows in ori.txt.
%
% Note:
%   The grain list is built from phase_sub_padded so that any optional top-pad
%   grains are included in the orientation file and ID map.
mask_grain_global = (phase_sub_padded >= 1) & (phase_sub_padded <= N_grains);

% Sorted list of global grain IDs present in the exported subdomain.
grain_ids_global  = unique(phase_sub_padded(mask_grain_global));
grain_ids_global  = sort(grain_ids_global(:));

% Number of solid grains present in the exported subdomain.
N_grains_sub = numel(grain_ids_global);

% Local pore ID used internally for bookkeeping.
% Pore elements are not written to the mesh.
pore_id_sub  = N_grains_sub + 1;

fprintf('Grains in domain: %d (pore phase in export = %d)\n', ...
    N_grains_sub, pore_id_sub);

% Create a lookup table from global grain ID to local grain ID.
%
% Example:
%   global2local(global_id) = local_id
%
% Grain IDs not present in the subdomain remain zero.
global2local = zeros(N_grains, 1, 'uint32');
global2local(grain_ids_global) = uint32(1:N_grains_sub);

% Convert the padded subdomain phase map from global grain IDs to compact
% local grain IDs.
%
% All entries are initialized as local pore ID. Solid grain entries are then
% overwritten with their corresponding local grain IDs.
phase_sub_local = uint32(pore_id_sub) * ones(size(phase_sub_padded), 'uint32');
phase_sub_local(mask_grain_global) = global2local(phase_sub_padded(mask_grain_global));


% ---------------------- WRITE GRAIN ID MAP -------------------------------
% Write a traceability file mapping local grain IDs used in the exported mesh
% and orientation file back to the original full-domain global grain IDs.
%
% Columns:
%   local_id  global_id
fid_map = fopen('grain_id_map_subdomain.txt', 'w');
fprintf(fid_map, 'local_id\tglobal_id\n');
fprintf(fid_map, '%d\t%d\n', [(1:N_grains_sub).', grain_ids_global].');
fclose(fid_map);

%% ---------------------- BUILD MESH NODE COORDINATES ----------------------
% Construct the node coordinate vectors for the exported Gmsh mesh.
%
% The selected subdomain is voxel-based. Each voxel becomes one hexahedral
% element, so the mesh requires one more node than the number of voxels in each
% direction:
%
%   Nx_nodes = Nx_sub + 1
%   Ny_nodes = Ny_sub + 1
%   Nz_nodes = Nz_sub_padded + 1
%
% xvec_sub, yvec_sub, and zvec_sub contain the left-node coordinates of the
% selected voxels. The final node in each direction is appended by adding one
% grid spacing.
%
% Units:
%   Coordinates remain in millimeters.
%
% Output convention:
%   The exported mesh is shifted so that the minimum corner of the selected
%   subdomain is at the origin.

% Build physical node coordinates in X and Y.
% Append the final right/top node by adding one grid spacing.
xnodes_phys = [xvec_sub;       xvec_sub(end)       + dx];
ynodes_phys = [yvec_sub;       yvec_sub(end)       + dy];

% Build physical node coordinates in Z for the original selected core.
% The optional Z-pad extension is handled below.
znodes_phys = [zvec_sub;       zvec_sub(end)       + dz];

% If top Z padding was added, append additional Z nodes to close the padded
% hexahedral elements.
%
% For pad_z_top layers, pad_z_top additional voxels require pad_z_top
% additional top nodes beyond the original top node.
if pad_z_top > 0
    dz_sub = zvec_sub(2) - zvec_sub(1);
    znodes_phys = [znodes_phys; ...
        zvec_sub(end) + dz_sub*(1:pad_z_top)' + dz_sub];
end

% Shift origin to the minimum corner of the selected subdomain.
% This keeps PRISMS-Plasticity input coordinates local to the exported mesh.
xnodes = xnodes_phys - xnodes_phys(1);
ynodes = ynodes_phys - ynodes_phys(1);
znodes = znodes_phys - znodes_phys(1);

% Mesh node counts.
Nx_nodes  = length(xnodes);
Ny_nodes  = length(ynodes);
Nz_nodes  = length(znodes);
num_nodes = Nx_nodes * Ny_nodes * Nz_nodes;

% Linear node-ID mapping used for Gmsh output.
%
% Node ordering:
%   kx varies fastest, then ky, then kz.
%
% This ordering must be consistent with the node-writing loop and element
% connectivity below.
nodeid = @(kx,ky,kz) (kz-1)*Ny_nodes*Nx_nodes + (ky-1)*Nx_nodes + kx;

%% ---------------------- WRITE GMSH NODES -------------------------------
% Write the node block of the Gmsh 2.2 ASCII mesh file.
%
% Gmsh node format:
%   node_id  x  y  z
%
% Coordinates are written in millimeters and have already been shifted so that
% the exported subdomain starts at the origin.
%
% Node ordering:
%   kx varies fastest, then ky, then kz.
%   This ordering must match the nodeid(kx,ky,kz) mapping used for element
%   connectivity.
fid = fopen(msh_filename, 'w');

% Gmsh 2.2 ASCII mesh header.
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');

% Write node count.
fprintf(fid, '$Nodes\n%d\n', num_nodes);

% Write nodal coordinates.
nid = 1;
for kz = 1:Nz_nodes
    for ky = 1:Ny_nodes
        for kx = 1:Nx_nodes
            fprintf(fid, '%d %.10g %.10g %.10g\n', ...
                nid, xnodes(kx), ynodes(ky), znodes(kz));
            nid = nid + 1;
        end
    end
end

fprintf(fid, '$EndNodes\n');
fclose(fid);

%% ---------------------- WRITE GMSH ELEMENTS -----------------------------
% Write the element block of the Gmsh 2.2 ASCII mesh file.
%
% Each non-pore voxel in phase_sub_local becomes one 8-node hexahedral element.
% Pore voxels are skipped and therefore do not appear in the exported mesh.
%
% Gmsh element format used here:
%
%   elem_id  5  2  physical_tag  elementary_tag  n1 n2 n3 n4 n5 n6 n7 n8
%
% where:
%   elem_id
%     Sequential element ID in the exported mesh.
%
%   5
%     Gmsh element type for an 8-node hexahedron.
%
%   2
%     Number of tags.
%
%   physical_tag, elementary_tag
%     Both are set to the compact local grain ID. This allows grain/material
%     identity to be recovered from the mesh element tags.
%
%   n1...n8
%     Node IDs of the hexahedral element.
%
% Element connectivity convention:
%   n1 = (kx,   ky,   kz)
%   n2 = (kx+1, ky,   kz)
%   n3 = (kx+1, ky+1, kz)
%   n4 = (kx,   ky+1, kz)
%   n5 = (kx,   ky,   kz+1)
%   n6 = (kx+1, ky,   kz+1)
%   n7 = (kx+1, ky+1, kz+1)
%   n8 = (kx,   ky+1, kz+1)
%
% This connectivity must be compatible with the element orientation expected by
% PRISMS-Plasticity/Gmsh.
fid = fopen(msh_filename, 'a');

% Count and store element lines before writing the $Elements block, because
% Gmsh requires the element count before the element records.
elements_written = 0;
eid              = 1;
element_lines    = cell(Nx_sub * Ny_sub * Nz_sub_padded, 1);

for kz = 1:Nz_sub_padded
    for ky = 1:Ny_sub
        for kx = 1:Nx_sub

            % Local grain ID for this voxel.
            gID = phase_sub_local(ky, kx, kz);

            % Skip pore voxels. Only solid elements are exported.            
            if gID == pore_id_sub
                continue;
            end

            % Hexahedral element node connectivity.
            n1 = nodeid(kx,   ky,   kz  );
            n2 = nodeid(kx+1, ky,   kz  );
            n3 = nodeid(kx+1, ky+1, kz  );
            n4 = nodeid(kx,   ky+1, kz  );
            n5 = nodeid(kx,   ky,   kz+1);
            n6 = nodeid(kx+1, ky,   kz+1);
            n7 = nodeid(kx+1, ky+1, kz+1);
            n8 = nodeid(kx,   ky+1, kz+1);

            % Store the Gmsh element record.
            % The two tags are both set to gID.
            element_lines{eid} = sprintf( ...
                '%d 5 2 %d %d %d %d %d %d %d %d %d %d\n', ...
                elements_written+1, gID, gID, ...
                n1, n2, n3, n4, n5, n6, n7, n8);

            elements_written = elements_written + 1;
            eid              = eid + 1;
        end
    end
end

% Remove unused preallocated cells.
element_lines = element_lines(1:elements_written);

% Write the Gmsh element block.
fprintf(fid, '$Elements\n%d\n', elements_written);
fprintf(fid, '%s', element_lines{:});
fprintf(fid, '$EndElements\n');
fclose(fid);


%% ---------------------- WRITE ORIENTATION FILE ---------------------------
% Write the PRISMS-Plasticity orientation file for the exported subdomain.
%
% The mesh uses compact local grain IDs rather than the original full-domain
% global grain IDs. Therefore, the orientation file must contain orientations
% in the same local grain-ID convention used by the mesh element tags.
%
% Input orientation convention:
%   ori_vec stores full-domain Rodrigues vectors as:
%
%     [r1_x; r1_y; r1_z; r2_x; r2_y; r2_z; ...]
%
% Output orientation convention:
%   ori.txt contains one row per local grain:
%
%     local_grain_id   r_x   r_y   r_z
%
% where local_grain_id corresponds to the compact ID used in micro.msh.
fprintf('Writing orientation file: %s\n', ori_filename);

% Convert the full orientation vector into an N_grains-by-3 matrix.
%
% Row g contains the Rodrigues vector for global grain ID g.
ori_all = reshape(ori_vec, 3, []).';

% Verify that the orientation array contains at least one Rodrigues vector for
% every full-domain grain.
if size(ori_all, 1) < N_grains
    error('ori_vec size mismatch: expected >= 3*N_grains entries.');
end

% Extract orientations for only the global grain IDs present in the exported
% subdomain.
%
% grain_ids_global(local_id) gives the original global grain ID corresponding
% to each compact local grain ID.
ori_sub = ori_all(grain_ids_global, :);

% Prepend compact local grain IDs.
%
% ori_mat columns:
%   1 : local grain ID
%   2 : Rodrigues r_x
%   3 : Rodrigues r_y
%   4 : Rodrigues r_z
ori_mat = [(1:N_grains_sub).', ori_sub];

% Write orientation file.
fid = fopen(ori_filename, 'w');
fprintf(fid, '**Grain ID, Rodrigues vector r_x, r_y, r_z\n');
fprintf(fid, '%d\t%.15g\t%.15g\t%.15g\n', ori_mat');
fclose(fid);

fprintf('Orientations written for %d grains.\n', N_grains_sub);

% Final status message.
disp('Finished exporting mesh and orientation.');
disp(['-> Mesh: ', msh_filename, ', Ori: ', ori_filename]);
