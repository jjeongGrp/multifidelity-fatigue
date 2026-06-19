%% Columnar/Equiaxed Microstructure Generation with XCT Pores
% -------------------------------------------------------------------------
% This script generates a 3D voxel-based synthetic polycrystalline
% microstructure inside an XCT-derived cubic domain.
%
% The workflow combines:
%   - Lognormally distributed grain diameters in the transverse XY plane
%   - Lognormally distributed grain aspect ratios in the Z direction
%   - Randomly seeded anisotropic grains
%   - Voxel-wise nearest-grain assignment using an ellipsoidal distance metric
%   - XCT-derived pore masking from a preprocessed cube-domain file
%   - Optional fiber-textured or random crystallographic orientation assignment
%   - IPF-based visualization of grain orientations and phase maps
%
% The grain structure is generated on the same voxel grid as the XCT cube.
% Pore voxels from the XCT-derived binary mask overwrite the grain phase map
% after grain assignment. The resulting phase map and orientations are exported 
% to HDF5 for downstream finite element or microstructure-based analyses.
%
% Inputs:
%   - Cube3p003mm_546vox_547grid_PoreDomain.mat
%       MAT-file containing the structure `cube_domain`, including:
%         * cube_domain.dx_mm       : voxel size in mm
%         * cube_domain.cube_mm     : cubic domain side length in mm
%         * cube_domain.x_grid_mm   : X grid coordinates
%         * cube_domain.y_grid_mm   : Y grid coordinates
%         * cube_domain.z_grid_mm   : Z grid coordinates
%         * cube_domain.pore_mask   : logical voxel mask for XCT pores
%
%   - input_structure_poly.h5
%       HDF5 file containing:
%         * /pix         : nodal phase/region IDs flattened with X varying fastest
%         * /orientation : grain orientations stored as Rodrigues vectors
%       The /pix dataset uses the padded domain converted from voxel phases
%       to nodal phases.
%
% Main Outputs:
%   - phase
%       3D voxel-wise phase map for the original cube.
%       Grain IDs are assigned as 1:N_grains.
%       Pore voxels are assigned pore_phase_id = N_grains + 1.
%
%   - phase2
%       Padded phase map with an XY air/pore ring and solid Z caps.
%
%   - ori_vec
%       Grain orientation array stored as Rodrigues vectors.
%       The vector is arranged as:
%           [r1_x; r1_y; r1_z; r2_x; r2_y; r2_z; ...]
%
%   - voxel_volume_summary.txt
%       Text summary of padded-domain dimensions, voxel counts, solid voxels,
%       air/pore voxels, and key phase IDs.
%
% Visualizations:
%   - Pole figure of sampled grain orientations
%   - XY, YZ, and XZ slice maps colored by IPF orientation
%   - 3D surface scatter plot of the voxelized microstructure
%
% Notes:
%   - The same phase ID is used for XCT pores and surrounding air/padding:
%       pore_phase_id = N_grains + 1
%   - Grain assignment is accelerated by processing the domain in blocks and
%     filtering the candidate grains considered in each block.
%   - Crystallographic orientation calculations require MTEX.
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
% License:       MIT License (See LICENSE file in the repository)
% -------------------------------------------------------------------------
clear; close all;


%% 1. Microstructure Parameters

% --- Grain size statistics ---
% Mean and standard deviation of the grain diameter in the XY plane.
% Units are millimeters.
mean_d = 0.0525;
std_d  = 0.18*mean_d;

% --- Grain aspect-ratio statistics ---
% Aspect ratio is defined as grain length in Z divided by grain width in XY.
% Individual grain aspect ratios are sampled from a lognormal distribution.
mean_aspect = 3.52;
std_aspect = 0.05 * mean_aspect;

% --- Orientation settings ---
% orientation_type controls the crystallographic texture:
%   'textured' : fiber texture with controlled angular spread
%   'random'   : uniformly random orientations
orientation_type = 'textured';
sigma_spread = 20;             % Angular spread in degrees for textured case



%% 2. Load XCT Cube Domain and Define Grid

% Load the preprocessed XCT cube domain.
% The file should contain a structure named `cube_domain`.
file_name_voxel = 'Cube3p003mm_546vox_547grid_PoreDomain.mat';
cube_file = file_name_voxel;
load(cube_file, 'cube_domain');

% Number of grain orientations to sample for pole-figure visualization.
K_grain = 2000;

% --- Domain and grid information ---
% dx is the voxel spacing in millimeters.
dx       = cube_domain.dx_mm;
dx_um    = dx * 1e3;

% cube_len is the side length of the cubic domain in millimeters.
cube_len    = cube_domain.cube_mm;
cube_len_um = cube_len * 1e3;

% Grid-point coordinate vectors.
% These describe the node/grid coordinates of the cube.
xvec = cube_domain.x_grid_mm(:);
yvec = cube_domain.y_grid_mm(:);
zvec = cube_domain.z_grid_mm(:);

Nx = numel(xvec);
Ny = numel(yvec);
Nz = numel(zvec);

% XCT-derived pore mask on the voxel/cell grid.
% True values indicate pore voxels.
pore_mask = cube_domain.pore_mask;

% Voxel dimensions of the pore/phase grid.
Nyv = size(pore_mask, 1);
Nxv = size(pore_mask, 2);
Nzv = size(pore_mask, 3);

% Total number of voxels and physical domain volume.
num_voxels = Nyv * Nxv * Nzv;
domain_vol  = cube_len^3;

% Voxel-center coordinate vectors.
% These coordinates are used for assigning each voxel to a grain.
x_cent = ((1:Nxv) - 0.5) * dx;
y_cent = ((1:Nyv) - 0.5) * dx;
z_cent = ((1:Nzv) - 0.5) * dx;

% 3D coordinate arrays for voxel centers.
[Yc, Xc, Zc] = ndgrid(y_cent, x_cent, z_cent);

%% 3. Generate Lognormal Grain Sizes, Aspect Ratios, and Seeds

% Estimate the number of grains required to fill the domain using an
% approximate ellipsoidal/cylindrical grain volume based on the mean
% transverse diameter and mean aspect ratio.
mean_grain_vol = pi*(mean_d/2)^2 * (mean_aspect*mean_d);
N_grains = ceil(domain_vol / mean_grain_vol);

% Sample grain diameters in the XY plane from a lognormal distribution.
sigma_log = sqrt(log((std_d/mean_d)^2 + 1));
mu_log = log(mean_d) - 0.5*sigma_log^2;
grain_diameters = lognrnd(mu_log, sigma_log, N_grains, 1);

% Transverse grain radius in XY.
r_xy = grain_diameters / 2;

% Sample grain aspect ratios from a lognormal distribution.
sigma_log_a = sqrt(log((std_aspect/mean_aspect)^2 + 1));
mu_log_a    = log(mean_aspect) - 0.5*sigma_log_a^2;
aspect_ratios = lognrnd(mu_log_a, sigma_log_a, N_grains, 1);

% Longitudinal grain radius in Z.
r_z  = aspect_ratios .* r_xy;

% Randomly seed grain centers throughout the cubic domain.
grain_seeds = [rand(N_grains,1)*cube_len, ...
    rand(N_grains,1)*cube_len, ...
    rand(N_grains,1)*cube_len];



%% 4. Anisotropic Grain Assignment Using Parallel Block Processing

% Assign each voxel to the nearest grain according to an anisotropic
% ellipsoidal distance metric:
%
%   d^2 = (dx/r_xy)^2 + (dy/r_xy)^2 + (dz/r_z)^2
%
% This produces elongated columnar/equiaxed grains depending on the sampled
% aspect ratios. The computation is performed in spatial blocks to reduce
% memory usage and parallelized using MATLAB's Parallel Computing Toolbox.

fprintf('\n=== PARALLEL + SPATIAL FILTERING MICROSTRUCTURE GENERATION ===\n');
main_timer = tic;

% Start a parallel pool if one is not already active.
if isempty(gcp('nocreate'))
    fprintf('Starting parallel pool...\n');
    pool_timer = tic;
    num_workers = feature('numcores');
    parpool('local', num_workers);
    fprintf('Parallel pool started in %.1fs with %d workers\n', toc(pool_timer), gcp().NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', gcp().NumWorkers);
end

% Preallocate the voxel-wise phase array.
% phase = 0 initially indicates unassigned voxels.
phase = zeros(Nyv, Nxv, Nzv, 'uint32');

% Precompute inverse squared radii for faster distance calculations.
inv_r_xy2 = single(1 ./ (r_xy.^2));
inv_r_z2 = single(1 ./ (r_z.^2));

% Store grain seeds in single precision to reduce memory use.
grain_seeds_single = single(grain_seeds);

% Build a coarse spatial-grid diagnostic structure.
% The current block assignment below uses direct bounding-box filtering
% rather than looking up grains from this structure.
fprintf('\nBuilding spatial acceleration structure...\n');
accel_timer = tic;

grid_res = 15;
x_edges = linspace(min(xvec), max(xvec), grid_res+1);
y_edges = linspace(min(yvec), max(yvec), grid_res+1);
z_edges = linspace(min(zvec), max(zvec), grid_res+1);

% Maximum approximate grain influence distance used for candidate filtering.
max_influence = max([r_xy(:); r_z(:)]) * 1.5;
total_cells = grid_res^3;

% Temporary flat cell array for parfor-compatible spatial-grid construction.
spatial_grid_temp = cell(total_cells, 1);

% Parallel construction using linear indexing (parfor-compatible)
parfor linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);

    cell_xmin = x_edges(i); cell_xmax = x_edges(i+1);
    cell_ymin = y_edges(j); cell_ymax = y_edges(j+1);
    cell_zmin = z_edges(k); cell_zmax = z_edges(k+1);

    % Find all grains whose centers are close enough to potentially affect this cell
    relevant_grains = find(...
        grain_seeds_single(:,1) >= (cell_xmin - max_influence) & ...
        grain_seeds_single(:,1) <= (cell_xmax + max_influence) & ...
        grain_seeds_single(:,2) >= (cell_ymin - max_influence) & ...
        grain_seeds_single(:,2) <= (cell_ymax + max_influence) & ...
        grain_seeds_single(:,3) >= (cell_zmin - max_influence) & ...
        grain_seeds_single(:,3) <= (cell_zmax + max_influence));

    spatial_grid_temp{linear_idx} = relevant_grains;
end

% Convert the flat parfor output into a 3D cell array indexed by spatial bin.
spatial_grid = cell(grid_res, grid_res, grid_res);
for linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);
    spatial_grid{i,j,k} = spatial_grid_temp{linear_idx};
end
clear spatial_grid_temp;

fprintf('Spatial acceleration structure built in %.1fs\n', toc(accel_timer));

% --- Block-processing setup ---
% The voxel grid is divided into smaller blocks to reduce peak memory usage.
block_size = 50;

n_blocks_x = ceil(Nxv / block_size);
n_blocks_y = ceil(Nyv / block_size);
n_blocks_z = ceil(Nzv / block_size);
total_blocks = n_blocks_x * n_blocks_y * n_blocks_z;

% Report block setup
fprintf('Processing %d blocks (%dx%dx%d) in parallel:\n', ...
    total_blocks, n_blocks_x, n_blocks_y, n_blocks_z);

% Group blocks into batches to allow progress reporting during parfor work.
num_workers = gcp().NumWorkers;
blocks_per_batch = max(1, floor(total_blocks / (num_workers * 4)));
n_batches = ceil(total_blocks / blocks_per_batch);

fprintf('Using %d batches of %d blocks each for progress tracking\n', ...
    n_batches, blocks_per_batch);

% Track progress and performance statistics.
completed_blocks = 0;
total_grains_processed = 0;
blocks_with_grains = 0;
batch_times = [];

% Process the domain one batch at a time.
% Within each batch, blocks are processed in parallel.
fprintf('\nStarting parallel processing with real-time progress:\n');
process_timer = tic;

for batch_idx = 1:n_batches
    batch_timer = tic;

    % Determine blocks in this batch
    batch_start = (batch_idx - 1) * blocks_per_batch + 1;
    batch_end = min(batch_idx * blocks_per_batch, total_blocks);
    current_batch_size = batch_end - batch_start + 1;

    % Create list of block indices for current batch
    block_list_batch = zeros(current_batch_size, 4);
    for local_idx = 1:current_batch_size
        block_num = batch_start + local_idx - 1;
        bz = ceil(block_num / (n_blocks_x * n_blocks_y));
        remaining = block_num - (bz-1) * n_blocks_x * n_blocks_y;
        by = ceil(remaining / n_blocks_x);
        bx = remaining - (by-1) * n_blocks_x;
        block_list_batch(local_idx, :) = [bx, by, bz, block_num];
    end

    % Pre-allocate cell arrays to collect block results
    batch_results = cell(current_batch_size, 1);
    batch_indices = cell(current_batch_size, 1);
    batch_grain_counts = zeros(current_batch_size, 1);

    % Extract grid parameters for parfor
    xvec_local = xvec;
    yvec_local = yvec;
    zvec_local = zvec;
    x_edges_local = x_edges;
    y_edges_local = y_edges;
    z_edges_local = z_edges;

    % Process current batch in parallel
    parfor local_idx = 1:current_batch_size
        bx = block_list_batch(local_idx, 1);
        by = block_list_batch(local_idx, 2);
        bz = block_list_batch(local_idx, 3);

        % Determine voxel-index bounds for the current block.
        ix1 = (bx-1)*block_size + 1;
        ix2 = min(bx*block_size, Nxv);
        iy1 = (by-1)*block_size + 1;
        iy2 = min(by*block_size, Nyv);
        iz1 = (bz-1)*block_size + 1;
        iz2 = min(bz*block_size, Nzv);

        % Determine physical coordinate bounds of the current block.
        block_xmin = xvec_local(ix1); block_xmax = xvec_local(ix2);
        block_ymin = yvec_local(iy1); block_ymax = yvec_local(iy2);
        block_zmin = zvec_local(iz1); block_zmax = zvec_local(iz2);

        % Compute the spatial-grid cell range overlapped by this block.
        % These indices are currently not used for candidate lookup because the
        % candidate grains are filtered directly from seed positions below.
        gx1 = max(1, floor((block_xmin - min(xvec_local)) / (max(xvec_local) - min(xvec_local)) * grid_res) + 1);
        gx2 = min(grid_res, ceil((block_xmax - min(xvec_local)) / (max(xvec_local) - min(xvec_local)) * grid_res));
        gy1 = max(1, floor((block_ymin - min(yvec_local)) / (max(yvec_local) - min(yvec_local)) * grid_res) + 1);
        gy2 = min(grid_res, ceil((block_ymax - min(yvec_local)) / (max(yvec_local) - min(yvec_local)) * grid_res));
        gz1 = max(1, floor((block_zmin - min(zvec_local)) / (max(zvec_local) - min(zvec_local)) * grid_res) + 1);
        gz2 = min(grid_res, ceil((block_zmax - min(zvec_local)) / (max(zvec_local) - min(zvec_local)) * grid_res));

        % Identify candidate grains whose seed locations are close enough to
        % potentially influence this block. This avoids testing every grain against
        % every voxel.
        relevant_grains = find(...
            grain_seeds_single(:,1) >= (block_xmin - max_influence) & ...
            grain_seeds_single(:,1) <= (block_xmax + max_influence) & ...
            grain_seeds_single(:,2) >= (block_ymin - max_influence) & ...
            grain_seeds_single(:,2) <= (block_ymax + max_influence) & ...
            grain_seeds_single(:,3) >= (block_zmin - max_influence) & ...
            grain_seeds_single(:,3) <= (block_zmax + max_influence));

        % If candidate grains are present, assign each voxel in the block to the
        % grain with the minimum normalized ellipsoidal distance.
        if ~isempty(relevant_grains)
            % Create voxel-center coordinate arrays for the current block.
            x_cent_local = x_cent;
            y_cent_local = y_cent;
            z_cent_local = z_cent;
            [Yb, Xb, Zb] = ndgrid(single(y_cent_local(iy1:iy2)), ...
                single(x_cent_local(ix1:ix2)), ...
                single(z_cent_local(iz1:iz2)));

            % Initialize local block phase map and minimum-distance array.
            block_phase = zeros(size(Yb), 'uint32');
            min_dist2   = inf(size(Yb), 'single');

            % Loop only over grains that can potentially influence this block.
            for k = relevant_grains'
                dX = Xb - grain_seeds_single(k,1);
                dY = Yb - grain_seeds_single(k,2);
                dZ = Zb - grain_seeds_single(k,3);

                % Compute anisotropic squared distance from each voxel to
                % the current grain seed.
                current_dist2 = (dX.^2) * inv_r_xy2(k) + (dY.^2) * inv_r_xy2(k) + (dZ.^2) * inv_r_z2(k);

                % Update voxel assignments where this grain is closer than
                % the previously stored nearest grain.
                update_mask = current_dist2 < min_dist2;
                block_phase(update_mask) = k;
                min_dist2(update_mask) = current_dist2(update_mask);
            end

            batch_results{local_idx} = block_phase;
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = length(relevant_grains);
        else
            batch_results{local_idx} = [];
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = 0;
        end
    end

    % Store batch results in main phase array
    batch_grains_processed = 0;
    batch_blocks_with_grains = 0;

    for local_idx = 1:current_batch_size
        if ~isempty(batch_results{local_idx})
            indices = batch_indices{local_idx};
            ix1 = indices(1); ix2 = indices(2);
            iy1 = indices(3); iy2 = indices(4);
            iz1 = indices(5); iz2 = indices(6);

            phase(iy1:iy2, ix1:ix2, iz1:iz2) = cast(batch_results{local_idx}, 'uint32');
            batch_blocks_with_grains = batch_blocks_with_grains + 1;
        end
        % Sum up grains processed in batch
        batch_grains_processed = batch_grains_processed + batch_grain_counts(local_idx);
    end

    % Update global counters
    completed_blocks = completed_blocks + current_batch_size;
    total_grains_processed = total_grains_processed + batch_grains_processed;
    blocks_with_grains = blocks_with_grains + batch_blocks_with_grains;

    % Record batch timing
    batch_time = toc(batch_timer);
    batch_times(end+1) = batch_time;

    % Calculate progress statistics
    progress_pct = 100 * completed_blocks / total_blocks;
    elapsed_total = toc(main_timer);
    avg_grains_per_block = batch_grains_processed / current_batch_size;

    % Estimate ETA for remaining batches
    if batch_idx > 1
        avg_batch_time = mean(batch_times);
        remaining_batches = n_batches - batch_idx;
        eta_seconds = remaining_batches * avg_batch_time;
    else
        eta_seconds = batch_time * (n_batches - 1);
    end


    % Print detailed progress
    fprintf('Batch %d/%d: Blocks %d-%d (%.1f%%) | %.1f grains/block avg | Elapsed: %.1fs | ETA: %.1fs (%.1f min)\n', ...
        batch_idx, n_batches, batch_start, batch_end, progress_pct, avg_grains_per_block, elapsed_total, eta_seconds, eta_seconds/60);

    if mod(batch_idx, max(1, floor(n_batches/5))) == 0 || batch_idx == n_batches
        overall_avg_grains = total_grains_processed / completed_blocks;
        processing_rate = completed_blocks / elapsed_total;
        speedup_factor = N_grains / overall_avg_grains;

        fprintf('  --> Detailed Stats: %.1f grains/block overall | %.1f blocks/sec | %.1fx speedup from spatial filtering\n', ...
            overall_avg_grains, processing_rate, speedup_factor);

        if batch_idx < n_batches
            fprintf('  --> Memory: %d/%d blocks processed | %d/%d non-empty | %.1f%% efficiency\n', ...
                completed_blocks, total_blocks, blocks_with_grains, completed_blocks, 100*blocks_with_grains/completed_blocks);
        end
    end
end

% Report final runtime and assignment completeness.
total_time = toc(main_timer);
avg_grains_per_block = total_grains_processed / completed_blocks;
speedup_factor = N_grains / avg_grains_per_block;
processing_rate = double(Nyv*Nxv*Nzv) / total_time;

fprintf('\n=== COMPLETION SUMMARY ===\n');
fprintf('Total time: %.1f seconds (%.2f minutes)\n', total_time, total_time/60);
fprintf('Processing rate: %s voxels/second\n', regexprep(sprintf('%.0f', processing_rate),'(\d)(?=(\d{3})+$)','$1,'));

% Check whether any voxels remained unassigned after grain partitioning.
unassigned = nnz(phase(:) == 0);
if unassigned > 0
    fprintf('\nWARNING: %d unassigned voxels detected!\n', unassigned);
else
    fprintf('SUCCESS: All %s voxels assigned successfully!\n', ...
        regexprep(sprintf('%.0f', double(Nyv)*double(Nxv)*double(Nzv)), '(\d)(?=(\d{3})+$)', '$1,'));
end

fprintf('\nParallel + Spatial filtering microstructure generation completed!\n');



%% 5. Apply XCT Pore Mask

% Overwrite grain assignments with the pore phase wherever the XCT-derived
% pore mask is true. Pore/void voxels are assigned a phase ID of N_grains + 1.
phase(pore_mask) = N_grains + 1;
fprintf('Phase assignment completed.\n');


%% 6. Add XY Air Padding and Solid Z Caps

% Add an XY padding ring and Z-direction end caps around the original volume.
%
% The XY padding is initialized as pore/air phase. The Z caps are generated
% by extruding the solid footprint from the bottom and top gauge-face slices.
%
% Note:
%   The same phase ID, pore_phase_id, is used for both XCT pores and the
%   surrounding air/padding region.
pad_xy = round((Nx-1)*0.02);
pad_z  = round((Nx-1)*0.02);

pore_phase_id = cast(N_grains + 1, 'uint32');

% Original phase-map dimensions.
[Nyv, Nxv, Nzv] = size(phase);

% Dimensions of the padded volume.
Nxv2 = Nxv + 2*pad_xy;
Nyv2 = Nyv + 2*pad_xy;
Nzv2 = Nzv + 2*pad_z;

% Initialize padded domain as pore/air.
phase2 = repmat(pore_phase_id, [Nyv2, Nxv2, Nzv2]);

% Insert the original microstructure into the center of the padded domain.
xs = (pad_xy + 1):(pad_xy + Nxv);
ys = (pad_xy + 1):(pad_xy + Nyv);
zs = (pad_z  + 1):(pad_z  + Nzv);
phase2(ys, xs, zs) = phase;

% Extract the bottom and top gauge-face slices used to define the cap shape.
ref_bot = phase2(ys, xs, pad_z+1);
ref_top = phase2(ys, xs, Nzv2-pad_z);

% Identify the solid footprint on each gauge-face slice.
foot_bot = (ref_bot ~= pore_phase_id);
foot_top = (ref_top ~= pore_phase_id);

% Confirm that each cap footprint is a single connected solid region.
cc_bot = bwconncomp(foot_bot);
cc_top = bwconncomp(foot_top);
if cc_bot.NumObjects > 1
    warning('Bottom cap footprint has %d disconnected regions.', cc_bot.NumObjects);
end
if cc_top.NumObjects > 1
    warning('Top cap footprint has %d disconnected regions.', cc_top.NumObjects);
end

% Determine representative grain IDs from the bottom and top solid slices.
% These IDs are reported in the summary file.
ref_bot_solid = ref_bot(foot_bot);
ref_top_solid = ref_top(foot_top);

if isempty(ref_bot_solid)
    error('Bottom gauge slice has no solid voxels in core window.');
end
if isempty(ref_top_solid)
    error('Top gauge slice has no solid voxels in core window.');
end

fill_id_bot = cast(mode(double(ref_bot_solid(:))), 'uint32');
fill_id_top = cast(mode(double(ref_top_solid(:))), 'uint32');

% Initialize the cap regions as air/pore phase before extrusion.
phase2(ys, xs, 1:pad_z)            = pore_phase_id;
phase2(ys, xs, Nzv2-pad_z+1:Nzv2) = pore_phase_id;

% Extrude the bottom gauge-face footprint through the lower Z cap.
for k = 1:pad_z
    cap_slice = pore_phase_id * ones(size(ref_bot), 'like', ref_bot);
    cap_slice(foot_bot) = ref_bot(foot_bot);
    phase2(ys, xs, k) = cap_slice;
end

% Extrude the top gauge-face footprint through the upper Z cap.
for k = (Nzv2 - pad_z + 1):Nzv2
    cap_slice = pore_phase_id * ones(size(ref_top), 'like', ref_top);
    cap_slice(foot_top) = ref_top(foot_top);
    phase2(ys, xs, k) = cap_slice;
end


%% 7. Write Voxel-Volume Summary

% Write a text summary describing the original and padded domain sizes,
% padding thicknesses, solid/air voxel counts, and key phase identifiers.
summary_file = fullfile(pwd, 'voxel_volume_summary.txt');
fid = fopen(summary_file, 'w');
if fid == -1
    error('Could not open summary file for writing: %s', summary_file);
end

% Compute voxel counts for the original domain, padded domain, and phase types.
total_voxels = numel(phase2);
interior_vox = Nyv * Nxv * Nzv;

% Approximate number of voxels associated with the XY padding ring.
pad_xy_vox = (Nyv2 * Nxv2 - Nyv * Nxv) * Nzv2;

% Number of voxels associated with both Z cap layers.
pad_z_vox    = Nyv2 * Nxv2 * 2 * pad_z;

% Count air/pore and solid voxels in the padded phase map.
air_vox      = nnz(phase2 == pore_phase_id);
solid_vox    = nnz(phase2 ~= pore_phase_id);

fprintf(fid, '============================================================\n');
fprintf(fid, '               Voxel Volume Summary                        \n');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Generated : %s\n', datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
fprintf(fid, '------------------------------------------------------------\n');
fprintf(fid, 'Padding Parameters\n');
fprintf(fid, '  pad_xy (XY air ring)  : %d voxels\n', pad_xy);
fprintf(fid, '  pad_z  (Z cap layers) : %d voxels per cap\n', pad_z);
fprintf(fid, '------------------------------------------------------------\n');
fprintf(fid, 'Domain Dimensions\n');
fprintf(fid, '  Original interior     : %d x %d x %d\n', Nyv, Nxv, Nzv);
fprintf(fid, '  Padded domain         : %d x %d x %d\n', Nyv2, Nxv2, Nzv2);
fprintf(fid, '------------------------------------------------------------\n');
fprintf(fid, 'Voxel Counts\n');
fprintf(fid, '  Total (padded domain) : %d  (100%%)\n', total_voxels);
fprintf(fid, '  Interior (original)   : %d  (%.1f%%)\n', interior_vox, 100*interior_vox/total_voxels);
fprintf(fid, '  XY air padding ring   : %d  (%.1f%%)\n', pad_xy_vox,   100*pad_xy_vox/total_voxels);
fprintf(fid, '  Z cap layers (both)   : %d  (%.1f%%)\n', pad_z_vox,    100*pad_z_vox/total_voxels);
fprintf(fid, '  Solid voxels          : %d  (%.1f%%)\n', solid_vox,    100*solid_vox/total_voxels);
fprintf(fid, '  Air voxels            : %d  (%.1f%%)\n', air_vox,      100*air_vox/total_voxels);
fprintf(fid, '------------------------------------------------------------\n');
fprintf(fid, 'Phase Info\n');
fprintf(fid, '  N_grains              : %d\n', N_grains);
fprintf(fid, '  pore_phase_id         : %d\n', pore_phase_id);
fprintf(fid, '  fill_id_bot           : %d\n', fill_id_bot);
fprintf(fid, '  fill_id_top           : %d\n', fill_id_top);
fprintf(fid, '============================================================\n');

fclose(fid);
fprintf('Volume summary written to: %s\n', summary_file);


%% 8. Assign Crystallographic Orientations

% Store one Rodrigues vector per grain.
% The output vector is arranged as three consecutive components per grain.
ori_vec = zeros(3*N_grains, 1);

if strcmp(orientation_type, 'textured')
    % Define cubic crystal symmetry for the material.
    cs = crystalSymmetry('m-3m');

    % Define the target fiber direction in the crystal frame.
    fiber_axis = normalize(vector3d(1,0,1));

    % Select a reference direction that is not parallel to the fiber axis.
    % This helps define a stable base rotation.
    ref = vector3d(0,1,0);
    if angle(fiber_axis, ref) < 1*degree
        ref = vector3d(1,0,0);
    end

    % Define a crystal-frame direction perpendicular to the fiber axis and
    % a corresponding sample-frame reference direction.
    c_perp = normalize(cross(fiber_axis, ref));     % perpendicular to fiber_axis
    s_perp = yvector;                               % define sample reference

    % Construct a base rotation that maps the crystal fiber axis to sample Z.
    % The secondary direction fixes the rotation about the fiber axis.
    R_base = rotation.map(fiber_axis, zvector, c_perp, s_perp);

    % Reference crystal directions used for rejection filtering.
    h111 = normalize(vector3d(1,1,1));

    % Reject orientations whose sample-Z inverse pole direction lies too
    % close to the <111> family.
    rejectCone111 = 30 * degree;
    max_attempts  = 200;

    for k = 1:N_grains

        % Apply a random spin about sample Z to create the fiber texture.
        phi    = 2*pi*rand();
        R_spin = rotation.byAxisAngle(zvector, phi);

        accepted = false;
        R_wobble = rotation.id;

        for attempt = 1:max_attempts

            % Apply a small random wobble about an axis in the sample XY plane.
            % The wobble controls the angular spread around the ideal fiber.
            wobble_ang = (sigma_spread * max(-3, min(3, randn()))) * degree;

            % Random unit axis in the sample XY plane for the wobble rotation.
            az   = 2*pi*rand();
            ax_w = normalize(vector3d(cos(az), sin(az), 0));

            % Candidate wobble rotation.
            R_wobble_candidate = rotation.byAxisAngle(ax_w, wobble_ang);

            % Candidate total orientation mapping crystal coordinates to sample
            % coordinates.
            R_total_candidate = R_spin * R_wobble_candidate * R_base;

            % Compute the crystal direction parallel to sample Z for the candidate
            % orientation.
            h_ipfz = normalize(inv(R_total_candidate) * zvector);

            % Compute the minimum angular distance to the symmetrically equivalent
            % <111> directions.
            h_ipfz_sym = symmetrise(h_ipfz, cs);
            angs111 = angle(h_ipfz_sym, h111);
            minAng111 = min(angs111);

            % Accept the candidate if it is outside the rejected <111> cone.
            if minAng111 > rejectCone111
                R_wobble = R_wobble_candidate;
                accepted = true;
                break;
            end
        end

        if ~accepted
            % If no candidate is accepted, use the ideal fiber orientation
            % without wobble as a fallback.
            R_wobble = rotation.id;
        end

        % Store the final orientation as a Rodrigues vector.
        R_total = R_spin * R_wobble * R_base;
        ori_k   = orientation(R_total, cs);

        rod = ori_k.Rodrigues;
        idx = 3*(k-1) + 1;
        ori_vec(idx:idx+2) = [rod.x; rod.y; rod.z];

    end

elseif strcmp(orientation_type, 'random')
    % Generate uniformly random orientations using random unit quaternions,
    % then convert each quaternion to a Rodrigues vector.
    for k = 1:N_grains
        % Generate a random unit quaternion using the Shoemake method.
        u1 = rand(); u2 = rand(); u3 = rand();

        q_w = sqrt(1-u1) * sin(2*pi*u2);
        q_x = sqrt(1-u1) * cos(2*pi*u2);
        q_y = sqrt(u1)   * sin(2*pi*u3);
        q_z = sqrt(u1)   * cos(2*pi*u3);

        % Convert quaternion to Rodrigues vector.
        % A large finite value is used when the scalar quaternion component
        % is close to zero.
        if abs(q_w) < 1e-8
            rod_vec = [q_x q_y q_z] * sign(q_w + 1e-15) * 1e6;
        else
            rod_vec = [q_x q_y q_z] / q_w;
        end

        idx = 3*(k-1)+1;
        ori_vec(idx:idx+2) = rod_vec;
    end
end

%% 9. Save to HDF5

% Save the padded phase map and grain orientations to an HDF5 file for
% downstream simulation input.
h5filename = 'input_structure_poly.h5';

% Delete any existing file with the same name so h5create can write fresh
% datasets without name conflicts.
if isfile(h5filename)
    delete(h5filename);
    fprintf('Previous file "%s" deleted.\n', h5filename);
end


% --- Convert voxel phases to nodal phases ---
% phase2 is defined on voxels/cells with dimensions (Y, X, Z).
% The output /pix dataset is defined on grid nodes, so each dimension is
% increased by one. Boundary nodes are assigned the phase of the nearest
% voxel by clamping indices at the upper edge.
[Nyv2, Nxv2, Nzv2] = size(phase2);  
Nx_nodes = Nxv2 + 1;
Ny_nodes = Nyv2 + 1;
Nz_nodes = Nzv2 + 1;

iy = min((1:Ny_nodes), Nyv2);
ix = min((1:Nx_nodes), Nxv2);
iz = min((1:Nz_nodes), Nzv2);

pix_node = zeros(Ny_nodes, Nx_nodes, Nz_nodes, 'uint32');
for k = 1:Nz_nodes
    pix_node(:,:,k) = phase2(iy, ix, iz(k));
end

% Convert from MATLAB array order (Y, X, Z) to output order (X, Y, Z).
% Flattening pix_xyz(:) then gives X as the fastest-varying index.
pix_xyz  = permute(pix_node, [2 1 3]);
pix_flat = int32(pix_xyz(:));


% Write nodal phase IDs to /pix.
h5create(h5filename, '/pix', [length(pix_flat) 1], 'Datatype', 'int32');
h5write(h5filename, '/pix', pix_flat);
h5writeatt(h5filename, '/pix', 'dimensions', int32([Nx_nodes Ny_nodes Nz_nodes]));

% Write grain orientations to /orientation.
% Orientations are stored as Rodrigues vectors with three consecutive
% components per grain.
h5create(h5filename, '/orientation', size(ori_vec), 'Datatype', 'double');
h5write(h5filename, '/orientation', ori_vec);
h5writeatt(h5filename, '/orientation', 'dimensions', int32(size(ori_vec)));

disp(['Saved HDF5 input structure to ', h5filename]);


n_phases = max(phase(:));
% Save basic parameters to a MATLAB .mat file (for reloading or checking parameters)
save('xct_poly_params_SS316L.mat','xvec','yvec','zvec','Nx','Ny','Nz', ...
    'cube_len','N_grains','n_phases');
% Save all current workspace variables to a .mat file (version 7.3 for large arrays)
save('xct_poly_params_SS316L_full.mat','-v7.3')

%% 10. Visualization

% Build an IPF color map for the generated grain orientations.
% The final color entry is reserved for the pore/air phase.
if strcmp(orientation_type, 'textured')

    % Convert stored Rodrigues vectors to MTEX orientation objects.
    cs = crystalSymmetry('m-3m');
    rods = reshape(ori_vec(1:3*N_grains), [3, N_grains])';
    ori_mtex = orientation.byRodrigues(rods, cs);

    % Generate IPF colors using the sample Z direction.
    RD = vector3d.Z;
    ipf_key = ipfHSVKey(cs);
    ipf_key.inversePoleFigureDirection = RD;
    cmap_ipf = ipf_key.orientation2color(ori_mtex);
    % Append dark gray color for pore/air voxels.
    cmap_full = [cmap_ipf; [0.15 0.15 0.15]];


elseif strcmp(orientation_type, 'random')
    % Convert randomly generated Rodrigues vectors to MTEX orientations and
    % assign IPF colors.
    cs = crystalSymmetry('m-3m');
    disp('Crystal symmetry successfully defined.');
    rods = reshape(ori_vec(1:3*N_grains), [3, N_grains])';
    ori_mtex = orientation('rodrigues', rods, cs);
    ipf_key = ipfColorKey(cs);
    ipf_key.inversePoleFigureDirection = vector3d.Z;
    cmap_ipf = ipf_key.orientation2color(ori_mtex);
    % Append dark gray color for pore/air voxels.
    cmap_full = [cmap_ipf; [0.15 0.15 0.15]];
end


%--------------------------------------------------------------------------
% Pole figure visualization
% Randomly sample up to K_grain orientations to estimate and plot the ODF.
N = numel(ori_vec) / 3;
K = min(K_grain, N);
% Select a random subset of grains for pole-figure plotting.
idx_grains = randperm(N, K);

% Extract Rodrigues vectors for the sampled grains.
ori_idx = bsxfun(@plus, 3*(idx_grains'-1), (1:3));
ori_idx = ori_idx';
ori_idx = ori_idx(:);
ori_vec_downsampled = ori_vec(ori_idx);
ori_vec_matrix = reshape(ori_vec_downsampled, 3, []).';

% Convert sampled Rodrigues vectors to MTEX orientation objects.
ori = orientation.byRodrigues(ori_vec_matrix, cs, ss);

% Estimate the orientation distribution function using a
% von Mises-Fisher kernel.
psi = SO3vonMisesFisherKernel('halfwidth',10*degree);
% Estimate the ODF from the sampled orientations using the kernel
odf = calcDensity(ori, 'kernel', psi);

% Plot standard pole figures for {100}, {110}, and {111}.
plotPDF(odf, [Miller(1,0,0,cs), Miller(1,1,0,cs), Miller(1,1,1,cs)], ...
    'antipodal', 'contourf', 'smooth', ...
    'colorrange', 'equal', ...
    'resolution', 1*degree);
mtexColorbar('title', 'MRD');
mtexColorMap WhiteJet

%--------------------------------------------------------------------------
% XY slice at the top surface of the original cube.
% These plots visualize the original XCT cube stored in `phase`, not the
% padded HDF5 export volume stored in `phase2`.
figure; clf
s = squeeze(phase(:,:,Nzv));
imagesc(x_cent, y_cent, s);
axis equal tight
set(gca, 'YDir', 'normal');
colormap(cmap_full);
set(gca, 'CLim', [1 size(cmap_full,1)]);
xlabel('X_1 [mm]'); ylabel('X_2 [mm]');
set(gca, 'FontSize', 15, 'LineWidth', 1.1, 'Box', 'on')
xticks(linspace(min(xvec), max(xvec), 6));
yticks(linspace(min(yvec), max(yvec), 6));
xlim([min(xvec), max(xvec)]);
ylim([min(yvec), max(yvec)]);

%--------------------------------------------------------------------------
% YZ slice at the left surface of the original cube.
% These plots visualize the original XCT cube stored in `phase`, not the
% padded HDF5 export volume stored in `phase2`.
figure; clf
s = squeeze(phase(:,1,:));
imagesc(yvec, zvec, s');
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'normal');
axis equal tight
colormap(cmap_full);
set(gca, 'CLim', [1 size(cmap_full,1)]);
xlabel('X_2 [mm]'); ylabel('X_3 [mm]');
set(gca, 'FontSize', 15, 'LineWidth', 1.1, 'Box', 'on')

xticks(linspace(min(yvec), max(yvec), 6));
yticks(linspace(min(zvec), max(zvec), 6));
xlim([min(yvec), max(yvec)]);
ylim([min(zvec), max(zvec)]);

%--------------------------------------------------------------------------
% XZ slice at the front surface of the original cube.
% These plots visualize the original XCT cube stored in `phase`, not the
% padded HDF5 export volume stored in `phase2`.
figure; clf
s = squeeze(phase(1,:,:));
imagesc(xvec, zvec, s');
set(gca, 'YDir', 'normal'); axis equal tight
colormap(cmap_full);
set(gca, 'CLim', [1 size(cmap_full,1)]);
xlabel('X_1 [mm]'); ylabel('X_3 [mm]');

set(gca, 'FontSize', 15, 'LineWidth', 1.1, 'Box', 'on')
xticks(linspace(min(xvec), max(xvec), 6));
yticks(linspace(min(zvec), max(zvec), 6));
xlim([min(xvec), max(xvec)]);
ylim([min(zvec), max(zvec)]);

%--------------------------------------------------------------------------
% 3D surface visualization of the microstructure
% These plots visualize the original XCT cube stored in `phase`, not the
% padded HDF5 export volume stored in `phase2`.

% Identify voxels located on the outer boundary of the original phase map.
[Nx,Ny,Nz] = size(phase);
% Create a logical mask for boundary/surface voxels.
surface_mask = false(Nyv, Nxv, Nzv);
surface_mask(1,:,:)   = true;  surface_mask(end,:,:) = true;
surface_mask(:,1,:)   = true;  surface_mask(:,end,:) = true;
surface_mask(:,:,1)   = true;  surface_mask(:,:,end) = true;

% Flatten coordinate and phase arrays for scatter plotting.
phase_flat = phase(:);
X_flat = Xc(:);
Y_flat = Yc(:);
Z_flat = Zc(:);
surface_inds = find(surface_mask(:));
surface_inds = surface_inds(1:1:end);
grain_inds = phase_flat(surface_inds);

% Retain only voxels with valid phase IDs for color lookup.
valid = grain_inds >= 1 & grain_inds <= size(cmap_full,1);
surface_inds = surface_inds(valid);
grain_inds = grain_inds(valid);
pt_colors = cmap_full(grain_inds, :);

% Plot colored surface voxels.
figure; clf;
hold on;
scatter3(X_flat(surface_inds), Y_flat(surface_inds), Z_flat(surface_inds), 9, ...
    pt_colors, ...
    'filled', ...
    'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0);

% Optionally plot explicit spherical pores if pore-center and pore-radius
% variables are available in the workspace. This is separate from the
% voxelized pore mask already included in `phase`.
if exist('n_pores','var') && n_pores > 0
    [xx,yy,zz]=sphere(30);
    for i=1:n_pores
        surf(x(i)+r(i)*xx, y(i)+r(i)*yy, z(i)+r(i)*zz, ...
            'FaceAlpha',0.7,'EdgeColor','none','FaceColor',[0.4,0.4,0.6]);
    end
end

hx = xlabel('X_1 (mm)'); set(hx, 'Rotation', 21);
hy = ylabel('X_2 (mm)'); set(hy, 'Rotation', -21);
zlabel('X_3 (mm)');
view([-42.71 22.97]);
axis equal tight; grid on; box on;

xticks(linspace(min(xvec), max(xvec), 6));
yticks(linspace(min(yvec), max(yvec), 6));
zticks(linspace(min(zvec), max(zvec), 6));
xlim([min(xvec), max(xvec)]);
ylim([min(yvec), max(yvec)]);
zlim([min(zvec), max(zvec)]);

set(gca,'FontSize',15,'LineWidth',1.2);
hAxes = gca;
hAxes.BoxStyle      = "full";
hAxes.ClippingStyle = "3dbox";
hold off


%% ---- END OF SCRIPT ----
