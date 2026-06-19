%% Synthetic Polycrystal Generation with Imported Single-Pore XCT Region
% -------------------------------------------------------------------------
% This script generates a 3D voxel-based synthetic polycrystalline
% microstructure around an imported single-pore region. Grain IDs are assigned
% using randomly seeded anisotropic grains with lognormally distributed grain
% diameters and aspect ratios. The imported pore mask overwrites the generated
% grain phase map.
%
% The script also assigns crystallographic orientations to each grain and can
% export a selected subdomain as PRISMS-Plasticity input files:
%
%   - micro.msh : Gmsh 2.2 hexahedral mesh
%   - ori.txt   : local grain ID and Rodrigues-vector orientation file
%
% Main workflow:
%   1. Load imported pore-region data from PoreRegion_SinglePore_Rank*.mat.
%   2. Define the simulation grid and pore mask.
%   3. Generate lognormal grain sizes and aspect ratios.
%   4. Randomly seed grains and assign each grid site to the nearest grain using
%      an anisotropic ellipsoidal distance metric.
%   5. Overwrite pore sites with the pore/void phase ID.
%   6. Add XY pore/air padding and solid Z caps.
%   7. Write a voxel-volume summary.
%   8. Assign crystallographic orientations as Rodrigues vectors.
%   9. Export a selected PRISMS-Plasticity subdomain mesh and orientation file.
%  10. Generate diagnostic visualizations.
%
% Inputs:
%   - PoreRegion_SinglePore_Rank*_101x101x101.mat
%       MAT-file containing `pore_export_data`, including coordinate vectors,
%       voxel size, domain size, and logical pore mask.
%
% Outputs:
%   - phase
%       3D phase map for the imported pore-region grid.
%
%   - phase2
%       Padded phase map with XY pore/air padding and solid Z caps.
%
%   - ori_vec
%       Grain orientations stored as Rodrigues vectors.
%
%   - voxel_volume_summary.txt
%       Text summary of padded-domain voxel counts and phase IDs.
%
%   - micro.msh
%       Gmsh 2.2 mesh file for the selected PRISMS-Plasticity subdomain.
%
%   - ori.txt
%       Orientation file for grains present in the exported subdomain.
%
% Notes:
%   - Solid grains are labeled 1:N_grains.
%   - The pore/air phase is labeled N_grains + 1.
%   - Arrays use MATLAB storage convention phase(y_index,x_index,z_index).
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
% -------------------------------------------------------------------------
% Define target grain-size, grain-shape, and orientation settings for the
% synthetic polycrystalline microstructure.
%
% Units:
%   Grain-size parameters are specified in millimeters.
%
% Phase convention used later:
%   solid grains = 1:N_grains
%   pore/void    = N_grains + 1
% -------------------------------------------------------------------------

% --- Grain size statistics ---
% Mean and standard deviation of the grain diameter in the transverse XY plane.
mean_d = 0.048;              % mm
std_d  = 0.15*mean_d;        % mm

% --- Grain aspect-ratio statistics ---
% Aspect ratio is defined as:
%
%   aspect ratio = grain length in Z / grain width in XY
%
% Individual grain aspect ratios are sampled from a lognormal distribution.
mean_aspect = 2.832;
std_aspect = 0.10 * mean_aspect;

% --- Orientation settings ---
% orientation_type controls crystallographic orientation assignment:
%
%   'textured' : fiber-textured orientations with controlled angular spread
%   'random'   : uniformly random orientations
%
% In the current configuration, random orientations are assigned.
orientation_type = 'random';




%% 2. Load Imported Single-Pore Region and Define Grid
% -------------------------------------------------------------------------
% Load the selected single-pore region exported from the XCT/pore-processing
% workflow.
%
% The input MAT-file is expected to contain:
%
%   pore_export_data
%
% with fields including:
%   - pore_export_data.metadata.voxel_size_um
%   - pore_export_data.metadata.domain_size_um
%   - pore_export_data.x_coords_um
%   - pore_export_data.y_coords_um
%   - pore_export_data.z_coords_um
%   - pore_export_data.pore_volume
%
% Array convention throughout this script:
%
%   phase(y_index, x_index, z_index)
%   pore_mask(y_index, x_index, z_index)
%
% Coordinate vectors are stored separately as:
%
%   xvec(x_index), yvec(y_index), zvec(z_index)
% -------------------------------------------------------------------------

% Select which ranked pore region to import.
%
% target_rank:
%   Rank/index of the pore region to load.
%
% grid_num:
%   Number of intervals/grid cells per direction used when exporting the pore
%   region.
%
% export_dim:
%   Number of grid entries per direction in the exported pore region.
%   For this data format, export_dim = grid_num + 1.
target_rank = 1;
grid_num = 100;
export_dim  = grid_num + 1;

% Build input MAT-file name and load the exported pore-region structure.
mat_filename = sprintf('PoreRegion_SinglePore_Rank%d_%dx%dx%d.mat', ...
    target_rank, export_dim, export_dim, export_dim);

load(mat_filename,'pore_export_data');

% Number of grain orientations to sample for pole-figure visualization.
% This is only used later if pole-figure plotting is enabled.
K_grain = 2000;

% --- Domain and grid information -----------------------------------------
% Grid spacing from metadata.
%
% dx_um:
%   Grid spacing in micrometers.
%
% dx:
%   Grid spacing converted to millimeters.
dx_um = pore_export_data.metadata.voxel_size_um;
dx = dx_um * 1e-3;

% Cubic domain side length.
%
% cube_len_um:
%   Domain side length in micrometers.
%
% cube_len:
%   Domain side length converted to millimeters.
cube_len_um = pore_export_data.metadata.domain_size_um(1);
cube_len = cube_len_um * 1e-3;

% Coordinate vectors for the imported pore-region grid.
% These are converted from micrometers to millimeters for consistency with
% grain-size and domain-length parameters.
xvec_um = pore_export_data.x_coords_um;
yvec_um = pore_export_data.y_coords_um;
zvec_um = pore_export_data.z_coords_um;

xvec = xvec_um * 1e-3;   % mm
yvec = yvec_um * 1e-3;   % mm
zvec = zvec_um * 1e-3;   % mm

% Number of grid entries in each coordinate direction.
Nx = length(xvec);
Ny = length(yvec);
Nz = length(zvec);

% Build coordinate grids using the script's array storage convention:
%
%   array(y_index, x_index, z_index)
%
% Therefore, ndgrid is called as ndgrid(yvec, xvec, zvec), producing arrays
% Y, X, and Z with size [Ny, Nx, Nz].
[Y, X, Z] = ndgrid(yvec, xvec, zvec);

% Approximate physical domain volume in mm^3.
domain_vol = cube_len^3;

% Total number of grid sites in the imported pore region.
num_voxels = Nx*Ny*Nz;

% Import logical pore mask.
%
% Convention:
%   pore_mask = true  -> pore/void
%   pore_mask = false -> solid/grain material
%
% Expected array size:
%   [Ny, Nx, Nz]
pore_mask = pore_export_data.pore_volume;

% Verify that the pore mask follows the expected [Y, X, Z] storage convention.
assert(isequal(size(pore_mask), [Ny, Nx, Nz]), ...
    'Convention mismatch: expected pore_mask size [Ny Nx Nz] but got %s', ...
    mat2str(size(pore_mask)));


%% 3. Generate Lognormal Grain Sizes, Aspect Ratios, and Seeds
% -------------------------------------------------------------------------
% Estimate the number of grains needed for the imported pore-region domain,
% then sample individual grain sizes and aspect ratios from lognormal
% distributions.
%
% Grain geometry model:
%   - r_xy defines the transverse grain radius in the XY plane.
%   - r_z defines the longitudinal grain radius in the Z direction.
%
% These radii are used later in the anisotropic ellipsoidal distance metric
% for voxel/grid-site assignment.
% -------------------------------------------------------------------------

% Estimate the number of grains required to fill the domain.
%
% The mean grain volume is approximated using a cylindrical/ellipsoidal
% column-like volume based on:
%   - mean transverse diameter = mean_d
%   - mean longitudinal length = mean_aspect * mean_d
mean_grain_vol = pi*(mean_d/2)^2 * (mean_aspect*mean_d);
N_grains = ceil(domain_vol / mean_grain_vol);

% Convert target mean/std diameter into lognormal distribution parameters.
sigma_log = sqrt(log((std_d/mean_d)^2 + 1));
mu_log = log(mean_d) - 0.5*sigma_log^2;

% Sample grain diameters in the XY plane.
grain_diameters = lognrnd(mu_log, sigma_log, N_grains, 1);

% Transverse grain radius in XY.
r_xy = grain_diameters / 2;

% Convert target mean/std aspect ratio into lognormal distribution parameters.
sigma_log_a = sqrt(log((std_aspect/mean_aspect)^2 + 1));
mu_log_a    = log(mean_aspect) - 0.5*sigma_log_a^2;

% Sample aspect ratio for each grain.
aspect_ratios = lognrnd(mu_log_a, sigma_log_a, N_grains, 1);

% Longitudinal grain radius in Z.
r_z  = aspect_ratios .* r_xy;

% Randomly seed grain centers throughout the cubic domain.
%
% Units:
%   grain_seeds are in millimeters.
%
% Columns:
%   1 = X coordinate
%   2 = Y coordinate
%   3 = Z coordinate
grain_seeds = [rand(N_grains,1)*cube_len, ...
    rand(N_grains,1)*cube_len, ...
    rand(N_grains,1)*cube_len];



%% 4. Anisotropic Grain Assignment Using Parallel Block Processing
% -------------------------------------------------------------------------
% Assign each grid site to the nearest grain using an anisotropic ellipsoidal
% distance metric:
%
%   d^2 = (DeltaX/r_xy)^2 + (DeltaY/r_xy)^2 + (DeltaZ/r_z)^2
%
% where:
%   r_xy = grain radius in the transverse XY plane
%   r_z  = grain radius in the Z direction
%
% The grain with the smallest normalized distance owns the grid site. This
% creates elongated or equiaxed grains depending on the sampled aspect ratios.
%
% The computation is performed block-by-block and in parallel to reduce memory
% usage and improve runtime.
%
% Important convention:
%   phase(y_index, x_index, z_index) stores the assigned grain/phase ID.
% -------------------------------------------------------------------------

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

% Preallocate the phase array.
%
% phase = 0 initially means unassigned.
%
% Final convention after grain assignment:
%   1:N_grains      = solid grain IDs
%   N_grains + 1    = pore/void phase, assigned later
phase = zeros(Ny, Nx, Nz, 'uint32');

% Precompute inverse squared radii for faster distance calculations.
inv_r_xy2 = single(1 ./ (r_xy.^2));
inv_r_z2 = single(1 ./ (r_z.^2));

% Store grain seed coordinates in single precision to reduce memory use during
% repeated distance calculations.
grain_seeds_single = single(grain_seeds);

% Build a coarse spatial-grid diagnostic structure.
%
% Note:
%   The current block-assignment implementation filters candidate grains
%   directly from seed positions using block bounding boxes. The spatial_grid
%   built below is retained as a diagnostic/future acceleration structure but
%   is not currently used for candidate lookup inside the parfor block.
fprintf('\nBuilding spatial acceleration structure...\n');
accel_timer = tic;

grid_res = 15;
x_edges = linspace(min(xvec), max(xvec), grid_res+1);
y_edges = linspace(min(yvec), max(yvec), grid_res+1);
z_edges = linspace(min(zvec), max(zvec), grid_res+1);

% Maximum approximate grain influence distance used for candidate filtering.
%
% The factor 1.5 provides a buffer around each block so that grains with large
% radii are still considered even if their seeds lie outside the block bounds.
max_influence = max([r_xy(:); r_z(:)]) * 1.5;
total_cells = grid_res^3;

% Temporary flat cell array for parfor-compatible spatial-grid construction.
spatial_grid_temp = cell(total_cells, 1);

% Build the diagnostic spatial grid in parallel using linear indexing.
% Each cell stores grain IDs whose seed locations lie close enough to possibly
% influence that coarse spatial cell.
parfor linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);

    cell_xmin = x_edges(i); cell_xmax = x_edges(i+1);
    cell_ymin = y_edges(j); cell_ymax = y_edges(j+1);
    cell_zmin = z_edges(k); cell_zmax = z_edges(k+1);

    % Find grains whose seed locations fall inside the expanded coarse-cell
    % bounding box.
    relevant_grains = find(...
        grain_seeds_single(:,1) >= (cell_xmin - max_influence) & ...
        grain_seeds_single(:,1) <= (cell_xmax + max_influence) & ...
        grain_seeds_single(:,2) >= (cell_ymin - max_influence) & ...
        grain_seeds_single(:,2) <= (cell_ymax + max_influence) & ...
        grain_seeds_single(:,3) >= (cell_zmin - max_influence) & ...
        grain_seeds_single(:,3) <= (cell_zmax + max_influence));

    spatial_grid_temp{linear_idx} = relevant_grains;
end

% Convert flat parfor output into a 3D cell array indexed by coarse spatial bin.
spatial_grid = cell(grid_res, grid_res, grid_res);
for linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);
    spatial_grid{i,j,k} = spatial_grid_temp{linear_idx};
end
clear spatial_grid_temp;

fprintf('Spatial acceleration structure built in %.1fs\n', toc(accel_timer));

% --- Block-processing setup ----------------------------------------------
% Divide the full [Ny, Nx, Nz] phase grid into smaller spatial blocks.
%
% Purpose:
%   Processing the entire 3D grid at once would require large temporary
%   coordinate and distance arrays. Block processing limits peak memory usage.
block_size = 50;

n_blocks_x = ceil(Nx / block_size);
n_blocks_y = ceil(Ny / block_size);
n_blocks_z = ceil(Nz / block_size);
total_blocks = n_blocks_x * n_blocks_y * n_blocks_z;

fprintf('Processing %d blocks (%dx%dx%d) in parallel:\n', ...
    total_blocks, n_blocks_x, n_blocks_y, n_blocks_z);

% Group blocks into batches.
%
% Each batch is processed using parfor. Batching allows progress reporting
% between parfor calls, since progress cannot be printed reliably from inside
% a parfor loop.
num_workers = gcp().NumWorkers;
blocks_per_batch = max(1, floor(total_blocks / (num_workers * 4)));
n_batches = ceil(total_blocks / blocks_per_batch);

fprintf('Using %d batches of %d blocks each for progress tracking\n', ...
    n_batches, blocks_per_batch);

% Initialize progress and performance counters.
completed_blocks = 0;
total_grains_processed = 0;
blocks_with_grains = 0;
batch_times = [];


fprintf('\nStarting parallel processing with real-time progress:\n');
process_timer = tic;

for batch_idx = 1:n_batches
    batch_timer = tic;

    % Determine global block-number range included in this batch.
    batch_start = (batch_idx - 1) * blocks_per_batch + 1;
    batch_end = min(batch_idx * blocks_per_batch, total_blocks);
    current_batch_size = batch_end - batch_start + 1;

    % Build a list of block indices for this batch.
    %
    % block_list_batch columns:
    %   1 = block index in X
    %   2 = block index in Y
    %   3 = block index in Z
    %   4 = global block number
    block_list_batch = zeros(current_batch_size, 4);
	
    for local_idx = 1:current_batch_size
        block_num = batch_start + local_idx - 1;
		
        bz = ceil(block_num / (n_blocks_x * n_blocks_y));
        remaining = block_num - (bz-1) * n_blocks_x * n_blocks_y;
        by = ceil(remaining / n_blocks_x);
        bx = remaining - (by-1) * n_blocks_x;
		
        block_list_batch(local_idx, :) = [bx, by, bz, block_num];
    end

    % Preallocate containers for results from each block in the batch.
    %
    % batch_results{local_idx}:
    %   local block phase assignment array.
    %
    % batch_indices{local_idx}:
    %   [ix1 ix2 iy1 iy2 iz1 iz2] index bounds for writing the block result
    %   back into the global phase array.
    %
    % batch_grain_counts(local_idx):
    %   number of candidate grains considered for that block.
    batch_results = cell(current_batch_size, 1);
    batch_indices = cell(current_batch_size, 1);
    batch_grain_counts = zeros(current_batch_size, 1);

    % Local copies for parfor classification.
    xvec_local = xvec;
    yvec_local = yvec;
    zvec_local = zvec;
	
    % These edge arrays are retained for diagnostic grid-index calculations
    % below. Candidate grains are still filtered directly from seed positions.
    x_edges_local = x_edges;
    y_edges_local = y_edges;
    z_edges_local = z_edges;

    % Process all blocks in the current batch in parallel.
    parfor local_idx = 1:current_batch_size
        bx = block_list_batch(local_idx, 1);
        by = block_list_batch(local_idx, 2);
        bz = block_list_batch(local_idx, 3);

        % Determine grid-index bounds for the current block.
        %
        % Array convention:
        %   phase(y_index, x_index, z_index)
        ix1 = (bx-1)*block_size + 1;
        ix2 = min(bx*block_size, Nx);
        iy1 = (by-1)*block_size + 1;
        iy2 = min(by*block_size, Ny);
        iz1 = (bz-1)*block_size + 1;
        iz2 = min(bz*block_size, Nz);

        % Determine physical coordinate bounds of the current block.
        block_xmin = xvec_local(ix1); block_xmax = xvec_local(ix2);
        block_ymin = yvec_local(iy1); block_ymax = yvec_local(iy2);
        block_zmin = zvec_local(iz1); block_zmax = zvec_local(iz2);

        % Compute coarse spatial-grid cell range overlapped by this block.
        %
        % Note:
        %   gx1/gx2/etc. are currently diagnostic only. Candidate grains are
        %   selected by the direct bounding-box filter below.
        gx1 = max(1, floor((block_xmin - min(xvec_local)) / (max(xvec_local) - min(xvec_local)) * grid_res) + 1);
        gx2 = min(grid_res, ceil((block_xmax - min(xvec_local)) / (max(xvec_local) - min(xvec_local)) * grid_res));
        gy1 = max(1, floor((block_ymin - min(yvec_local)) / (max(yvec_local) - min(yvec_local)) * grid_res) + 1);
        gy2 = min(grid_res, ceil((block_ymax - min(yvec_local)) / (max(yvec_local) - min(yvec_local)) * grid_res));
        gz1 = max(1, floor((block_zmin - min(zvec_local)) / (max(zvec_local) - min(zvec_local)) * grid_res) + 1);
        gz2 = min(grid_res, ceil((block_zmax - min(zvec_local)) / (max(zvec_local) - min(zvec_local)) * grid_res));

        % Identify candidate grains whose seed locations are close enough to
        % potentially influence this block.
        %
        % The block bounding box is expanded by max_influence to avoid checking
        % all grains while still including grains with large radii.
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
            % Create coordinate arrays for the current block.
            %
            % Array size:
            %   [iy2-iy1+1, ix2-ix1+1, iz2-iz1+1]
            %
            % Coordinate convention:
            %   Xb corresponds to x-coordinate,
            %   Yb corresponds to y-coordinate,
            %   Zb corresponds to z-coordinate.
            [Yb, Xb, Zb] = ndgrid(single(yvec_local(iy1:iy2)), ...
                single(xvec_local(ix1:ix2)), ...
                single(zvec_local(iz1:iz2)));

            % Initialize local phase assignment and nearest-distance arrays.
            block_phase = zeros(size(Yb), 'uint32');
            min_dist2   = inf(size(Yb), 'single');

            % Assign each grid site in the block to the candidate grain with
            % the minimum anisotropic normalized distance.
            for k = relevant_grains'
                dX = Xb - grain_seeds_single(k,1);
                dY = Yb - grain_seeds_single(k,2);
                dZ = Zb - grain_seeds_single(k,3);

                % Anisotropic squared distance to grain seed k.
                current_dist2 = (dX.^2) * inv_r_xy2(k) + (dY.^2) * inv_r_xy2(k) + (dZ.^2) * inv_r_z2(k);

                % Update assignment where grain k is the closest candidate so far.
                update_mask = current_dist2 < min_dist2;
                block_phase(update_mask) = k;
                min_dist2(update_mask) = current_dist2(update_mask);
            end

            batch_results{local_idx} = block_phase;
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = length(relevant_grains);
        else
            % No candidate grains were found for this block.
            % This should be rare because seeds are generated throughout the domain.
            batch_results{local_idx} = [];
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = 0;
        end
    end

    % Copy completed block results from this batch back into the global phase
    % array.
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
		
        % Count total candidate-grain evaluations for this batch.
        batch_grains_processed = batch_grains_processed + batch_grain_counts(local_idx);
    end

    % Update cumulative progress counters.
    completed_blocks = completed_blocks + current_batch_size;
    total_grains_processed = total_grains_processed + batch_grains_processed;
    blocks_with_grains = blocks_with_grains + batch_blocks_with_grains;

    % Record batch timing.
    batch_time = toc(batch_timer);
    batch_times(end+1) = batch_time;

    % Compute progress statistics.
    progress_pct = 100 * completed_blocks / total_blocks;
    elapsed_total = toc(main_timer);
    avg_grains_per_block = batch_grains_processed / current_batch_size;

    % Estimate remaining runtime.
    if batch_idx > 1
        avg_batch_time = mean(batch_times);
        remaining_batches = n_batches - batch_idx;
        eta_seconds = remaining_batches * avg_batch_time;
    else
        eta_seconds = batch_time * (n_batches - 1);
    end


    % Print batch-level progress summary.
    fprintf('Batch %d/%d: Blocks %d-%d (%.1f%%) | %.1f grains/block avg | Elapsed: %.1fs | ETA: %.1fs (%.1f min)\n', ...
        batch_idx, n_batches, batch_start, batch_end, progress_pct, avg_grains_per_block, elapsed_total, eta_seconds, eta_seconds/60);

    % Print detailed performance statistics periodically.
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
processing_rate = double(Nx*Ny*Nz) / total_time;

fprintf('\n=== COMPLETION SUMMARY ===\n');
fprintf('Total time: %.1f seconds (%.2f minutes)\n', total_time, total_time/60);
fprintf('Processing rate: %s voxels/second\n', regexprep(sprintf('%.0f', processing_rate),'(\d)(?=(\d{3})+$)','$1,'));

% Check whether any grid sites remained unassigned after grain partitioning.
unassigned = nnz(phase(:) == 0);

if unassigned > 0
    fprintf('\nWARNING: %d unassigned voxels detected!\n', unassigned);
else
    fprintf('SUCCESS: All %s voxels assigned successfully!\n', ...
	    regexprep(sprintf('%.0f',double(Nx*Ny*Nz)),'(\d)(?=(\d{3})+$)','$1,'));
end

fprintf('\nParallel + Spatial filtering microstructure generation completed!\n');



%% 5. Apply Imported Pore Mask
% -------------------------------------------------------------------------
% Overwrite the synthetic grain assignments with the imported pore/void phase.
%
% Phase convention:
%   solid grains = 1:N_grains
%   pore/void    = N_grains + 1
%
% The pore mask is imported from pore_export_data.pore_volume and uses the
% same [Y, X, Z] storage convention as phase.
% -------------------------------------------------------------------------

phase(pore_mask) = N_grains + 1;
fprintf('Phase assignment completed.\n');


%% 6. Add XY Pore/Air Padding and Solid Z Caps
% -------------------------------------------------------------------------
% Expand the generated phase map by adding:
%
%   1. An XY pore/air padding ring around the original imported domain.
%   2. Solid top and bottom Z caps generated by extruding the solid footprints
%      from the first and last interior slices.
%
% Purpose:
%   The padding/caps create a simulation-ready domain with surrounding pore/air
%   in the lateral directions and solid material at the Z ends.
%
% Important:
%   The same phase ID is used for physical pores and artificial pore/air
%   padding:
%
%       pore_phase_id = N_grains + 1
%
% Array convention:
%   phase2(y_index, x_index, z_index)
% -------------------------------------------------------------------------

% Padding thicknesses in grid sites.
%
% pad_xy:
%   Number of pore/air padding layers added on each side in X and Y.
%
% pad_z:
%   Number of solid cap layers added at the bottom and top in Z.
pad_xy = round((Nx-1)*0.02);
pad_z  = round((Nx-1)*0.02);

% Pore/air phase ID.
pore_phase_id = cast(N_grains + 1, 'uint32');

% Original phase-map dimensions.
[Nyv, Nxv, Nzv] = size(phase);

% Padded phase-map dimensions.
Nxv2 = Nxv + 2*pad_xy;
Nyv2 = Nyv + 2*pad_xy;
Nzv2 = Nzv + 2*pad_z;

% Initialize the expanded domain as pore/air.
phase2 = repmat(pore_phase_id, [Nyv2, Nxv2, Nzv2]);

% Insert original generated microstructure into the center of the padded domain.
xs = (pad_xy + 1):(pad_xy + Nxv);
ys = (pad_xy + 1):(pad_xy + Nyv);
zs = (pad_z  + 1):(pad_z  + Nzv);

phase2(ys, xs, zs) = phase;

% Extract the bottom and top interior slices used to define the solid cap
% footprints.
ref_bot = phase2(ys, xs, pad_z+1);
ref_top = phase2(ys, xs, Nzv2-pad_z);

% Identify solid footprints on the bottom and top interior slices.
% Pore/air sites are excluded from the cap footprint.
foot_bot = (ref_bot ~= pore_phase_id);
foot_top = (ref_top ~= pore_phase_id);

% Check whether each solid footprint is connected.
% Disconnected footprints may indicate that the imported pore geometry cuts
% through the face or that the selected domain has multiple solid islands.
cc_bot = bwconncomp(foot_bot);
cc_top = bwconncomp(foot_top);

if cc_bot.NumObjects > 1
    warning('Bottom cap footprint has %d disconnected regions.', cc_bot.NumObjects);
end

if cc_top.NumObjects > 1
    warning('Top cap footprint has %d disconnected regions.', cc_top.NumObjects);
end

% Extract solid grain IDs on the bottom and top interior slices.
% These representative IDs are written to the summary file for traceability.
ref_bot_solid = ref_bot(foot_bot);
ref_top_solid = ref_top(foot_top);

if isempty(ref_bot_solid)
    error('Bottom gauge slice has no solid voxels in core window.');
end

if isempty(ref_top_solid)
    error('Top gauge slice has no solid voxels in core window.');
end

% Most frequent grain IDs on the bottom and top interior slices.
% Note:
%   These IDs are reported only for summary/diagnostic purposes. The cap
%   extrusion below copies the full face grain pattern, not just these mode IDs.
fill_id_bot = cast(mode(double(ref_bot_solid(:))), 'uint32');
fill_id_top = cast(mode(double(ref_top_solid(:))), 'uint32');

% Initialize Z cap regions to pore/air before solid-footprint extrusion.
phase2(ys, xs, 1:pad_z)            = pore_phase_id;
phase2(ys, xs, Nzv2-pad_z+1:Nzv2) = pore_phase_id;

% Extrude the bottom interior face pattern into the lower Z cap.
% Solid footprint sites copy their grain IDs from ref_bot. Pore/air sites remain
% pore/air.
for k = 1:pad_z
    cap_slice = pore_phase_id * ones(size(ref_bot), 'like', ref_bot);
    cap_slice(foot_bot) = ref_bot(foot_bot);
    phase2(ys, xs, k) = cap_slice;
end

% Extrude the top interior face pattern into the upper Z cap.
% Solid footprint sites copy their grain IDs from ref_top. Pore/air sites remain
% pore/air.
for k = (Nzv2 - pad_z + 1):Nzv2
    cap_slice = pore_phase_id * ones(size(ref_top), 'like', ref_top);
    cap_slice(foot_top) = ref_top(foot_top);
    phase2(ys, xs, k) = cap_slice;
end


%% 7. Write Voxel/Grid-Site Volume Summary
% -------------------------------------------------------------------------
% Write a text summary of the original and padded domain dimensions, padding
% thicknesses, phase IDs, and solid/pore-air counts.
%
% Note:
%   The term "voxel" here refers to one entry of the phase array. For this
%   imported pore-region format, these may also be interpreted as grid sites.
%
% Important:
%   The geometric counts pad_xy_vox and pad_z_vox are not mutually exclusive
%   categories if interpreted as padding regions; their percentages should not
%   be summed.
% -------------------------------------------------------------------------
summary_file = fullfile(pwd, 'voxel_volume_summary.txt');
fid = fopen(summary_file, 'w');
if fid == -1
    error('Could not open summary file for writing: %s', summary_file);
end

% Total number of entries in the padded phase map.
total_voxels = numel(phase2);

% Number of entries in the original unpadded phase map.
interior_vox = Nxv * Nyv * Nzv;

% Number of entries in the XY pore/air padding ring over the full padded
% Z height.
pad_xy_vox = total_voxels - (Nxv * Nyv * Nzv2);

% Number of entries in the top and bottom Z padding layers over the full
% padded XY footprint.
pad_z_vox    = Nxv2 * Nyv2 * 2 * pad_z;

% Count pore/air and solid entries in the padded phase map.
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
% -------------------------------------------------------------------------
% Assign one crystallographic orientation to each generated solid grain.
%
% Output:
%   ori_vec
%     Column vector containing one Rodrigues vector per grain:
%
%       [r1_x; r1_y; r1_z; r2_x; r2_y; r2_z; ...]
%
% Orientation modes:
%   'textured'
%     Generates a fiber-textured orientation distribution using MTEX rotations.
%
%   'random'
%     Generates uniformly random orientations using random unit quaternions.
%
% Note:
%   Only solid grain phases 1:N_grains receive orientations. The pore/air
%   phase, N_grains+1, does not receive an orientation.
% -------------------------------------------------------------------------

% Store three Rodrigues components per grain.
ori_vec = zeros(3*N_grains, 1);

if strcmp(orientation_type, 'textured')
    % Define cubic crystal symmetry for the material.
    cs = crystalSymmetry('m-3m');

    % Target fiber direction in the crystal frame.
    fiber_axis = normalize(vector3d(1,0,1));

    % Select a crystal-frame reference direction that is not parallel to the
    % fiber axis. This is used to define a stable base rotation.
    ref = vector3d(0,1,0);
    if angle(fiber_axis, ref) < 1*degree
        ref = vector3d(1,0,0);
    end

    % Define a crystal-frame direction perpendicular to the fiber axis and a
    % sample-frame reference direction.
    %
    % These two direction mappings fix the base rotation, including rotation
    % about the fiber axis.
    c_perp = normalize(cross(fiber_axis, ref));
    s_perp = yvector;

    % Construct base rotation mapping:
    %
    %   crystal fiber_axis -> sample Z
    %   crystal c_perp     -> sample Y
    R_base = rotation.map(fiber_axis, zvector, c_perp, s_perp);

    % Reference <111> crystal direction used for rejection filtering.
    h111 = normalize(vector3d(1,1,1));

    % Reject candidate orientations whose sample-Z inverse-pole direction lies
    % too close to the <111> family.
    rejectCone111 = 30 * degree;
    max_attempts  = 200;

    for k = 1:N_grains

        % Apply a random spin about sample Z to create the fiber texture.
        phi    = 2*pi*rand();
        R_spin = rotation.byAxisAngle(zvector, phi);

        accepted = false;
        R_wobble = rotation.id;

        for attempt = 1:max_attempts

            % Apply a small random wobble about a random axis in the sample
            % XY plane. The wobble controls angular spread around the ideal
            % fiber orientation.
            %
            % The random normal value is clipped to ±3 sigma.
            wobble_ang = (sigma_spread * max(-3, min(3, randn()))) * degree;

            % Random unit axis in the sample XY plane for the wobble rotation.
            az   = 2*pi*rand();
            ax_w = normalize(vector3d(cos(az), sin(az), 0));

            % Candidate wobble rotation.
            R_wobble_candidate = rotation.byAxisAngle(ax_w, wobble_ang);

            % Candidate total orientation mapping crystal coordinates to sample
            % coordinates.
            R_total_candidate = R_spin * R_wobble_candidate * R_base;

            % Compute the crystal direction parallel to sample Z for this
            % candidate orientation.
            h_ipfz = normalize(inv(R_total_candidate) * zvector);

            % Compute the minimum angular distance between the candidate
            % inverse-pole direction and the symmetrically equivalent <111>
            % family.
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
            % If no candidate is accepted, fall back to the ideal fiber
            % orientation without wobble.
            R_wobble = rotation.id;
        end

        % Final grain orientation.
        R_total = R_spin * R_wobble * R_base;
        ori_k   = orientation(R_total, cs);

        % Store orientation as Rodrigues vector.
        rod = ori_k.Rodrigues;
        idx = 3*(k-1) + 1;
        ori_vec(idx:idx+2) = [rod.x; rod.y; rod.z];

    end

elseif strcmp(orientation_type, 'random')
    % Generate uniformly random orientations using random unit quaternions.
    %
    % The Shoemake method samples a uniform distribution on the unit quaternion
    % sphere, which corresponds to uniformly random rotations in SO(3).
    for k = 1:N_grains
	
        % Random unit quaternion components.
        u1 = rand(); u2 = rand(); u3 = rand();

        q_w = sqrt(1-u1) * sin(2*pi*u2);
        q_x = sqrt(1-u1) * cos(2*pi*u2);
        q_y = sqrt(u1)   * sin(2*pi*u3);
        q_z = sqrt(u1)   * cos(2*pi*u3);

        % Convert quaternion to Rodrigues vector:
        %
        %   r = q_vector / q_scalar
        %
        % If q_w is close to zero, the Rodrigues vector magnitude becomes very
        % large. Use a large finite value to avoid Inf/NaN values.
        if abs(q_w) < 1e-8
            rod_vec = [q_x q_y q_z] * sign(q_w + 1e-15) * 1e6;
        else
            rod_vec = [q_x q_y q_z] / q_w;
        end

        % Store Rodrigues vector in ori_vec.
        idx = 3*(k-1)+1;
        ori_vec(idx:idx+2) = rod_vec;
    end
end


%% 9. Export PRISMS-Plasticity Mesh and Orientation Files
% -------------------------------------------------------------------------
% Export a selected subdomain as PRISMS-Plasticity input files:
%
%   - micro.msh : Gmsh 2.2 ASCII hexahedral mesh
%   - ori.txt   : local grain ID and Rodrigues-vector orientation file
%
% Inputs from previous sections:
%   phase
%     Unpadded phase map for the imported pore-region grid, stored as
%     phase(y_index, x_index, z_index).
%
%   ori_vec
%     Rodrigues orientation vectors for all generated solid grains.
%
%   N_grains
%     Number of solid grains in the generated microstructure.
%
%   xvec, yvec, zvec
%     Coordinate vectors in millimeters.
%
% Phase convention:
%   solid grains = 1:N_grains
%   pore/void    = N_grains + 1
%
% Pore elements are skipped during mesh export.
%
% Note:
%   This PRISMS export uses the unpadded phase map `phase`, not the padded
%   phase map `phase2`.
% -------------------------------------------------------------------------

%% ---------------------- EXPORT SETTINGS ---------------------------------
% Output file names for PRISMS-Plasticity input.
%
% micro.msh:
%   Gmsh 2.2 ASCII mesh file containing 8-node hexahedral elements.
%
% ori.txt:
%   Orientation file containing compact local grain IDs and Rodrigues vectors
%   for grains present in the exported subdomain.
msh_filename = 'micro.msh';
ori_filename = 'ori.txt';

disp(['Writing GMSH mesh file: ', msh_filename])

%% ---------------------- SELECT SUBDOMAIN --------------------------------
% Define the physical bounds of the PRISMS-Plasticity subdomain.
%
% Units:
%   All coordinate ranges are in millimeters.
%
% Coordinate convention:
%   x_range = [xmin, xmax]
%   y_range = [ymin, ymax]
%   z_range = [zmin, zmax]
%
% Selection convention used by this code:
%   Grid sites are selected using inclusive coordinate bounds:
%
%       range_min <= coordinate <= range_max
%
% This matches the find(...) statements below.
%
% Important:
%   If xvec/yvec/zvec represent voxel-center or grid-site coordinates, the
%   selected phase_sub directly corresponds to those selected grid sites.
%
%   If xvec/yvec/zvec represent node coordinates, inclusive selection may not
%   correspond to a strictly voxel-contained subdomain. In that case, adjust
%   the selection rule to use left-node vectors and [range_min, range_max).
x_range = [0.0275, 0.5225]; % mm
y_range = [0.0275, 0.5225]; % mm
z_range = [0.0275, 0.5225]; % mm

% Select coordinate indices inside the requested bounds.
ix = find(xvec >= x_range(1) & xvec <= x_range(2));
iy = find(yvec >= y_range(1) & yvec <= y_range(2));
iz = find(zvec >= z_range(1) & zvec <= z_range(2));

if isempty(ix) || isempty(iy) || isempty(iz)
    error('Requested x/y/z range selects an empty subdomain. Check units of xvec/yvec/zvec vs x_range.');
end

% Report requested physical size.
fprintf('Requested domain size: %.5f x %.5f x %.5f mm^3\n', ...
    diff(x_range), diff(y_range), diff(z_range));

%% ---------------------- EXTRACT SUBDOMAIN PHASE MAP ---------------------
% Extract coordinate vectors for the selected subdomain.
xvec_sub  = xvec(ix);
yvec_sub  = yvec(iy);
zvec_sub  = zvec(iz);

% Extract phase map for the selected subdomain.
%
% Array convention:
%   phase_sub(y_index, x_index, z_index)
phase_sub = phase(iy, ix, iz);

% Subdomain grid-site/element counts.
Ny_sub = length(yvec_sub);
Nx_sub = length(xvec_sub);
Nz_sub = length(zvec_sub);

fprintf('Core subdomain voxels: Nx=%d x Ny=%d x Nz=%d\n', ...
    Nx_sub, Ny_sub, Nz_sub);

%% ---------------------- LOCAL GRAIN-ID REMAPPING -------------------------
% Build a compact local grain-ID system for the exported PRISMS-Plasticity
% subdomain.
%
% Full-domain phase convention:
%   solid grains = 1:N_grains
%   pore/void    = N_grains + 1
%
% Export/subdomain convention:
%   - Only solid grains present in the selected subdomain are retained.
%   - Present solid grains are remapped to compact local IDs:
%
%         1:N_grains_sub
%
%   - A local pore ID is defined as:
%
%         pore_id_sub = N_grains_sub + 1
%
%     but pore elements are skipped when writing the Gmsh mesh.
%
% Why remap?
%   PRISMS-Plasticity orientation files are easier to manage when grain IDs are
%   compact and correspond directly to the rows in ori.txt.
%
% Note:
%   The grain list is built from phase_sub, so only grains present in the
%   selected exported subdomain are included in ori.txt.
% -------------------------------------------------------------------------

% Identify solid grain entries in the selected subdomain.
mask_grain_global = (phase_sub >= 1) & (phase_sub <= N_grains);

% Sorted list of full-domain/global grain IDs present in this subdomain.
grain_ids_global  = unique(phase_sub(mask_grain_global));
grain_ids_global  = sort(grain_ids_global(:));

% Number of solid grains present in the exported subdomain.
N_grains_sub = numel(grain_ids_global);

% Local pore ID used internally for bookkeeping.
% Pore elements are not written to the mesh.
pore_id_sub  = N_grains_sub + 1;

fprintf('Grains in domain: %d (pore phase in export = %d)\n', ...
    N_grains_sub, pore_id_sub);

% Create a lookup table from global grain ID to compact local grain ID.
%
% Example:
%   global2local(global_id) = local_id
%
% Grain IDs not present in this subdomain remain zero.
global2local = zeros(N_grains, 1, 'uint32');
global2local(grain_ids_global) = uint32(1:N_grains_sub);

% Convert the selected subdomain phase map from global grain IDs to compact
% local grain IDs.
%
% All entries are initialized as the local pore ID. Solid grain entries are
% then overwritten with their corresponding local grain IDs.
phase_sub_local = uint32(pore_id_sub) * ones(size(phase_sub),'uint32');
phase_sub_local(mask_grain_global) = global2local(phase_sub(mask_grain_global));


%% ---------------------- BUILD MESH NODE COORDINATES ----------------------
% Construct node coordinate vectors for the exported Gmsh mesh.
%
% The selected phase_sub array contains:
%
%   Nx_sub x Ny_sub x Nz_sub
%
% hexahedral elements/grid cells. Therefore, the mesh requires:
%
%   Nx_nodes = Nx_sub + 1
%   Ny_nodes = Ny_sub + 1
%   Nz_nodes = Nz_sub + 1
%
% nodes in each direction.
%
% Coordinate construction used here:
%   Node coordinates are generated uniformly from the requested physical bounds:
%
%       x_range, y_range, z_range
%
%   For example, Nx_sub elements in X require Nx_sub+1 nodes between
%   x_range(1) and x_range(2).
%
% Units:
%   Coordinates are in millimeters.
%
% Output convention:
%   The exported mesh is shifted so that the minimum corner of the selected
%   subdomain is at the origin.
% -------------------------------------------------------------------------

% Build uniformly spaced physical node coordinates from the requested bounds.
xnodes_phys = linspace(x_range(1), x_range(2), Nx_sub + 1);
ynodes_phys = linspace(y_range(1), y_range(2), Ny_sub + 1);
znodes_phys = linspace(z_range(1), z_range(2), Nz_sub + 1);

% Shift the coordinate origin to the minimum corner of the selected subdomain.
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
%   ix varies fastest, then iy, then iz.
%
% This ordering must match both:
%   - the node-writing loop below, and
%   - the hexahedral element connectivity in the element-writing section.
nodeid = @(ix,iy,iz) (iz-1)*Ny_nodes*Nx_nodes + (iy-1)*Nx_nodes + ix;

%% ---------------------- WRITE GMSH NODES --------------------------------
% Write the node block of the Gmsh 2.2 ASCII mesh file.
%
% Gmsh node format:
%
%   node_id  x  y  z
%
% Coordinates are written in millimeters and have already been shifted so that
% the exported subdomain starts at the origin.
%
% Node ordering:
%   ix varies fastest, then iy, then iz.
% -------------------------------------------------------------------------

fid = fopen(msh_filename, 'w');

% Gmsh 2.2 ASCII mesh header.
fprintf(fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');

% Write node count and node records.
fprintf(fid, '$Nodes\n%d\n', num_nodes);

% Write nodal coordinates.
nid = 1;
for iz=1:Nz_nodes
    for iy=1:Ny_nodes
        for ix=1:Nx_nodes
            fprintf(fid, '%d %.10g %.10g %.10g\n', ...
			    nid, xnodes(ix), ynodes(iy), znodes(iz));
            nid = nid+1;
        end
    end
end

fprintf(fid, '$EndNodes\n');
fclose(fid);

%% ---------------------- WRITE GMSH ELEMENTS ------------------------------
% Write the element block of the Gmsh 2.2 ASCII mesh file.
%
% Each non-pore entry in phase_sub_local becomes one 8-node hexahedral element.
% Pore entries are skipped and therefore do not appear in the exported mesh.
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
%     Number of element tags.
%
%   physical_tag, elementary_tag
%     Both are set to the compact local grain ID. This allows the grain identity
%     to be recovered from the element tags.
%
%   n1...n8
%     Node IDs of the hexahedral element.
%
% Element connectivity convention:
%   n1 = (ix,   iy,   iz)
%   n2 = (ix+1, iy,   iz)
%   n3 = (ix+1, iy+1, iz)
%   n4 = (ix,   iy+1, iz)
%   n5 = (ix,   iy,   iz+1)
%   n6 = (ix+1, iy,   iz+1)
%   n7 = (ix+1, iy+1, iz+1)
%   n8 = (ix,   iy+1, iz+1)
%
% This node ordering should be verified against the element orientation expected
% by PRISMS-Plasticity/Gmsh for the intended material model and boundary
% conditions.
% -------------------------------------------------------------------------

fid = fopen(msh_filename, 'a');

% Count and store element records before writing the $Elements block, because
% Gmsh requires the total number of elements before the element records.
%
% Preallocate for the maximum possible number of exported elements. The actual
% number may be smaller because pore entries are skipped.
elements_written = 0;
eid              = 1;
element_lines = cell(Nx*Ny*Nz,1);

for iz = 1:Nz_sub
    for iy = 1:Ny_sub
        for ix = 1:Nx_sub

            % Compact local grain ID for this subdomain cell.
            gID = phase_sub_local(iy, ix, iz);

            % Skip pore/void entries. Only solid grain elements are exported.           
            if gID == pore_id_sub
                continue;
            end

            % Hexahedral element node connectivity.
            n1 = nodeid(ix  ,iy  ,iz  );
            n2 = nodeid(ix+1,iy  ,iz  );
            n3 = nodeid(ix+1,iy+1,iz  );
            n4 = nodeid(ix  ,iy+1,iz  );
            n5 = nodeid(ix  ,iy  ,iz+1);
            n6 = nodeid(ix+1,iy  ,iz+1);
            n7 = nodeid(ix+1,iy+1,iz+1);
            n8 = nodeid(ix  ,iy+1,iz+1);

            % Store the Gmsh element record.
            % Both Gmsh tags are set to the local grain ID.
            element_lines{eid} = ...
                sprintf('%d 5 2 %d %d %d %d %d %d %d %d %d %d\n', ...
                elements_written+1, gID, gID, n1, n2, n3, n4, n5, n6, n7, n8);
            elements_written = elements_written+1;
            eid = eid+1;
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
% global grain IDs. Therefore, the orientation file must use the same compact
% local grain-ID convention as the mesh element tags.
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
% where local_grain_id corresponds to the compact grain ID used in micro.msh.
% -------------------------------------------------------------------------

fprintf('Writing orientation file: %s\n', ori_filename);

% Convert the full orientation vector into an N-by-3 matrix.
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
%   1 = local grain ID
%   2 = Rodrigues r_x
%   3 = Rodrigues r_y
%   4 = Rodrigues r_z
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


%% ---------------------- MICROSTRUCTURE PORE VISUALIZATION ----------------
% Visualize the pore/void phase in the full generated microstructure.
%
% This diagnostic plot renders the pore phase as a blue isosurface and draws
% the bounding box of the full domain. It is used for quick verification of
% pore location and geometry.
%
% Note:
%   This plot uses the full phase array, not the PRISMS subdomain phase_sub.
figure; clf; hold on;

% Pore/void mask in the full generated phase map.
inclu = (phase == (N_grains + 1));

% Render the pore/void isosurface using the full-domain coordinate vectors.
if any(inclu(:))
    fv = isosurface(xvec, yvec, zvec, double(inclu), 0.5);

    p = patch(fv);
    set(p, 'FaceColor', [0 0 1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1.0, ...
        'DiffuseStrength', 0.8, ...
        'SpecularStrength', 0.5, ...
        'SpecularExponent', 10);

    isonormals(xvec, yvec, zvec, double(inclu), p);
end

% Draw full-domain bounding-box edges for spatial reference.
plot3([xvec(1) xvec(end) xvec(end) xvec(1) xvec(1)], ...
    [yvec(1) yvec(1) yvec(end) yvec(end) yvec(1)], ...
    [zvec(1) zvec(1) zvec(1)  zvec(1)  zvec(1)], 'k-', 'LineWidth', 1.0);
plot3([xvec(1) xvec(end) xvec(end) xvec(1) xvec(1)], ...
    [yvec(1) yvec(1) yvec(end) yvec(end) yvec(1)], ...
    [zvec(end) zvec(end) zvec(end) zvec(end) zvec(end)], 'k-', 'LineWidth', 1.0);
for xnow = [xvec(1), xvec(end)]
    for ynow = [yvec(1), yvec(end)]
        plot3([xnow xnow],[ynow ynow],[zvec(1) zvec(end)], 'k-', 'LineWidth', 1.0)
    end
end

view([-42.71 22.97]);
axis equal tight; box on
set(gca, 'Color', [0.9 0.9 0.9]);
camlight('headlight');
lighting gouraud;

hx = xlabel('X_1 (mm)'); set(hx, 'Rotation', 21);
hy = ylabel('X_2 (mm)'); set(hy, 'Rotation', -21);
zlabel('X_3 (mm)');
xticks(linspace(min(xvec), max(xvec), 6));
yticks(linspace(min(yvec), max(yvec), 6));
zticks(linspace(min(zvec), max(zvec), 6));
xlim([min(xvec), max(xvec)]);
ylim([min(yvec), max(yvec)]);
zlim([min(zvec), max(zvec)]);

set(gca, 'FontSize', 15, 'LineWidth', 1.3);
hAxes = gca;
hAxes.BoxStyle    = "full";
hAxes.ClippingStyle = "3dbox";
hold off;
legend('Pore/void phase','Domain boundary','Location','northeast')
