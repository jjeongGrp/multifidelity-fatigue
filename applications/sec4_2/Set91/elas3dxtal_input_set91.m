%% Melt-Pool-Informed Microstructure Generation with XCT Pores
% -------------------------------------------------------------------------
% Generates a 3D voxel-based polycrystalline microstructure within an
% XCT-derived cubic domain.
%
% Workflow:
%   - Builds melt-pool core, boundary, and transition regions
%   - Generates region-specific lognormal grain sizes and aspect ratios
%   - Assigns voxels to anisotropic grains using an ellipsoidal metric
%   - Applies the XCT pore mask
%   - Assigns fiber-textured or random grain orientations
%   - Adds XY air padding and Z caps
%   - Exports phase, orientation, and region data
%
% Input:
%   Cube3p003mm_546vox_547grid_PoreDomain.mat
%     `cube_domain` must contain:
%       dx_mm, cube_mm, x_grid_mm, y_grid_mm, z_grid_mm, and pore_mask.
%
% Outputs:
%   input_structure_poly.h5
%     /pix              - Padded nodal phase IDs, X varying fastest
%     /orientation      - Grain Rodrigues vectors
%     /grain_region     - Region ID for each grain
%     /melt_pool_region - Unpadded voxel-wise region labels
%
%   voxel_volume_summary.txt      - Domain and phase summary
%   xct_poly_params_SS316L.mat    - Selected model parameters
%   xct_poly_params_SS316L_full.mat - Full workspace
%
% Main variables:
%   phase        - Original voxel phase map
%   phase2       - Padded phase map
%   region_label - Melt-pool region map:
%                    1 = core, 2 = boundary, 3 = transition
%   ori_vec      - Grain orientations stored as Rodrigues vectors
%
% Notes:
%   - Grain IDs are 1:N_grains.
%   - XCT pores and external air use phase ID N_grains + 1.
%   - Grain assignment uses parallel block processing and spatial filtering.
%   - Orientation calculations and IPF plots require MTEX.
%   - The Parallel Computing, Statistics and Machine Learning, and Image
%     Processing Toolboxes are required.
% -------------------------------------------------------------------------
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
clear;


%% 1. Microstructure Parameters

% --- Grain size statistics ---
% Mean and standard deviation of the grain diameter in the XY plane.
% Units are millimeters.
P60_d_ebsd = 48.617e-3;
std_d_ebsd = 28.775e-3;


% --- Grain aspect-ratio statistics ---
% 60th percentile aspect ratio and its standard deviation.
P60_ar_ebsd = 2.6377;
std_ar_ebsd = 0.8155;

% Orientation settings
%   'textured' : fiber texture with controlled angular spread
%   'random'   : uniformly random orientations
orientation_type = 'textured';


%% 2. Load XCT Cube Domain and Define Grid

% Load the preprocessed XCT cube domain.
file_name_voxel = 'Cube3p003mm_546vox_547grid_PoreDomain.mat';
cube_file = file_name_voxel;
load(cube_file, 'cube_domain');

% Number of grain orientations to sample for pole-figure visualization.
K_grain = 2000;

% Domain and grid information
% dx is the voxel spacing in millimeters.
dx       = cube_domain.dx_mm;
dx_um    = dx * 1e3;

% cube_len is the side length of the cubic domain in millimeters.
cube_len    = cube_domain.cube_mm;
cube_len_um = cube_len * 1e3;

% Grid-point coordinate vectors.
xvec = cube_domain.x_grid_mm(:);
yvec = cube_domain.y_grid_mm(:);
zvec = cube_domain.z_grid_mm(:);

Nx = numel(xvec);
Ny = numel(yvec);
Nz = numel(zvec);

% XCT-derived pore mask on the voxel/cell grid.
pore_mask = cube_domain.pore_mask;

% Voxel dimensions of the pore/phase grid.
Nyv = size(pore_mask, 1);
Nxv = size(pore_mask, 2);
Nzv = size(pore_mask, 3);

% Total number of voxels and physical domain volume.
num_voxels = Nyv * Nxv * Nzv;
domain_vol  = cube_len^3;

% Voxel-center coordinate vectors.
x_cent = ((1:Nxv) - 0.5) * dx;
y_cent = ((1:Nyv) - 0.5) * dx;
z_cent = ((1:Nzv) - 0.5) * dx;

% 3D coordinate arrays for voxel centers.
[Yc, Xc, Zc] = ndgrid(y_cent, x_cent, z_cent);

%% 2B. Define Melt-Pool-Informed Spatial Regions
% Region IDs:
%   1 = melt-pool core / columnar grain region
%   2 = melt-pool boundary/overlap / fine equiaxed grain region
%   3 = transition/background region
REG_CORE     = uint8(1);
REG_BOUNDARY = uint8(2);
REG_TRANS    = uint8(3);

% -------------------------------------------------------------------------
% LPBF process/melt-pool geometry parameters
% -------------------------------------------------------------------------
mp.layer_thickness = 0.020;
mp.hatch_spacing   = 0.090;

% Melt-pool width
mp.pool_width      = 0.140;

% Melt-pool depth
mp.pool_depth      = 0.080;

% Pool length along scan direction
mp.pool_length     = 0.560;

% Normalized semi-elliptical regions:
mp.core_rho        = 0.50;
mp.boundary_rho    = 0.90;

% Treat overlapping melt pools/remelted areas as fine-equiaxed/boundary.
mp.overlap_as_boundary = false;

% hatch shift between adjacent layers.
mp.alternate_hatch_shift = true;
mp.hatch_shift_fraction  = 0.5;

% -------------------------------------------------------------------------
% Scan rotation strategy
% -------------------------------------------------------------------------
% Defines the layer-wise scan direction used by build_meltpool_region_map.
% Available options for mp.scan_strategy:
%
%   'fixed'
%       Same scan angle for every layer
%       Uses mp.scan_angle_deg
%       mp.scan_strategy = 'fixed';
%       mp.scan_angle_deg = 0;
%
%   'alternate90'
%       Alternates between 0 deg and 90 deg:
%           Layer 0: 0 deg
%           Layer 1: 90 deg
%           Layer 2: 0 deg
%           Layer 3: 90 deg
%           mp.scan_strategy = 'alternate90';
%
%   'rotate67'
%       Rotates scan direction by 67 deg each layer:
%           Layer 0: 0 deg
%           Layer 1: 67 deg
%           Layer 2: 134 deg
%           Layer 3: 21 deg
%           Layer 4: 88 deg
%           ...
%           mp.scan_strategy = 'rotate67';
%
%   'rotate_increment'
%       Rotates by a user-defined angle increment each layer.
%       Uses mp.scan_angle_increment_deg.
%       mp.scan_strategy = 'rotate_increment';
%       mp.scan_angle_increment_deg = 45;
%
%   'custom'
%       Uses a user-defined repeating list of scan angles.
%       Uses mp.scan_angles_deg.
%       mp.scan_strategy = 'custom';
%       mp.scan_angles_deg = [0 67 134 21 88 155 42 109 176 63];
%
% -------------------------------------------------------------------------
mp.scan_strategy = 'rotate67';

% scan-grid origin
mp.scan_origin = [0, 0];

% Build melt-pool region map.
[region_label, rho_min_map, overlap_count] = build_meltpool_region_map( ...
    y_cent, x_cent, z_cent, mp, REG_CORE, REG_BOUNDARY, REG_TRANS);


%% 3. Generate Region-Specific Grain Sizes, Aspect Ratios, and Seeds
%--------------------------------------------------------------------------
% Region-specific microstructural statistics.
% Grain-bearing solid regions:
%   1. REG_CORE     : melt-pool core / columnar grains
%   2. REG_BOUNDARY : melt-pool boundary/overlap / fine equiaxed grains
%   3. REG_TRANS    : transition/background grains
%--------------------------------------------------------------------------
% Solid material mask
solid_mask = ~pore_mask;

n_solid_voxels = nnz(solid_mask);

if n_solid_voxels == 0
    error('No solid voxels found. Check pore_mask.');
end

% Compute solid-region volume fractions
frac_core  = nnz(region_label == REG_CORE     & solid_mask) / n_solid_voxels;
frac_bnd   = nnz(region_label == REG_BOUNDARY & solid_mask) / n_solid_voxels;
frac_trans = nnz(region_label == REG_TRANS    & solid_mask) / n_solid_voxels;

w_vol = [frac_core, frac_bnd, frac_trans];
w_vol = w_vol / sum(w_vol);

% Region-specific relative statistics
sigma_d  = [0.50, 0.45, 0.35];
rel_d    = [2.50, 0.70, 1.00];

sigma_ar = [0.22, 0.25, 0.22];
rel_ar   = [1.10, 0.80, 1.00];

% calibration weights
w_cal = w_vol ./ rel_d.^2;
w_cal = w_cal / sum(w_cal);

% Mixture distribution function handles
mixture_pdf = @(x, params, weights) ...
    sum(cell2mat(arrayfun(@(i) weights(i) * ...
    lognpdf(x, params(i,1), params(i,2)), ...
    1:size(params,1), 'UniformOutput', false)'), 1);

mixture_cdf = @(x, params, weights) ...
    sum(cell2mat(arrayfun(@(i) weights(i) * ...
    logncdf(x, params(i,1), params(i,2)), ...
    1:size(params,1), 'UniformOutput', false)'), 1);

mixture_percentile = @(p, params, weights) ...
    fzero(@(x) mixture_cdf(x, params, weights) - p, ...
    [1e-9, 1e3]);

% Grain size mixture calibration
params_d_base = [
    log(rel_d(1)) - 0.5*sigma_d(1)^2, sigma_d(1);
    log(rel_d(2)) - 0.5*sigma_d(2)^2, sigma_d(2);
    log(rel_d(3)) - 0.5*sigma_d(3)^2, sigma_d(3)
    ];

P60_base_d = mixture_percentile(0.60, params_d_base, w_cal);
target_P60_d_um = P60_d_ebsd * 1e3;
k_d_cal_um = target_P60_d_um / P60_base_d;

mean_d_regions_um = k_d_cal_um * rel_d;
mean_d_regions    = mean_d_regions_um * 1e-3;

mu_d_regions_um = log(mean_d_regions_um) - 0.5*sigma_d.^2;
params_d = [mu_d_regions_um(:), sigma_d(:)];
P60_mix_check = mixture_percentile(0.60, params_d, w_cal);

% Exact mixture mean/std in µm.
mean_components_d_um = exp(params_d(:,1) + 0.5*params_d(:,2).^2);

var_components_d_um = ...
    (exp(params_d(:,2).^2) - 1) .* ...
    exp(2*params_d(:,1) + params_d(:,2).^2);

mean_mix = sum(w_cal(:) .* mean_components_d_um(:));
second_moment_d = sum(w_cal(:) .* ...
    (var_components_d_um(:) + mean_components_d_um(:).^2));
std_mix = sqrt(second_moment_d - mean_mix^2);

% --- Aspect ratio mixture calibration ---
w_ar_cal = w_cal;

params_ar_base = [
    log(rel_ar(1)) - 0.5*sigma_ar(1)^2, sigma_ar(1);
    log(rel_ar(2)) - 0.5*sigma_ar(2)^2, sigma_ar(2);
    log(rel_ar(3)) - 0.5*sigma_ar(3)^2, sigma_ar(3)
    ];

P60_base_ar = mixture_percentile(0.60, params_ar_base, w_ar_cal);
k_ar_cal = P60_ar_ebsd / P60_base_ar;

mean_ar_regions = k_ar_cal * rel_ar;
mu_ar_regions   = log(mean_ar_regions) - 0.5*sigma_ar.^2;

params_ar = [mu_ar_regions(:), sigma_ar(:)];
P60_mix_check_ar = mixture_percentile(0.60, params_ar, w_ar_cal);

mean_components_ar = exp(params_ar(:,1) + 0.5*params_ar(:,2).^2);
var_components_ar = ...
    (exp(params_ar(:,2).^2) - 1) .* ...
    exp(2*params_ar(:,1) + params_ar(:,2).^2);

mean_mix_ar = sum(w_ar_cal(:) .* mean_components_ar(:));
second_moment_ar = sum(w_ar_cal(:) .* ...
    (var_components_ar(:) + mean_components_ar(:).^2));
std_mix_ar = sqrt(second_moment_ar - mean_mix_ar^2);
std_d_regions  = mean_d_regions  .* sqrt(exp(sigma_d.^2)  - 1);
std_ar_regions = mean_ar_regions .* sqrt(exp(sigma_ar.^2) - 1);

% --- Populate stats structure ---
stats = struct();

% Region 1: melt-pool core / columnar grains
stats(double(REG_CORE)).name        = 'melt_pool_core_columnar';

stats(double(REG_CORE)).mean_d      = mean_d_regions(1);
stats(double(REG_CORE)).std_d       = std_d_regions(1);

stats(double(REG_CORE)).mean_aspect = mean_ar_regions(1);
stats(double(REG_CORE)).std_aspect  = std_ar_regions(1);

stats(double(REG_CORE)).texture     = 'strong_fiber';
stats(double(REG_CORE)).sigma_spread = 10;
stats(double(REG_CORE)).random_orientation_fraction = 0.02;

% Region 2: melt-pool boundary / overlap / fine equiaxed grains
stats(double(REG_BOUNDARY)).name        = 'melt_pool_boundary_fine_equiaxed';

stats(double(REG_BOUNDARY)).mean_d      = mean_d_regions(2);
stats(double(REG_BOUNDARY)).std_d       = std_d_regions(2);

stats(double(REG_BOUNDARY)).mean_aspect = mean_ar_regions(2);
stats(double(REG_BOUNDARY)).std_aspect  = std_ar_regions(2);

stats(double(REG_BOUNDARY)).texture     = 'weak_fiber';
stats(double(REG_BOUNDARY)).sigma_spread = 30;
stats(double(REG_BOUNDARY)).random_orientation_fraction = 0.30;

% Region 3: transition/background
stats(double(REG_TRANS)).name        = 'transition';

stats(double(REG_TRANS)).mean_d      = mean_d_regions(3);
stats(double(REG_TRANS)).std_d       = std_d_regions(3);

stats(double(REG_TRANS)).mean_aspect = mean_ar_regions(3);
stats(double(REG_TRANS)).std_aspect  = std_ar_regions(3);

stats(double(REG_TRANS)).texture     = 'moderate_fiber';
stats(double(REG_TRANS)).sigma_spread = 18;
stats(double(REG_TRANS)).random_orientation_fraction = 0.15;

% Generate region-specific grain seeds, sizes, and aspect ratios
grain_seeds     = [];
grain_diameters = [];
aspect_ratios   = [];
grain_regions   = uint8([]);
active_grain_regs = uint8([REG_CORE, REG_BOUNDARY, REG_TRANS]);

for reg = active_grain_regs

    reg_id = double(reg);
    region_mask = (region_label == reg) & solid_mask;
    region_voxels = nnz(region_mask);
    region_vol    = region_voxels * dx^3;
    if region_voxels == 0
        warning('Region %d has no non-pore voxels. Skipping.', reg_id);
        continue;
    end

    mean_d_reg = stats(reg_id).mean_d;
    std_d_reg  = stats(reg_id).std_d;
    mean_aspect_reg = stats(reg_id).mean_aspect;
    std_aspect_reg  = stats(reg_id).std_aspect;
    mean_grain_vol_reg = (4/3) * pi * (mean_d_reg/2)^2 * ...
        (mean_aspect_reg * mean_d_reg / 2);
    N_reg = max(1, ceil(region_vol / mean_grain_vol_reg));

    MAX_GRAINS_PER_REGION = 5000000;
    if N_reg > MAX_GRAINS_PER_REGION
        warning('Region %d: N_reg=%d exceeds cap, truncating to %d.', ...
            reg_id, N_reg, MAX_GRAINS_PER_REGION);
        N_reg = MAX_GRAINS_PER_REGION;
    end

    % Sample lognormal transverse grain diameters.
    sigma_log_d_reg = sqrt(log((std_d_reg/mean_d_reg)^2 + 1));
    mu_log_d_reg    = log(mean_d_reg) - 0.5*sigma_log_d_reg^2;
    d_reg = lognrnd(mu_log_d_reg, sigma_log_d_reg, N_reg, 1);

    % Sample lognormal aspect ratios.
    sigma_log_a_reg = sqrt(log((std_aspect_reg/mean_aspect_reg)^2 + 1));
    mu_log_a_reg    = log(mean_aspect_reg) - 0.5*sigma_log_a_reg^2;
    a_reg = lognrnd(mu_log_a_reg, sigma_log_a_reg, N_reg, 1);

    % Sample seed locations from solid voxels in this region.
    idx_reg = find(region_mask);
    if isempty(idx_reg)
        warning('Region %d has no valid non-pore voxels for seed sampling. Skipping.', reg_id);
        continue;
    end

    sampled_idx = idx_reg(randi(numel(idx_reg), N_reg, 1));
    sampled_idx = sampled_idx(:);
    [iy_s, ix_s, iz_s] = ind2sub(size(region_label), sampled_idx);
    seeds_reg = [x_cent(ix_s(:)), y_cent(iy_s(:)), z_cent(iz_s(:))];
    seeds_reg = reshape(seeds_reg, [], 3);

    % Append to global grain arrays.
    grain_seeds     = [grain_seeds; seeds_reg];
    grain_diameters = [grain_diameters; d_reg];
    aspect_ratios   = [aspect_ratios; a_reg];
    grain_regions   = [grain_regions; repmat(reg, N_reg, 1)];

end

N_grains = size(grain_seeds, 1);
if N_grains == 0
    error('No grains were generated. Check region_label, pore_mask, and region parameters.');
end
r_xy = grain_diameters / 2;
r_z = aspect_ratios .* r_xy;


%% 4. Anisotropic Grain Assignment Using Parallel Block Processing
fprintf('\n=== PARALLEL + SPATIAL FILTERING MICROSTRUCTURE GENERATION ===\n');
main_timer = tic;

% Start a parallel pool if one is not already active
if isempty(gcp('nocreate'))
    fprintf('Starting parallel pool...\n');
    pool_timer = tic;
    num_workers = feature('numcores');
    parpool('local', num_workers);
    fprintf('Parallel pool started in %.1fs with %d workers\n', ...
        toc(pool_timer), gcp().NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', gcp().NumWorkers);
end


% Preallocate voxel-wise phase array
phase = zeros(Nyv, Nxv, Nzv, 'uint32');
solid_mask = ~pore_mask;

% Region-aware grain assignment settings.
use_region_penalty = true;
% Penalty added to normalized ellipsoidal distance
region_mismatch_penalty = single(1.0);

% Local copies used inside parfor
region_label_local  = region_label;
grain_regions_local = grain_regions;
pore_mask_local     = pore_mask;
solid_mask_local    = solid_mask;

% Precompute inverse squared radii for faster distance calculations
inv_r_xy2 = single(1 ./ (r_xy.^2));
inv_r_z2  = single(1 ./ (r_z.^2));

% Store grain seeds in single precision to reduce memory use.
grain_seeds_single = single(grain_seeds);

% Build spatial acceleration structure
fprintf('\nBuilding spatial acceleration structure...\n');
accel_timer = tic;

grid_res = 15;
x_edges = linspace(min(xvec), max(xvec), grid_res+1);
y_edges = linspace(min(yvec), max(yvec), grid_res+1);
z_edges = linspace(min(zvec), max(zvec), grid_res+1);
max_influence = max([prctile(r_xy, 99); prctile(r_z, 99)]) * 3.0;

total_cells = grid_res^3;
spatial_grid_temp = cell(total_cells, 1);

parfor linear_idx = 1:total_cells

    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);

    cell_xmin = x_edges(i);
    cell_xmax = x_edges(i+1);
    cell_ymin = y_edges(j);
    cell_ymax = y_edges(j+1);
    cell_zmin = z_edges(k);
    cell_zmax = z_edges(k+1);

    relevant_grains = find( ...
        grain_seeds_single(:,1) >= (cell_xmin - max_influence) & ...
        grain_seeds_single(:,1) <= (cell_xmax + max_influence) & ...
        grain_seeds_single(:,2) >= (cell_ymin - max_influence) & ...
        grain_seeds_single(:,2) <= (cell_ymax + max_influence) & ...
        grain_seeds_single(:,3) >= (cell_zmin - max_influence) & ...
        grain_seeds_single(:,3) <= (cell_zmax + max_influence));

    spatial_grid_temp{linear_idx} = relevant_grains;

end

% Convert flat parfor output into 3D cell array.
spatial_grid = cell(grid_res, grid_res, grid_res);
for linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);
    spatial_grid{i,j,k} = spatial_grid_temp{linear_idx};
end

clear spatial_grid_temp;

% Block-processing setup
block_size = 50;
n_blocks_x = ceil(Nxv / block_size);
n_blocks_y = ceil(Nyv / block_size);
n_blocks_z = ceil(Nzv / block_size);
total_blocks = n_blocks_x * n_blocks_y * n_blocks_z;

num_workers = gcp().NumWorkers;
blocks_per_batch = max(1, floor(total_blocks / (num_workers * 4)));
n_batches = ceil(total_blocks / blocks_per_batch);

% Track progress and performance statistics.
completed_blocks = 0;
total_grains_processed = 0;
blocks_with_grains = 0;
batch_times = zeros(1, n_batches);

% Process the domain one batch at a time.
fprintf('\nStarting parallel processing with real-time progress:\n');
process_timer = tic;

% Main batch loop
for batch_idx = 1:n_batches

    batch_timer = tic;
    batch_start = (batch_idx - 1) * blocks_per_batch + 1;
    batch_end   = min(batch_idx * blocks_per_batch, total_blocks);
    current_batch_size = batch_end - batch_start + 1;

    % Create list of block indices for current batch
    block_list_batch = zeros(current_batch_size, 4);

    for local_idx = 1:current_batch_size
        block_num = batch_start + local_idx - 1;
        [bx, by, bz] = ind2sub([n_blocks_x, n_blocks_y, n_blocks_z], block_num);
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
    x_cent_local = x_cent;
    y_cent_local = y_cent;
    z_cent_local = z_cent;

    % Domain limits from node/grid coordinates.
    x_domain_min = min(xvec);
    x_domain_max = max(xvec);
    y_domain_min = min(yvec);
    y_domain_max = max(yvec);
    z_domain_min = min(zvec);
    z_domain_max = max(zvec);

    dx_local = dx;

    % Process current batch in parallel
    parfor local_idx = 1:current_batch_size

        bx = block_list_batch(local_idx, 1);
        by = block_list_batch(local_idx, 2);
        bz = block_list_batch(local_idx, 3);

        ix1 = (bx-1)*block_size + 1;
        ix2 = min(bx*block_size, Nxv);
        iy1 = (by-1)*block_size + 1;
        iy2 = min(by*block_size, Nyv);
        iz1 = (bz-1)*block_size + 1;
        iz2 = min(bz*block_size, Nzv);

        if numel(xvec_local) == Nxv + 1
            block_xmin = xvec_local(ix1);
            block_xmax = xvec_local(ix2 + 1);
        else
            block_xmin = x_cent_local(ix1) - dx_local/2;
            block_xmax = x_cent_local(ix2) + dx_local/2;
        end

        if numel(yvec_local) == Nyv + 1
            block_ymin = yvec_local(iy1);
            block_ymax = yvec_local(iy2 + 1);
        else
            block_ymin = y_cent_local(iy1) - dx_local/2;
            block_ymax = y_cent_local(iy2) + dx_local/2;
        end

        if numel(zvec_local) == Nzv + 1
            block_zmin = zvec_local(iz1);
            block_zmax = zvec_local(iz2 + 1);
        else
            block_zmin = z_cent_local(iz1) - dx_local/2;
            block_zmax = z_cent_local(iz2) + dx_local/2;
        end

        % Compute spatial-grid cell range overlapped by this block.
        gx1 = floor((block_xmin - x_domain_min) / ...
            (x_domain_max - x_domain_min) * grid_res) + 1;
        gx2 = ceil((block_xmax - x_domain_min) / ...
            (x_domain_max - x_domain_min) * grid_res);
        gy1 = floor((block_ymin - y_domain_min) / ...
            (y_domain_max - y_domain_min) * grid_res) + 1;
        gy2 = ceil((block_ymax - y_domain_min) / ...
            (y_domain_max - y_domain_min) * grid_res);
        gz1 = floor((block_zmin - z_domain_min) / ...
            (z_domain_max - z_domain_min) * grid_res) + 1;
        gz2 = ceil((block_zmax - z_domain_min) / ...
            (z_domain_max - z_domain_min) * grid_res);
        gx1 = max(1, min(grid_res, gx1));
        gx2 = max(1, min(grid_res, gx2));
        gy1 = max(1, min(grid_res, gy1));
        gy2 = max(1, min(grid_res, gy2));
        gz1 = max(1, min(grid_res, gz1));
        gz2 = max(1, min(grid_res, gz2));
        gx2 = max(gx1, gx2);
        gy2 = max(gy1, gy2);
        gz2 = max(gz1, gz2);

        n_cells_local = (gx2-gx1+1) * (gy2-gy1+1) * (gz2-gz1+1);
        cand_cells = cell(n_cells_local, 1);

        ccell = 0;
        for gx = gx1:gx2
            for gy = gy1:gy2
                for gz = gz1:gz2
                    ccell = ccell + 1;
                    cand_cells{ccell} = spatial_grid{gx,gy,gz};
                end
            end
        end

        cand_cells = cand_cells(~cellfun('isempty', cand_cells));
        if isempty(cand_cells)
            relevant_grains = [];
        else
            relevant_grains = unique(vertcat(cand_cells{:}));
        end

        block_shape = [iy2-iy1+1, ix2-ix1+1, iz2-iz1+1];
        block_phase = zeros(block_shape, 'uint32');

        Pb = pore_mask_local(iy1:iy2, ix1:ix2, iz1:iz2);
        solid_block_mask = ~Pb;
        if ~any(solid_block_mask(:))
            batch_results{local_idx} = block_phase;
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = 0;
            continue;
        end

        % Assign solid voxels to nearest relevant grain
        if ~isempty(relevant_grains)
            [Yb, Xb, Zb] = ndgrid( ...
                single(y_cent_local(iy1:iy2)), ...
                single(x_cent_local(ix1:ix2)), ...
                single(z_cent_local(iz1:iz2)));

            % Melt-pool region labels for the current block.
            Rb = region_label_local(iy1:iy2, ix1:ix2, iz1:iz2);

            % Initialize local block phase map and minimum-distance array.
            min_dist2 = inf(size(Yb), 'single');

            for k = relevant_grains'
                dX = Xb - grain_seeds_single(k,1);
                dY = Yb - grain_seeds_single(k,2);
                dZ = Zb - grain_seeds_single(k,3);

                current_dist2 = ...
                    (dX.^2) * inv_r_xy2(k) + ...
                    (dY.^2) * inv_r_xy2(k) + ...
                    (dZ.^2) * inv_r_z2(k);

                % Region-aware penalty
                if use_region_penalty
                    mismatch = Rb ~= grain_regions_local(k);
                    current_dist2 = current_dist2 + ...
                        region_mismatch_penalty * single(mismatch);
                end

                update_mask = solid_block_mask & current_dist2 < min_dist2;
                block_phase(update_mask) = uint32(k);
                min_dist2(update_mask) = current_dist2(update_mask);
            end

            batch_results{local_idx} = block_phase;
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = length(relevant_grains);

        else

            % Solid voxels may remain unassigned
            batch_results{local_idx} = block_phase;
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
            ix1 = indices(1);
            ix2 = indices(2);
            iy1 = indices(3);
            iy2 = indices(4);
            iz1 = indices(5);
            iz2 = indices(6);
            phase(iy1:iy2, ix1:ix2, iz1:iz2) = ...
                cast(batch_results{local_idx}, 'uint32');

            if batch_grain_counts(local_idx) > 0
                batch_blocks_with_grains = batch_blocks_with_grains + 1;
            end

        end
        % Sum up grains processed in batch
        batch_grains_processed = batch_grains_processed + batch_grain_counts(local_idx);

    end

    % Progress reporting
    completed_blocks = completed_blocks + current_batch_size;
    total_grains_processed = total_grains_processed + batch_grains_processed;
    blocks_with_grains = blocks_with_grains + batch_blocks_with_grains;

    % Record batch timing
    batch_time = toc(batch_timer);
    batch_times(batch_idx) = batch_time;

    % Calculate progress statistics
    progress_pct = 100 * completed_blocks / total_blocks;
    elapsed_total = toc(main_timer);

    avg_grains_per_block = batch_grains_processed / current_batch_size;

    % Estimate ETA for remaining batches
    if batch_idx > 1
        avg_batch_time = mean(batch_times(1:batch_idx));
        remaining_batches = n_batches - batch_idx;
        eta_seconds = remaining_batches * avg_batch_time;
    else
        eta_seconds = batch_time * (n_batches - 1);
    end


    % Print detailed progress
    fprintf(['Batch %d/%d: Blocks %d-%d (%.1f%%) | ', ...
        '%.1f grains/block avg | Elapsed: %.1fs | ETA: %.1fs (%.1f min)\n'], ...
        batch_idx, n_batches, batch_start, batch_end, progress_pct, ...
        avg_grains_per_block, elapsed_total, eta_seconds, eta_seconds/60);

    if mod(batch_idx, max(1, floor(n_batches/5))) == 0 || batch_idx == n_batches

        overall_avg_grains = total_grains_processed / completed_blocks;
        processing_rate = completed_blocks / elapsed_total;
        filtering_efficiency = 1 - (overall_avg_grains / N_grains);

        fprintf(['  --> Detailed Stats: %.1f grains/block overall | ', ...
            '%.1f blocks/sec | %.1f%% grains filtered out\n'], ...
            overall_avg_grains, processing_rate, 100 * filtering_efficiency);

        if batch_idx < n_batches

            fprintf(['  --> Memory: %d/%d blocks processed | ', ...
                '%d/%d non-empty | %.1f%% efficiency\n'], ...
                completed_blocks, total_blocks, blocks_with_grains, ...
                completed_blocks, 100*blocks_with_grains/completed_blocks);

        end

    end

end

% Assign XCT pore phase after solid grain assignment
pore_phase_id = uint32(N_grains + 1);
phase(pore_mask) = pore_phase_id;

fprintf('\nParallel + spatial filtering microstructure generation completed!\n');


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

% Region-specific texture spread assigned to each grain.
sigma_spread_by_grain = zeros(N_grains, 1);
random_orientation_fraction_by_grain = zeros(N_grains, 1);

for k = 1:N_grains
    reg_id = double(grain_regions(k));
    sigma_spread_by_grain(k) = stats(reg_id).sigma_spread;
    random_orientation_fraction_by_grain(k) = ...
        stats(reg_id).random_orientation_fraction;
end

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

        % Some grains, especially in melt-pool boundary/overlap regions, may be
        % assigned random orientations to represent weaker texture.
        if rand() < random_orientation_fraction_by_grain(k)

            % --- Inline: generate uniformly random orientation as Rodrigues vector ---
            % Uses the Shoemake (1992) method to produce a uniform random quaternion,
            % then converts to a Rodrigues vector.
            u1_r = rand();
            u2_r = rand();
            u3_r = rand();

            q_w_r = sqrt(1 - u1_r) * sin(2*pi*u2_r);
            q_x_r = sqrt(1 - u1_r) * cos(2*pi*u2_r);
            q_y_r = sqrt(u1_r)     * sin(2*pi*u3_r);
            q_z_r = sqrt(u1_r)     * cos(2*pi*u3_r);

            ROD_MAX = tan(89.9 * pi/180);
            if abs(q_w_r) < 1e-8
                rod_vec = [q_x_r; q_y_r; q_z_r] / ...
                    max(abs([q_x_r q_y_r q_z_r])) * ROD_MAX;
            else
                rod_vec = [q_x_r; q_y_r; q_z_r] / q_w_r;
            end
            % -------------------------------------------------------------------------

            idx = 3*(k-1) + 1;
            ori_vec(idx:idx+2) = rod_vec(:);

            continue;
        end

        % Apply a random spin about sample Z to create the fiber texture.
        phi    = 2*pi*rand();
        R_spin = rotation.byAxisAngle(zvector, phi);

        accepted = false;
        R_wobble = rotation.id;

        for attempt = 1:max_attempts

            % Apply a small random wobble about an axis in the sample XY plane.
            % The wobble controls the angular spread around the ideal fiber.
            sigma_k = sigma_spread_by_grain(k);
            wobble_ang = (sigma_k * max(-3, min(3, randn()))) * degree;

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
            h_ipfz = normalize(R_total_candidate \ zvector);

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
        ROD_MAX = tan(89.9 * pi/180);
        if abs(q_w) < 1e-8
            rod_vec = [q_x q_y q_z] / max(abs([q_x q_y q_z])) * ROD_MAX;
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

% Write grain region labels.
h5create(h5filename, '/grain_region', size(grain_regions), 'Datatype', 'uint8');
h5write(h5filename, '/grain_region', grain_regions);
h5writeatt(h5filename, '/grain_region', 'description', ...
    '1=core/columnar, 2=boundary/equiaxed/overlap, 3=transition');

% Write voxel-wise melt-pool region labels.
% Convert MATLAB array order [Y, X, Z] to output order [X, Y, Z].
region_xyz  = permute(region_label, [2 1 3]);
region_flat = uint8(region_xyz(:));

h5create(h5filename, '/melt_pool_region', [length(region_flat) 1], ...
    'Datatype', 'uint8');
h5write(h5filename, '/melt_pool_region', region_flat);
h5writeatt(h5filename, '/melt_pool_region', ...
    'dimensions', int32([Nxv Nyv Nzv]));
h5writeatt(h5filename, '/melt_pool_region', 'description', ...
    '1=core/columnar, 2=boundary/equiaxed/overlap, 3=transition');

disp(['Saved HDF5 input structure to ', h5filename]);

n_phases = max(phase(:));
% Save basic parameters to a MATLAB .mat file
save('xct_poly_params_SS316L.mat','xvec','yvec','zvec','Nx','Ny','Nz', ...
    'cube_len','N_grains','n_phases', ...
    'region_label','grain_regions','stats','mp', ...
    'grain_diameters','aspect_ratios');
% Save all current workspace variables to a .mat file (version 7.3 for large arrays)
save('xct_poly_params_SS316L_full.mat','-v7.3')

%% ---- END OF SCRIPT ----

function [region_label, rho_min_map, overlap_count] = build_meltpool_region_map( ...
    y_cent, x_cent, z_cent, mp, REG_CORE, REG_BOUNDARY, REG_TRANS)

% BUILD_MELTPOOL_REGION_MAP
Ny = numel(y_cent);
Nx = numel(x_cent);
Nz = numel(z_cent);

% Melt-pool geometry
half_width = mp.pool_width / 2;
depth      = mp.pool_depth;

if isfield(mp, 'pool_length')
    half_length = mp.pool_length / 2;
else
    half_length = 3.0 * half_width;
    fprintf('  mp.pool_length not set — using default 3x pool_width = %.3f mm\n', ...
        2 * half_length);
end

% Search factor
rho_search = max([1.0, mp.core_rho, mp.boundary_rho]);
search_factor = sqrt(rho_search);

% --- Domain limits ---
xmin = min(x_cent);
xmax = max(x_cent);
ymin = min(y_cent);
ymax = max(y_cent);
zmin = min(z_cent);
zmax = max(z_cent);

n_layers = ceil(zmax / mp.layer_thickness) + 2;

if isfield(mp, 'pool_x_spacing')
    pool_u_spacing = mp.pool_x_spacing;
else
    pool_u_spacing = half_length;
end

if isfield(mp, 'scan_origin')
    x_origin = mp.scan_origin(1);
    y_origin = mp.scan_origin(2);
else
    x_origin = 0;
    y_origin = 0;
end

if ~isfield(mp, 'overlap_as_boundary')
    mp.overlap_as_boundary = false;
end

if ~isfield(mp, 'alternate_hatch_shift')
    mp.alternate_hatch_shift = false;
end

if ~isfield(mp, 'hatch_shift_fraction')
    mp.hatch_shift_fraction = 0.5;
end


rho_min_map   = inf(Ny, Nx, Nz, 'single');
overlap_count = zeros(Ny, Nx, Nz, 'uint8');
x_corners = [xmin xmax xmax xmin] - x_origin;
y_corners = [ymin ymin ymax ymax] - y_origin;

for ell = 0:n_layers
    z_top = ell * mp.layer_thickness;
    if z_top - search_factor * depth > zmax || z_top < zmin
        continue;
    end

    theta_deg = get_layer_scan_angle_deg(mp, ell);
    ct = cosd(theta_deg);
    st = sind(theta_deg);
    u_corners =  x_corners * ct + y_corners * st;
    v_corners = -x_corners * st + y_corners * ct;
    umin = min(u_corners);
    umax = max(u_corners);
    vmin = min(v_corners);
    vmax = max(v_corners);

    if isfield(mp, 'alternate_hatch_shift') && mp.alternate_hatch_shift && mod(ell, 2) == 1
        hatch_offset = mp.hatch_shift_fraction * mp.hatch_spacing;
    else
        hatch_offset = 0;
    end

    if mod(ell, 2) == 1
        u_scan_offset = pool_u_spacing / 2;
    else
        u_scan_offset = 0;
    end

    h_min = floor((vmin - search_factor * half_width  - hatch_offset) / mp.hatch_spacing) - 2;
    h_max = ceil( (vmax + search_factor * half_width  - hatch_offset) / mp.hatch_spacing) + 2;
    p_min = floor((umin - search_factor * half_length - u_scan_offset) / pool_u_spacing) - 2;
    p_max = ceil( (umax + search_factor * half_length + u_scan_offset) / pool_u_spacing) + 2;

    for h = h_min:h_max
        v0 = h * mp.hatch_spacing + hatch_offset;
        for pu = p_min:p_max
            u0 = pu * pool_u_spacing + u_scan_offset;
            x0 = x_origin + u0 * ct - v0 * st;
            y0 = y_origin + u0 * st + v0 * ct;
            x_radius = search_factor * (abs(half_length * ct) + abs(half_width * st));
            y_radius = search_factor * (abs(half_length * st) + abs(half_width * ct));
            ix_range = find(x_cent >= x0 - x_radius & x_cent <= x0 + x_radius);
            iy_range = find(y_cent >= y0 - y_radius & y_cent <= y0 + y_radius);
            iz_range = find(z_cent >= z_top - search_factor * depth & z_cent <= z_top);

            if isempty(ix_range) || isempty(iy_range) || isempty(iz_range)
                continue;
            end

            [Y2, X2] = ndgrid(y_cent(iy_range), x_cent(ix_range));
            Xr = X2 - x_origin;
            Yr = Y2 - y_origin;
            U2 =  Xr * ct + Yr * st;
            V2 = -Xr * st + Yr * ct;
            du2 = ((U2 - u0) / half_length).^2;
            dv2 = ((V2 - v0) / half_width ).^2;

            rho_xy = du2 + dv2;
            dz = (z_cent(iz_range) - z_top) / depth;
            dz2 = reshape(dz.^2, [1, 1, numel(iz_range)]);

            rho_sub = single(reshape(rho_xy, [numel(iy_range), numel(ix_range), 1]) + dz2);
            inside = rho_sub <= 1.0;
            rho_prev = rho_min_map(iy_range, ix_range, iz_range);
            rho_min_map(iy_range, ix_range, iz_range) = min(rho_prev, rho_sub);
            overlap_count(iy_range, ix_range, iz_range) = ...
                overlap_count(iy_range, ix_range, iz_range) + uint8(inside);

        end
    end
end

region_label = REG_TRANS * ones(Ny, Nx, Nz, 'uint8');
region_label(rho_min_map <= mp.core_rho) = REG_CORE;
boundary_mask = rho_min_map > mp.core_rho & rho_min_map <= mp.boundary_rho;
region_label(boundary_mask) = REG_BOUNDARY;
if mp.overlap_as_boundary
    region_label(overlap_count >= 2) = REG_BOUNDARY;
end

end


function theta_deg = get_layer_scan_angle_deg(mp, ell)
% GET_LAYER_SCAN_ANGLE_DEG Returns scan angle for layer ell in degrees.

if ~isfield(mp, 'scan_strategy')
    theta_deg = 0;
    return;
end

strategy = lower(mp.scan_strategy);
switch strategy
    case {'fixed', 'uni', 'unidirectional', 'x', 'fixed0'}
        if isfield(mp, 'scan_angle_deg')
            theta_deg = mp.scan_angle_deg;
        else
            theta_deg = 0;
        end

    case {'fixed90', 'y'}
        theta_deg = 90;

    case {'alternate90', 'alternating90', '0_90', '90'}
        theta_deg = 90 * mod(ell, 2);

    case {'rotate67', 'rotating67', 'increment67'}
        theta_deg = ell * 67;

    case {'rotate_increment', 'increment', 'rotating'}
        if ~isfield(mp, 'scan_angle_increment_deg')
            error('mp.scan_angle_increment_deg must be set for rotate_increment strategy.');
        end
        theta_deg = ell * mp.scan_angle_increment_deg;

    case {'custom', 'custom_angles', 'angle_list'}
        if ~isfield(mp, 'scan_angles_deg')
            error('mp.scan_angles_deg must be set for custom scan strategy.');
        end
        idx = mod(ell, numel(mp.scan_angles_deg)) + 1;
        theta_deg = mp.scan_angles_deg(idx);

    otherwise
        error('Unknown mp.scan_strategy: %s', mp.scan_strategy);
end

theta_deg = mod(theta_deg, 180);

end