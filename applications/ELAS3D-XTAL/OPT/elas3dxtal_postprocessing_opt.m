%% elas3d-xtal Full-Field Stress Postprocessing and Hotspot Analysis
% -------------------------------------------------------------------------
% This script postprocesses full-field elastic stress results from the
% elas3d-xtal finite element solver for a padded polycrystalline
% microstructure with XCT-derived pores.
%
% The script reads the solver output HDF5 file containing von Mises stress,
% reconstructs the padded simulation grid, identifies solid and pore/air
% regions from the input phase map, and performs fatigue-hotspot analysis near
% internal pore surfaces. Hotspots are ranked using stress-, volume-,
% defect-size-, and shape-aware severity metrics.
%
% Main workflow:
%   1. Load original microstructure/grid parameters.
%   2. Reconstruct the padded nodal coordinate system used by the Fortran solver.
%   3. Read von Mises stress from fullfield_poly.h5.
%   4. Read phase labels from input_structure_poly.h5.
%   5. Identify solid, pore/air, core, and analysis-region masks.
%   6. Apply solid-aware critical-distance stress averaging.
%   7. Build a pore-surface band around internal pores.
%   8. Detect connected high-stress hotspot islands.
%   9. Associate each hotspot island with the nearest/touching pore component.
%  10. Compute defect-aware hotspot severity scores.
%  11. Rank the most critical hotspots.
%  12. Estimate crystal-plasticity subdomain box sizes around hotspots.
%  13. Write a text log and save the full MATLAB workspace.
%  14. Generate 2D slice plots and a 3D hotspot visualization.
%
% Inputs:
%   - xct_poly_params_SS316L.mat
%       MAT-file containing original cube/grid parameters from microstructure
%       generation, including:
%         * xvec, yvec, zvec : original unpadded coordinate vectors, mm
%         * Nx, Ny, Nz       : original grid dimensions
%         * cube_len         : original cube side length, mm
%
%   - input_structure_poly.h5
%       HDF5 input file used by the Fortran solver. This script reads:
%         * /pix
%             Phase-label array over the padded nodal simulation grid.
%             Solid grain phases are labeled 1:(N_grains).
%             The pore/air phase is assumed to be the maximum phase ID.
%
%   - fullfield_poly.h5
%       HDF5 output file written by the Fortran solver. This script reads:
%         * /vm
%             von Mises stress field over the padded simulation grid.
%             Values are assumed to be in Pa and are converted to MPa.
%
% Main Outputs:
%   - hotspot_analysis_log.txt
%       Text log summarizing critical-distance averaging, pore-band thresholding,
%       hotspot island properties, severity ranking, and crystal-plasticity
%       subdomain box sizing.
%
%   - xct_poly_post_SS316L_full.mat
%       MATLAB workspace file containing all computed masks, fields, hotspot
%       metrics, rankings, and plotting variables.
%
%   - MATLAB figures
%       Generated figures include:
%         * zoomed YZ slices through selected hotspots
%         * zoomed XY slices through selected hotspots
%         * 3D surface scatter plot with hotspot markers
%
% Important analysis definitions:
%   - Pore/air phase:
%       pore_phase_id = max(unique(/pix))
%       This phase includes both physical XCT pores and artificial exterior
%       air/padding regions.
%
%   - Solid mask:
%       all phase IDs except pore_phase_id
%
%   - ROI:
%       the analysis region after excluding Z caps and optional interface buffers
%
%   - Pore-surface band:
%       solid grid sites/voxels within a specified shell distance from internal
%       pores
%
%   - Hotspot island:
%       connected component of high averaged von Mises stress satisfying the
%       selected threshold and connectivity criteria
%
%   - Critical-distance averaging:
%       Gaussian smoothing of the stress field normalized by a smoothed solid mask
%       so that pore/air regions do not artificially reduce nearby solid stress
%
% Units:
%   - Coordinates are converted to micrometers, µm.
%   - HDF5 von Mises stress is read in Pa and converted to MPa.
%   - Volumes are reported in grid-site/voxel counts and µm^3.
%   - Defect length metrics are reported in µm.
%
% Notes:
%   - The script assumes the MATLAB preprocessing/generation and Fortran solver
%     used the same padding values pad_xy and pad_z.
%   - The Fortran solver operates on the padded nodal grid. Therefore, /pix and
%     /vm are reshaped using dimensions nx x ny x nz reconstructed here.
%   - Pore/air padding and Z caps can be excluded from the hotspot search using
%     the user settings below.
%   - The sample_type setting controls the severity-score weighting:
%       'EOS' : gas/equiaxed pore dominated
%       'KH'  : keyhole pore dominated
%       'LOF' : lack-of-fusion defect dominated
%
% Required MATLAB functions/toolboxes:
%   - Image Processing Toolbox:
%       imgaussfilt3, imdilate, bwconncomp, bwareaopen, regionprops3
%   - MATLAB built-in HDF5 support:
%       h5read
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


%% ---------------------- USER SETTINGS -----------------------------------
% User-defined analysis parameters.
%
% These values must be consistent with the microstructure-generation script
% and the Fortran elas3d-xtal simulation. In particular, pad_xy, pad_z, and
% dx_um must match the values used when creating input_structure_poly.h5 and
% fullfield_poly.h5.

% --- Padding and grid spacing ---
% pad_xy:
%   Number of artificial air/pore padding voxels added to each side in X and Y.
%
% pad_z:
%   Number of solid cap layers added to each side in Z.
%
% These values are used to reconstruct the padded nodal grid and to exclude
% padding/cap regions from hotspot detection.
pad_xy = 11;               % voxels; must match microstructure generation
pad_z  = 11;               % voxels; must match microstructure generation

% Voxel/grid spacing.
% The Fortran output fields are interpreted on a regular grid with this spacing.
dx_um  = 5.5;              % micrometers per grid step
dx_mm  = dx_um/1000;       % millimeters per grid step

% --- Z-direction ROI exclusion ---
% Exclude the artificial Z caps and, optionally, an additional buffer into the
% gauge region. This helps avoid artificial stress concentrations near the
% cap-gauge interfaces.
%
% z_interface_buf_vox:
%   Extra number of voxels to exclude from the gauge region adjacent to each
%   Z cap. Typical tuning range: 0-25 voxels.
z_interface_buf_vox = 0;      % additional voxels excluded beyond Z caps

% --- Hotspot candidate controls ---
% vm_threshold_MPa:
%   Manual von Mises stress threshold used when use_pct_threshold = false.
%
% N_hotspots:
%   Number of top-ranked hotspot islands to report and plot.
vm_threshold_MPa = 1000;      % manual hotspot threshold, MPa
N_hotspots       = 1;         % number of ranked hotspots to report

% --- Critical-distance averaging and pore-band parameters ---
% mean_d_um:
%   Mean grain diameter used to define the critical averaging length.
%
% crit_dist_frac:
%   Fraction of mean grain diameter used as the characteristic averaging
%   distance. The resulting critical distance is:
%
%       crit_dist_um = crit_dist_frac * mean_d_um
%
% shell_thick_um:
%   Thickness of the solid pore-surface band around internal pores.
%
% pct_thr:
%   Percentile threshold applied to the averaged stress values inside the
%   pore-surface band. Used when use_pct_threshold = true.
mean_d_um        = 52.5;      % mean grain diameter, micrometers
crit_dist_frac   = 0.125;     % fraction of mean grain diameter
shell_thick_um   = 3*dx_um;   % pore-band thickness, micrometers
pct_thr          = 95.0;      % pore-band stress percentile threshold

% ---------------------- SAMPLE TYPE CONFIGURATION -----------------------
% sample_type controls the weighting used in the hotspot severity score.
%
% Options:
%   'EOS' : gas/equiaxed pore dominated morphology
%   'KH'  : keyhole pore dominated morphology
%   'LOF' : lack-of-fusion dominated morphology; includes aspect-ratio weighting
sample_type = 'EOS';

% --- Loading-axis configuration ---
% loading_axis defines the direction normal to the projected defect area used
% in the Murakami-style sqrt(area) metric.
%
% Examples:
%   loading_axis = 'Z' -> projected area is measured on the X-Y plane
%   loading_axis = 'X' -> projected area is measured on the Y-Z plane
%   loading_axis = 'Y' -> projected area is measured on the X-Z plane
loading_axis = 'Z';

% --- Plot controls ---
% font_sz:
%   Base font size for figures.
%
% use_turbo:
%   If true, use MATLAB's turbo colormap. If false, use parula.
%
% str_lim:
%   Upper colorbar limit for von Mises stress plots, in MPa.
font_sz = 15;
use_turbo = true;
str_lim   = 1000;


% --- Fatigue hotspot island controls ---
% use_pore_band_only:
%   If true, hotspot candidates are restricted to the solid pore-surface band.
%   This is recommended when fatigue initiation is expected near pores/defects.
%
% use_pct_threshold:
%   If true, use the pore-band percentile threshold vm_thr.
%   If false, use the manual threshold vm_threshold_MPa.
%
% min_island_vox:
%   Minimum connected-component size retained as a hotspot island. Smaller
%   components are removed as likely numerical specks.
%
% connectivity:
%   3D connectivity used for connected-component labeling.
%   Valid values are typically 6, 18, or 26.
use_pore_band_only   = true;   % restrict candidates to pore band
use_pct_threshold    = true;   % use percentile-based threshold
min_island_vox       = 10;     % minimum hotspot island size, voxels
connectivity         = 26;     % connected-component connectivity

% --- Padding exclusion controls ---
% Exclude artificial padding and cap regions from hotspot detection.
%
% exclude_pad_xy:
%   Removes the artificial XY air/pore padding ring from the search domain.
%
% exclude_pad_z:
%   Removes the artificial Z cap layers from the search domain.
exclude_pad_xy       = true;
exclude_pad_z        = true;

% --- Surface-breaking hotspot boost ---
% If enabled, hotspot islands touching the artificial XY air padding are treated
% as surface-breaking or near-surface candidates and receive a severity boost.
%
% surface_boost_factor:
%   Multiplicative factor applied to the severity score of such islands.
enable_surface_boost = true;
surface_boost_factor = 2.0;

% --- Crystal-plasticity subdomain box sizing ---
% Estimate a local CP simulation box size around each selected hotspot.
%
% If enabled, the box size is based on an equivalent spherical diameter computed
% from the hotspot island volume:
%
%   D_eq = (6*V/pi)^(1/3)
%   box  = max(cp_box_min_um, cp_box_factor*D_eq)
%
% This box is intended to surround the stressed hotspot region, not necessarily
% the physical pore alone.
enable_cp_box_sizing = true;
cp_box_factor        = 3.0;    % multiplier applied to hotspot equivalent diameter
cp_box_min_um        = 90;     % minimum CP box size, micrometers

% --- Visualization window controls ---
% zoom_fullwidth:
%   Physical width of the zoomed 2D hotspot slice plots.
%
% num_ticks:
%   Number of axis ticks used in zoomed plots and 3D visualizations.
zoom_fullwidth = 495.0;   % zoom-window full width, micrometers
num_ticks = 6;

%% ---------------------- OUTPUT LOG SETUP --------------------------------
% Create a text log file that records the main hotspot-analysis settings,
% thresholds, ranked hotspot properties, and CP box-size estimates.
%
% The log is written progressively throughout the script using fid.
log_filename = 'hotspot_analysis_log.txt';
fid = fopen(log_filename, 'w');
if fid == -1
    error('Could not open log file: %s', log_filename);
end


%% ---------------------- LOAD ORIGINAL GRID PARAMETERS -------------------
% Load grid and geometry parameters saved by the microstructure-generation
% script. These describe the original, unpadded XCT cube/domain.
%
% Expected variables:
%   xvec, yvec, zvec : original coordinate vectors, mm
%   Nx, Ny, Nz       : original grid dimensions
%   cube_len         : original cube side length, mm
load('xct_poly_params_SS316L.mat','xvec','yvec','zvec','Nx','Ny','Nz','cube_len');

% Original unpadded grid dimensions.
Nx0 = numel(xvec);  
Ny0 = numel(yvec);  
Nz0 = numel(zvec);

% Original core-domain length in micrometers.
L = cube_len * 1000;

% Reconstruct padded nodal dimensions used by the Fortran solver.
%
% The MATLAB generation script starts from a voxel/cell-centered phase map and
% writes a nodal /pix field for the Fortran code. Therefore:
%
%   nodal dimension = original voxel count + 1 + 2*padding
%
% Since Nx0, Ny0, and Nz0 are grid/node counts for the original cube, the
% number of original voxels is Nx0-1, Ny0-1, and Nz0-1.
nx = (Nx0-1) + 2*pad_xy + 1;
ny = (Ny0-1) + 2*pad_xy + 1;
nz = (Nz0-1) + 2*pad_z  + 1;

% Rebuild padded nodal coordinate vectors.
%
% Coordinates are first constructed in millimeters and then converted to
% micrometers. Index 1 corresponds to the negative padding offset:
%
%   xvec_um(1) = -pad_xy*dx_um
%   yvec_um(1) = -pad_xy*dx_um
%   zvec_um(1) = -pad_z *dx_um
xvec_mm = ((0:nx-1) - pad_xy) * dx_mm;
yvec_mm = ((0:ny-1) - pad_xy) * dx_mm;
zvec_mm = ((0:nz-1) - pad_z ) * dx_mm;

xvec_um = xvec_mm * 1000;
yvec_um = yvec_mm * 1000;
zvec_um = zvec_mm * 1000;

% Original unpadded cube length in micrometers.
L_um = cube_len * 1000;

%% ---------------------- LOAD FIELDS FROM HDF5 ---------------------------
% Load von Mises stress from the Fortran full-field output file.
%
% /vm is assumed to be stored as a 1D array over the padded nodal grid.
% Values are read in Pa and converted to MPa later.
vm = h5read('fullfield_poly.h5', '/vm');           % Pa, length nx*ny*nz
assert(numel(vm) == nx*ny*nz, 'vm length does not match nx*ny*nz.');

% Reshape the 1D stress vector into a 3D array.
%
% Convention used here:
%   vm3d(i,j,k) corresponds to coordinates:
%     xvec_um(i), yvec_um(j), zvec_um(k)
%
% This must match the flattening convention used when writing/reading the
% Fortran HDF5 data.
vm3d = reshape(vm, nx, ny, nz);

% Load the input phase-label field used by the Fortran solver.
%
% /pix contains the padded nodal phase IDs. The maximum phase ID is assumed to
% be the pore/air phase.
pix_input = h5read('input_structure_poly.h5','/pix');
assert(numel(pix_input) == nx*ny*nz, 'input pix length mismatch with nx*ny*nz.');

% Reshape phase labels to the same padded grid convention as vm3d.
phase3d_in = reshape(pix_input, nx, ny, nz);

% Identify the pore/air phase ID from the input phase field.
%
% Convention:
%   solid grain phases = 1:(pore_phase_id-1)
%   pore/air phase     = pore_phase_id
all_phases_in  = unique(phase3d_in(:));
pore_phase_id  = max(all_phases_in);
N_grains_input = pore_phase_id - 1;

% Construct binary masks for pore/air and solid regions.
mask_pore  = (phase3d_in == pore_phase_id);
mask_solid = ~mask_pore;

%% ---------------------- DEFINE ANALYSIS ROI ------------------------------
% Define the Z-direction region of interest used for hotspot analysis.
%
% The padded simulation domain includes artificial solid Z caps. These caps,
% and optionally an additional interface buffer, are excluded to avoid selecting
% hotspots caused by cap/gauge transition artifacts.
%
% z_excl_bottom and z_excl_top are the total numbers of slices excluded from
% the bottom and top of the padded domain, respectively.
z_excl_bottom = pad_z + z_interface_buf_vox;
z_excl_top    = pad_z + z_interface_buf_vox;

% First and last Z indices retained for analysis.
iz_keep_start = 1 + z_excl_bottom;
iz_keep_end   = nz - z_excl_top;

if iz_keep_end <= iz_keep_start
    error('ROI Z-range empty. Reduce z_interface_buf_vox or pad_z.');
end

% mask_roi is true only inside the retained Z-analysis region.
% XY padding is not removed here; it is handled later by mask_core_early.
mask_roi = false(nx, ny, nz);
mask_roi(:,:,iz_keep_start:iz_keep_end) = true;

%% ---------------------- CRITICAL-DISTANCE AVERAGING ----------------------
% Apply solid-aware spatial averaging to the von Mises stress field.
%
% Purpose:
%   Fatigue initiation is often governed by stresses averaged over a finite
%   microstructural length scale rather than a single voxel/node value. Here,
%   the averaging length is defined as a fraction of the mean grain diameter.
%
% Method:
%   1. Set stress to zero outside solid material and outside the Z ROI.
%   2. Smooth the stress field with a 3D Gaussian filter.
%   3. Smooth the solid mask with the same Gaussian filter.
%   4. Divide the smoothed stress by the smoothed solid mask.
%
% This normalized filtering avoids artificially lowering near-pore stresses by
% averaging solid stress with pore/air zeros.
%
% crit_dist_um:
%   Physical critical distance used for averaging.
%
% sigma_crit_vox:
%   Gaussian smoothing width in voxel units. The division by 2.0 converts the
%   chosen critical length into an approximate Gaussian sigma.
crit_dist_um = crit_dist_frac * mean_d_um;
sigma_crit_vox = (crit_dist_um / dx_um) / 2.0;

% Prepare stress field for solid-aware filtering.
% Pore/air regions and regions outside the Z ROI are set to zero before
% Gaussian smoothing.
vm_temp2 = vm3d;
vm_temp2(~mask_solid) = 0;
vm_temp2(~mask_roi)   = 0;

% Smooth stress numerator and solid-mask denominator using the same Gaussian
% filter. Padding with zero prevents artificial wraparound at domain boundaries.
num_c = imgaussfilt3(vm_temp2,            sigma_crit_vox, 'Padding', 0);
den_c = imgaussfilt3(single(mask_solid),  sigma_crit_vox, 'Padding', 0);

% Avoid division by very small denominator values outside/near non-solid regions.
den_c(den_c < 1e-6) = 1;

% Solid-aware averaged stress field.
vm3d_cavg = num_c ./ den_c;

% Force pore/air regions to zero after averaging.
vm3d_cavg(~mask_solid) = 0;

% Apply the final Z ROI mask and convert from Pa to MPa.
vm3d_cavg_roi = vm3d_cavg;
vm3d_cavg_roi(~mask_roi) = 0;
vm_MPa_roi = vm3d_cavg_roi / 1e6;


%% ---------------------- PORE-SURFACE BAND ROI ----------------------------
% Build a solid-region band surrounding internal pores.
%
% Purpose:
%   Fatigue hotspots are expected to form near defects, especially near pore
%   surfaces. This section identifies solid voxels/nodes within a prescribed
%   distance from internal pores, while excluding artificial padding and cap
%   regions.
%
% shell_thick_vox:
%   Pore-band thickness expressed in voxels.
shell_thick_vox = max(1, round(shell_thick_um / dx_um));

% Create a spherical structuring element used to dilate the internal pore mask.
% The dilation radius is shell_thick_vox.
[si,sj,sk] = ndgrid(-shell_thick_vox:shell_thick_vox, ...
    -shell_thick_vox:shell_thick_vox, ...
    -shell_thick_vox:shell_thick_vox);
se_sphere = (si.^2 + sj.^2 + sk.^2) <= shell_thick_vox^2;


% --- Define core analysis mask before pore-band construction ---------------
% mask_core_early excludes artificial padding/cap regions from pore-band
% construction. This prevents the exterior air padding from being interpreted
% as physical XCT pores.
mask_core_early = true(nx, ny, nz);

% Exclude artificial XY air/pore padding if requested.
if exclude_pad_xy
    mask_core_early(1:pad_xy,:,:)         = false;
    mask_core_early(end-pad_xy+1:end,:,:) = false;
    mask_core_early(:,1:pad_xy,:)         = false;
    mask_core_early(:,end-pad_xy+1:end,:) = false;
end

% Exclude artificial Z caps if requested.
if exclude_pad_z
    mask_core_early(:,:,1:pad_z)          = false;
    mask_core_early(:,:,end-pad_z+1:end)  = false;
end

% Also restrict to the Z ROI defined above.
mask_core_early = mask_core_early & mask_roi;

% Identify internal physical pores only.
% This removes artificial air padding and excluded cap/interface regions.
mask_internal_pore = mask_pore & mask_core_early;

% Dilate internal pores and intersect with solid material to obtain a
% solid-only pore-surface band.
dilated_pore       = imdilate(mask_internal_pore, se_sphere);
mask_pore_band     = dilated_pore & mask_solid & mask_roi & mask_core_early;

% Extract averaged von Mises stress values inside the pore-surface band.
% Zero values are removed before percentile thresholding.
band_vm_stress = vm_MPa_roi(mask_pore_band);
band_vm_stress = band_vm_stress(band_vm_stress > 0);

% Percentile-based stress threshold for hotspot detection.
vm_thr = prctile(band_vm_stress, pct_thr);

%% ---------------------- HOTSPOT ISLAND DETECTION -------------------------
% Detect connected high-stress hotspot islands inside the analysis region.
%
% A hotspot island is defined as a connected component of solid material where
% the critical-distance-averaged von Mises stress exceeds the selected threshold.
%
% Candidate region:
%   If use_pore_band_only = true, candidates are restricted to the solid
%   pore-surface band.
%
%   If use_pore_band_only = false, candidates may occur anywhere in the solid
%   analysis ROI.
%
% Threshold:
%   If use_pct_threshold = true, the threshold is the pore-band percentile vm_thr.
%
%   If use_pct_threshold = false, the threshold is the manual value
%   vm_threshold_MPa.
mask_core = mask_core_early;

% Candidate base region.
if use_pore_band_only
    % Restrict hotspot search to solid material near internal pore surfaces.    
    base_mask = mask_pore_band;
else
    % Allow hotspot search throughout the solid Z ROI.    
    base_mask = mask_solid & mask_roi;
end

% Enforce exclusion of artificial padding/cap regions.
base_mask = base_mask & mask_core;

% Select stress threshold.
if use_pct_threshold
    thr_use = vm_thr;               % pore-band percentile threshold, MPa
else
    thr_use = vm_threshold_MPa;     % manual threshold, MPa
end

% High-stress candidate mask.
high_mask = base_mask & (vm_MPa_roi >= thr_use);

% Remove very small connected components that are likely numerical artifacts.
if min_island_vox > 1
    high_mask = bwareaopen(high_mask, min_island_vox, connectivity);
end

% Label connected high-stress regions.
CC = bwconncomp(high_mask, connectivity);
if CC.NumObjects == 0
    error('No hotspot islands found. Lower thr_use, reduce min_island_vox, or increase pore band thickness.');
end

% Compute geometric and stress-weighted properties for each hotspot island.
%
% regionprops3 with intensity image vm_MPa_roi returns:
%   Volume              : island size in voxels
%   MeanIntensity       : mean averaged von Mises stress, MPa
%   MaxIntensity        : peak averaged von Mises stress, MPa
%   WeightedCentroid    : stress-weighted centroid, [x y z] = [column row slice]
%   BoundingBox         : bounding box in image coordinates
%   PrincipalAxisLength : approximate principal-axis lengths of the island
props = regionprops3(CC, vm_MPa_roi, ...
    'Volume', 'MeanIntensity', 'MaxIntensity', 'WeightedCentroid', 'BoundingBox','PrincipalAxisLength');

%% ---------------------- DEFECT / PORE COMPONENTS ------------------------
% Identify physical internal pore components and associate each detected
% high-stress hotspot island with a nearby pore.
%
% Important distinction:
%   - A hotspot island is a connected region of high averaged stress in the
%     surrounding solid material.
%   - A pore/defect component is a connected region of pore voxels/nodes.
%
% These are not the same object. Hotspot island volume should not be interpreted
% as pore volume. This section computes actual pore metrics and links each
% hotspot island to a pore component for defect-aware severity scoring.

fprintf('\n[Defect/Pore Association]\n');

% Use only internal/core pores.
% mask_internal_pore was constructed to exclude artificial XY air padding,
% Z caps, and regions outside the analysis ROI.
mask_defect = mask_internal_pore;

% Label connected internal pore components.
CC_pores = bwconncomp(mask_defect, 26);

if CC_pores.NumObjects == 0
    warning('No internal pore components found. Defect metrics will fall back to hotspot metrics.');
end

% Compute geometric properties of each pore component.
%
% regionprops3 convention:
%   Centroid and WeightedCentroid are reported as:
%     [x y z] = [column row slice]
%
% PrincipalAxisLength gives approximate lengths of the principal axes of the
% connected component in voxel units.
if CC_pores.NumObjects > 0
    pore_props = regionprops3(CC_pores, ...
        'Volume', ...
        'Centroid', ...
        'BoundingBox', ...
        'PrincipalAxisLength');
else
    pore_props = table();
end

% Convert connected-component labels into a label matrix for fast lookup.
% pore_label(m) = pore component ID at that grid site, or 0 if not pore.
pore_label = labelmatrix(CC_pores);

% Number of internal pore components.
n_pores = CC_pores.NumObjects;

% Allocate pore metric arrays.
%
% pore_V_vox:
%   Pore component volume in grid sites/voxels.
%
% pore_V_um3:
%   Pore component volume in micrometers cubed.
%
% pore_Deq_um:
%   Equivalent spherical diameter:
%     Deq = (6*V/pi)^(1/3)
%
% pore_sqrtArea_um:
%   Murakami-style length scale sqrt(projected area), where projection is
%   taken normal to loading_axis.
%
% pore_aspect_ratio:
%   Ratio of largest to smallest principal axis length.
pore_V_vox          = zeros(n_pores,1);
pore_V_um3          = zeros(n_pores,1);
pore_Deq_um         = zeros(n_pores,1);
pore_sqrtArea_um    = zeros(n_pores,1);
pore_aspect_ratio   = ones(n_pores,1);



for p = 1:n_pores
    % Linear indices belonging to pore component p.    
    vox_p = CC_pores.PixelIdxList{p};

    % Volume metrics.    
    pore_V_vox(p) = numel(vox_p);
    pore_V_um3(p) = pore_V_vox(p) * dx_um^3;

    % Equivalent spherical diameter based on pore volume.
    pore_Deq_um(p) = (6 * pore_V_um3(p) / pi)^(1/3);

    % Convert linear indices to array subscripts.
    % Array convention:
    %   phase3d_in(i,j,k) corresponds to X,Y,Z.
    [ip, jp, kp] = ind2sub([nx, ny, nz], vox_p);

    % Compute projected pore area normal to the loading direction.
    % The projected area is estimated by counting unique occupied grid cells
    % in the projection plane and multiplying by dx_um^2.    
    switch upper(loading_axis)
        case 'X'
            % Projection normal to X: project onto Y-Z plane.
            proj_pairs = unique([jp, kp], 'rows');

        case 'Y'
            % Projection normal to Y: project onto X-Z plane.
            proj_pairs = unique([ip, kp], 'rows');

        case 'Z'
            % Projection normal to Z: project onto X-Y plane.
            proj_pairs = unique([ip, jp], 'rows');

        otherwise
            error('Unknown loading_axis: %s. Use X, Y, or Z.', loading_axis);
    end

    area_proj_um2 = size(proj_pairs,1) * dx_um^2;

    % Murakami-style defect length scale.
    pore_sqrtArea_um(p) = sqrt(area_proj_um2);

    % Pore aspect ratio from principal-axis lengths.
    % If the shape measurement is invalid or degenerate, use a neutral factor 1.
    PAL_pore = pore_props.PrincipalAxisLength(p,:);

    if all(isfinite(PAL_pore)) && min(PAL_pore) > 0
        pore_aspect_ratio(p) = max(PAL_pore) / max(min(PAL_pore), eps);
    else
        pore_aspect_ratio(p) = 1.0;
    end
end

fprintf('  Pore components found: %d\n', n_pores);

% Associate each hotspot island with one pore component.
%
% Primary method:
%   1. Dilate the hotspot island by assoc_radius_vox.
%   2. Find pore labels intersecting the dilated island.
%   3. Assign the hotspot to the pore with the largest intersection/contact
%      count.
%
% Fallback method:
%   If no pore intersects the dilated hotspot but pore components exist, assign
%   the hotspot to the nearest pore centroid.
%
% If no pore components exist:
%   The hotspot is left unassociated, and downstream severity metrics fall back
%   to hotspot-based length/shape metrics.
assoc_radius_vox = max(1, shell_thick_vox);

% Spherical structuring element for hotspot-pore association.
[ai, aj, ak] = ndgrid(-assoc_radius_vox:assoc_radius_vox, ...
    -assoc_radius_vox:assoc_radius_vox, ...
    -assoc_radius_vox:assoc_radius_vox);

se_assoc = (ai.^2 + aj.^2 + ak.^2) <= assoc_radius_vox^2;

% Association results for each hotspot island.
assoc_pore_id      = zeros(CC.NumObjects,1);   % associated pore ID, or 0
assoc_method       = strings(CC.NumObjects,1); % association method used
assoc_contact_vox  = zeros(CC.NumObjects,1);   % contact/intersection count

% Convert pore centroids from regionprops3 image coordinates to array-index
% coordinates used by this script.
%
% regionprops3 returns:
%   Centroid = [x y z] = [column row slice]
%
% Script array indices are:
%   [i j k] = [row-like x-index, column-like y-index, z-index]
%
% Therefore:
%   i <- y
%   j <- x
%   k <- z
if n_pores > 0
    pore_centroid_raw = pore_props.Centroid;
    pore_centroid_ijk = [pore_centroid_raw(:,2), ...
        pore_centroid_raw(:,1), ...
        pore_centroid_raw(:,3)];
else
    pore_centroid_ijk = [];
end

for r = 1:CC.NumObjects
    % Binary mask for hotspot island r.    
    island = false(nx, ny, nz);
    island(CC.PixelIdxList{r}) = true;

    % Dilate hotspot to find nearby/touching pore components.
    island_dilated = imdilate(island, se_assoc);

    % Pore component labels intersecting the dilated hotspot.    
    labels_touching = pore_label(island_dilated);
    labels_touching = labels_touching(labels_touching > 0);

    if ~isempty(labels_touching)
        % Assign to the pore with the largest intersection/contact count.
        counts = accumarray(double(labels_touching(:)), 1, [n_pores, 1]);
        [assoc_contact_vox(r), assoc_pore_id(r)] = max(counts);
        assoc_method(r) = "dilated-contact";

    elseif n_pores > 0
        % Fallback: assign to nearest pore centroid.
        wc = props.WeightedCentroid(r,:);

        % WeightedCentroid is [x y z] = [column row slice].
        % Convert to script index convention [i j k].
        hotspot_centroid_ijk = [wc(2), wc(1), wc(3)];

        d2 = sum((pore_centroid_ijk - hotspot_centroid_ijk).^2, 2);
        [~, assoc_pore_id(r)] = min(d2);
        assoc_contact_vox(r) = 0;
        assoc_method(r) = "nearest-centroid";

    else
        % No pore components are available for association.
        assoc_pore_id(r) = 0;
        assoc_contact_vox(r) = 0;
        assoc_method(r) = "none";
    end
end


%% ---------------------- DEFECT-AWARE SEVERITY SCORE ----------------------
% Compute a ranking score for each hotspot island.
%
% The severity score combines:
%   S1 : distributed hotspot energy-like term
%   S2 : defect-size and peak-stress term
%   S3 : peak-weighted hotspot energy-like term
%   S4 : defect/hotspot shape term, used for LOF defects
%
% Key distinction:
%   - Hotspot volume describes the stressed solid region.
%   - Defect metrics describe the associated pore/defect geometry.
%
% Therefore:
%   S1 and S3 use hotspot island volume.
%   S2 uses the associated pore projected sqrt(area), when available.
%   S4 uses associated pore aspect ratio, when available.
%
% Each score component is min-max normalized before weighting.

sample_type = upper(string(sample_type));

% --- Extract hotspot island properties -----------------------------------
% regionprops3 returns Volume in voxel/grid-site counts and stress intensities
% from vm_MPa_roi in MPa.
V_vox      = props.Volume;                % hotspot island volume, voxels
V_um3      = V_vox * dx_um^3;             % hotspot island volume, µm^3
sigma_mean = props.MeanIntensity;         % mean hotspot stress, MPa
sigma_max  = props.MaxIntensity;          % peak hotspot stress, MPa

% Characteristic hotspot length used as fallback if no pore is associated.
L_hot_um = V_um3.^(1/3);

% Estimate hotspot aspect ratio from principal-axis lengths.
% This is used as a fallback shape metric when no associated pore is available.
PAL_hot = props.PrincipalAxisLength;

hotspot_aspect_ratio = ones(CC.NumObjects,1);

for r = 1:CC.NumObjects
    pal_r = PAL_hot(r,:);

    if all(isfinite(pal_r)) && min(pal_r) > 0
        hotspot_aspect_ratio(r) = max(pal_r) / max(min(pal_r), eps);
    else
        hotspot_aspect_ratio(r) = 1.0;
    end
end

% --- Defect-aware metrics associated with each hotspot --------------------
% Initialize with hotspot-based fallback values.
% These are overwritten when a valid associated pore is available.
defect_Deq_um       = L_hot_um;
defect_sqrtArea_um  = L_hot_um;
defect_aspect_ratio = hotspot_aspect_ratio;

for r = 1:CC.NumObjects
    pid = assoc_pore_id(r);

    if pid > 0
        defect_Deq_um(r)       = pore_Deq_um(pid);
        defect_sqrtArea_um(r)  = pore_sqrtArea_um(pid);
        defect_aspect_ratio(r) = pore_aspect_ratio(pid);
    end
end

% --- Individual severity components --------------------------------------

% S1: Distributed hotspot energy-like proxy.
% Uses hotspot island volume because this term represents the volume of highly
% stressed solid material.
S1 = V_um3 .* sigma_mean.^2;

% S2: Defect-size plus peak-stress proxy.
% Uses associated pore sqrt(projected area), following Murakami-style defect
% size reasoning. If no pore is associated, a hotspot length fallback is used.
S2 = sigma_max .* defect_sqrtArea_um;

% S3: Conservative peak-weighted hotspot energy proxy.
% Uses hotspot island volume and peak stress.
S3 = V_um3 .* sigma_max.^2;

% S4: Shape/aspect-ratio factor.
% Mainly relevant for LOF defects, where elongated or flat defects are more
% crack-like than rounded pores.
S4 = defect_aspect_ratio;

% --- Safe min-max normalization ------------------------------------------
% Each score component is normalized to [0,1]. If all candidates have the same
% value for a component, minmax_safe returns ones so that the component remains
% neutral rather than zeroing all severity scores.
S1_norm = minmax_safe(S1);
S2_norm = minmax_safe(S2);
S3_norm = minmax_safe(S3);
S4_norm = minmax_safe(S4);

% --- Weight selection by defect/sample type -------------------------------
% The severity-score weights are heuristic and can be tuned for a specific
% material, defect class, or calibration dataset.
%
% EOS:
%   Rounded/equiaxed gas-pore morphology. Distributed hotspot energy receives
%   slightly higher weight.
%
% KH:
%   Keyhole-pore morphology. Balanced weighting between distributed stress and
%   peak-weighted hotspot severity.
%
% LOF:
%   Lack-of-fusion morphology. Includes an additional shape/aspect-ratio term
%   because LOF defects are often elongated and crack-like.
switch sample_type

    case "EOS"
        % Rounded/equiaxed gas pores.
        % Distributed energy is slightly emphasized.
        w1 = 0.40;
        w2 = 0.30;
        w3 = 0.30;

        severity = w1*S1_norm + w2*S2_norm + w3*S3_norm;

        weight_description = sprintf( ...
            'EOS weights: S1=%.2f, S2=%.2f, S3=%.2f', ...
            w1, w2, w3);

    case "KH"
        % Keyhole pores.
        % Balanced between distributed energy and peak-weighted energy.
        w1 = 0.35;
        w2 = 0.30;
        w3 = 0.35;

        severity = w1*S1_norm + w2*S2_norm + w3*S3_norm;

        weight_description = sprintf( ...
            'KH weights: S1=%.2f, S2=%.2f, S3=%.2f', ...
            w1, w2, w3);

    case "LOF"
        % Lack-of-fusion defects.
        % More crack-like, so include shape/aspect ratio.
        %
        % S1: hotspot energy
        % S2: defect projected sqrt(area) + peak stress
        % S3: peak-weighted hotspot energy
        % S4: defect aspect ratio / crack-like morphology
        w1 = 0.20;
        w2 = 0.30;
        w3 = 0.35;
        w4 = 0.15;

        assert(abs(w1 + w2 + w3 + w4 - 1.0) < 1e-9, ...
            'LOF weights must sum to 1.0');

        severity = w1*S1_norm + w2*S2_norm + w3*S3_norm + w4*S4_norm;

        weight_description = sprintf( ...
            'LOF weights: S1=%.2f, S2=%.2f, S3=%.2f, S4(shape)=%.2f', ...
            w1, w2, w3, w4);

    otherwise
        error('Unknown sample_type: %s. Use EOS, LOF, or KH.', sample_type);
end


%% ---- Surface-breaking boost --------------------------------------------
% Optionally boost the severity score for hotspot islands that touch the
% artificial XY air padding.
%
% Interpretation:
%   If a hotspot island touches the external air/pore padding, it may represent
%   a surface-breaking or near-surface defect. Such defects can be more critical
%   for fatigue and are optionally ranked more aggressively.
%
% Note:
%   This boost affects ranking only. It does not modify stress values.

if enable_surface_boost && exclude_pad_xy

    % Build a mask for the artificial XY air padding region.
    % Z caps are not included here; this boost is intended for lateral
    % surface-breaking features.    
    mask_air = false(nx, ny, nz);

    % XY padding air.
    mask_air(1:pad_xy,:,:)         = true;
    mask_air(end-pad_xy+1:end,:,:) = true;
    mask_air(:,1:pad_xy,:)         = true;
    mask_air(:,end-pad_xy+1:end,:) = true;

    % Track whether each hotspot island touches the artificial air region.    
    touches_air = false(CC.NumObjects, 1);

    % 6-neighbor plus center structuring element used to test direct surface
    % contact without allowing diagonal-only contact to dominate.
    se6 = false(3,3,3);
    se6(2,2,1) = true;  % -Z
    se6(2,2,3) = true;  % +Z
    se6(2,1,2) = true;  % -Y
    se6(2,3,2) = true;  % +Y
    se6(1,2,2) = true;  % -X
    se6(3,2,2) = true;  % +X
    se6(2,2,2) = true;  % center

    for r = 1:CC.NumObjects
        % Binary mask for the current hotspot island.        
        island = false(nx, ny, nz);
        island(CC.PixelIdxList{r}) = true;

        % Dilate by the 6-neighbor stencil and test for contact with air padding.
        island_nbhd = imdilate(island, se6);

        touches_air(r) = any(mask_air(island_nbhd), 'all');
    end

    % Apply multiplicative severity boost to surface-touching islands.
    severity_boosted = severity;
    severity_boosted(touches_air) = severity_boosted(touches_air) * surface_boost_factor;

else
    % No surface boost requested.
    touches_air = false(CC.NumObjects, 1);
    severity_boosted = severity;
end


% Sort hotspot islands from most to least severe.
% If surface boosting is enabled, severity_boosted is used for ranking.
[sorted_score, sort_idx] = sort(severity_boosted, 'descend');

% Limit output to the requested number of hotspots.
num_hotspots = min(N_hotspots, numel(sorted_score));

% Extract IDs, properties, and scores for the top-ranked hotspots.
top_idx   = sort_idx(1:num_hotspots);
top_props = props(top_idx, :);
top_score = sorted_score(1:num_hotspots);


%% ---- Get weighted-centroid voxel indices --------------------------------
% Convert hotspot weighted centroids from regionprops3 coordinates to array
% indices used by this script.
%
% regionprops3 WeightedCentroid convention:
%   [x y z] = [column row slice]
%
% Script array convention:
%   vm_MPa_roi(i,j,k) = value at [X_1, X_2, X_3]
%
% Therefore:
%   i = row-like coordinate = WeightedCentroid(:,2)
%   j = column-like coordinate = WeightedCentroid(:,1)
%   k = slice coordinate = WeightedCentroid(:,3)
%
% The resulting indices are rounded to the nearest grid site and clamped to
% valid array bounds.

ci_raw = round(top_props.WeightedCentroid(:,2));
cj_raw = round(top_props.WeightedCentroid(:,1));
ck_raw = round(top_props.WeightedCentroid(:,3));

% Clamp rounded centroids to valid padded-grid bounds.
ci_raw = max(1, min(nx, ci_raw));
cj_raw = max(1, min(ny, cj_raw));
ck_raw = max(1, min(nz, ck_raw));


% The stress-weighted centroid can occasionally round to a pore/air location
% or outside the ROI, especially for irregular islands near defects.
%
% If the rounded centroid is invalid, replace it with the peak-stress voxel
% inside the same connected hotspot island. This ensures the selected point is
% physically located in solid material inside the analysis ROI.
ci = zeros(num_hotspots, 1);
cj = zeros(num_hotspots, 1);
ck = zeros(num_hotspots, 1);

was_snapped = false(num_hotspots, 1);

for t = 1:num_hotspots

    % Connected-component ID of this ranked hotspot.
    island_id = top_idx(t);
  
    % Linear indices of all voxels/nodes in this hotspot island.
    island_linear_idx = CC.PixelIdxList{island_id};

    % Convert island voxels to subscripts for fallback peak selection.
    [ii_all, jj_all, kk_all] = ind2sub([nx, ny, nz], island_linear_idx);

    % Candidate centroid index.
    raw_lin = sub2ind([nx, ny, nz], ci_raw(t), cj_raw(t), ck_raw(t));

    % Accept centroid only if it lies in solid material and inside the ROI.
    centroid_is_valid = mask_solid(raw_lin) && mask_roi(raw_lin);

    if centroid_is_valid
        % Use the rounded weighted centroid.
        ci(t) = ci_raw(t);
        cj(t) = cj_raw(t);
        ck(t) = ck_raw(t);
    else
        % Fallback: use the maximum-stress voxel within this hotspot island.
        island_stress = vm_MPa_roi(island_linear_idx);
        [~, peak_local_idx] = max(island_stress);

        ci(t) = ii_all(peak_local_idx);
        cj(t) = jj_all(peak_local_idx);
        ck(t) = kk_all(peak_local_idx);

        was_snapped(t) = true;

        fprintf(['  [WARN] Hotspot #%d centroid was invalid/pore/outside ROI — ' ...
            'snapped to peak-stress voxel at (%d,%d,%d)\n'], ...
            t, ci(t), cj(t), ck(t));
    end
end

ci = max(1, min(nx, ci));
cj = max(1, min(ny, cj));
ck = max(1, min(nz, ck));


%% ---- Coordinates in padded simulation frame -----------------------------
% Convert selected hotspot indices to physical coordinates in the padded
% simulation frame.
%
% In this frame:
%   xvec_um(1) = -pad_xy*dx_um
%   yvec_um(1) = -pad_xy*dx_um
%   zvec_um(1) = -pad_z *dx_um
%
% These coordinates are consistent with the stress and phase arrays used for
% plotting in this script.
h_x = xvec_um(ci);
h_y = yvec_um(cj);
h_z = zvec_um(ck);

% Backward-compatible variable names used by later plotting/logging sections.
top_i    = ci;
top_j    = cj;
top_k    = ck;
top_vals = top_props.MaxIntensity;


%% ---- Coordinates in unpadded cube frame ---------------------------------
% Convert hotspot coordinates from the padded simulation frame to the original
% unpadded cube frame.
%
% Padded-frame origin:
%   The first index lies at negative padding distance.
%
% Unpadded cube-frame origin:
%   The original XCT cube starts at 0.
%
% These coordinates may be useful when mapping hotspot locations back to the
% original XCT/generation domain.

h_x_cube = h_x + pad_xy * dx_um;
h_y_cube = h_y + pad_xy * dx_um;
h_z_cube = h_z + pad_z  * dx_um;


% Print a short console summary. A more detailed summary is written later to
% the log file.
fprintf('\nFatigue hotspot islands using thr=%.2f MPa\n', thr_use);
fprintf('Found %d islands, reporting top %d\n', CC.NumObjects, num_hotspots);


%% ---- CP subdomain box sizing --------------------------------------------
% Estimate crystal-plasticity subdomain box size around each selected hotspot.
%
% The box size is based on the hotspot island volume, not the physical pore
% volume. This is intentional: the CP subdomain should cover the stressed solid
% region surrounding the defect.
%
% Equivalent hotspot diameter:
%   D_eq = (6*V_hotspot/pi)^(1/3)
%
% CP box rule:
%   box_um = max(cp_box_min_um, cp_box_factor * D_eq)
%
% box_vox converts the physical box size to a voxel/grid-step count.

if enable_cp_box_sizing

    % Hotspot island volume in µm^3.
    V_top_um3 = top_props.Volume * dx_um^3;

    % Equivalent spherical diameter of the selected hotspot island.
    D_eq_um = (6 * V_top_um3 / pi).^(1/3);

    % CP box size in physical units and grid steps.
    box_um  = max(cp_box_min_um, cp_box_factor * D_eq_um);
    box_vox = ceil(box_um / dx_um);

else
    % CP box sizing disabled.
    D_eq_um = nan(num_hotspots, 1);
    box_um  = nan(num_hotspots, 1);
    box_vox = nan(num_hotspots, 1);
end


%% ---------------------- LOG: CRITICAL-DISTANCE AVERAGING ----------------
% Write the critical-distance averaging settings to both the console and the
% text log file.
%
% These values document the length scale used to smooth/average the von Mises
% stress field before hotspot detection.
str = sprintf('[Critical Distance Averaging]\n');
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  mean_d_um      = %.2f um\n', mean_d_um);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  crit_dist_frac = %.2f\n', crit_dist_frac);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  crit_dist_um   = %.2f um\n', crit_dist_um);
fprintf('%s', str); fprintf(fid, '%s', str);

%% ---------------------- LOG: PORE-SURFACE BAND --------------------------
% Write pore-surface band parameters and the percentile-based stress threshold.
%
% vm_thr is computed from the critical-distance-averaged von Mises stresses
% inside mask_pore_band.
str = sprintf('\n[Pore Surface Band]\n');
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  shell_thick_um  = %.2f um\n', shell_thick_um);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  shell_thick_vox = %d vox\n', shell_thick_vox);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  Percentile threshold on pore-band: p%.1f = %.2f MPa\n', pct_thr, vm_thr);
fprintf('%s', str); fprintf(fid, '%s', str);

%% ---------------------- LOG: FATIGUE HOTSPOT ISLANDS --------------------
% Write hotspot-detection settings and a ranked summary of selected hotspot
% islands.
%
% Notes:
%   - Scores are severity scores after optional surface-breaking boost.
%   - Mean/Max stress values are from the critical-distance-averaged von Mises
%     field, vm_MPa_roi.
%   - Coordinates are reported in the padded simulation frame.
str = sprintf('\n[Fatigue Hotspot Islands]\n');
fprintf('%s', str); fprintf(fid, '%s', str);

% Document whether the selected threshold came from pore-band percentile
% thresholding or from the manual threshold setting.
if use_pct_threshold
    thr_str = 'pore-band percentile';
else
    thr_str = 'manual';
end

str = sprintf('  Threshold used      = %.2f MPa (%s)\n', thr_use, thr_str);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  use_pore_band_only  = %d\n', use_pore_band_only);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  min_island_vox      = %d\n', min_island_vox);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  connectivity        = %d\n', connectivity);
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('  Total islands found = %d | Reporting top %d\n', CC.NumObjects, num_hotspots);
fprintf('%s', str); fprintf(fid, '%s', str);

% Table header for selected hotspot islands.
str = sprintf('\n  %-4s %-12s %-8s %-10s %-10s %-10s %-10s %-10s %-16s\n', ...
    '#', 'Score', 'Vol(vox)', 'Mean(MPa)', 'Max(MPa)', ...
    'x(um)', 'y(um)', 'z(um)', 'idx(i,j,k)');
fprintf('%s', str); fprintf(fid, '%s', str);

str = sprintf('  %s\n', repmat('-', 1, 100));
fprintf('%s', str); fprintf(fid, '%s', str);

% Write one row per reported hotspot.
for t = 1:num_hotspots
    ii = top_i(t); jj = top_j(t); kk = top_k(t);

    % Whether the global hotspot island touches the artificial XY air padding.
    % This string is currently computed but not printed in the table below.
    % Add it to the table if surface-breaking status should be logged.    
    if touches_air(top_idx(t))
        surf_str = 'YES*';
    else
        surf_str = 'no';
    end

    global_island_id = top_idx(t);

    str = sprintf('  #%-3d %-12.3e %-8d %-10.1f %-10.1f %-10.1f %-10.1f %-10.1f (%d,%d,%d)\n', ...
        t, top_score(t), top_props.Volume(t), ...
        top_props.MeanIntensity(t), top_props.MaxIntensity(t), ...
        h_x(t), h_y(t), h_z(t), ii, jj, kk);  

    fprintf('%s', str); 
    fprintf(fid, '%s', str);
end

%% ---------------------- LOG: CP BOX SIZING ------------------------------
% Write CP subdomain box-size estimates for selected hotspots.
%
% The CP box size is based on the equivalent diameter of the hotspot island
% volume, not directly on the associated pore volume.
if enable_cp_box_sizing
    str = sprintf('\n[CP Subdomain Box Sizing]\n');
    fprintf('%s', str); fprintf(fid, '%s', str);
    str = sprintf('  Rule: box = max(%.1f um, %.1f x D_eq)\n', cp_box_min_um, cp_box_factor);
    fprintf('%s', str); fprintf(fid, '%s', str);

    str = sprintf('  %-4s %-12s %-12s %-12s\n', '#', 'D_eq(um)', 'Box(um)', 'Box(vox)');
    fprintf('%s', str); fprintf(fid, '%s', str);
    str = sprintf('  %s\n', repmat('-', 1, 44));
    fprintf('%s', str); fprintf(fid, '%s', str);

    for t = 1:num_hotspots
        str = sprintf('  #%-3d %-12.1f %-12.1f %-12d\n', ...
            t, D_eq_um(t), box_um(t), box_vox(t));
        fprintf('%s', str); fprintf(fid, '%s', str);
    end
end

%% ---------------------- CLOSE LOG AND SAVE WORKSPACE --------------------
% Finalize the analysis log and close the file handle.
str = sprintf('\n========================================\n');
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf(' Analysis complete: %s\n', datetime("now"));
fprintf('%s', str); fprintf(fid, '%s', str);
str = sprintf('========================================\n');
fprintf('%s', str); fprintf(fid, '%s', str);

fclose(fid);
fprintf('Log saved to: %s\n', log_filename);

% Save all current workspace variables to a .mat file
save('xct_poly_post_SS316L_full.mat','-v7.3')



%--------------------------------------------------------------------------
% HOTSPOT YZ-SLICE PLOTS
zoom_halfwidth_y_um   = zoom_fullwidth*0.5;
zoom_halfwidth_z_um   = zoom_fullwidth*0.5;
zoom_ticks            = num_ticks;

z_roi_idx = iz_keep_start:iz_keep_end;
z_roi_vec = zvec_um(z_roi_idx);

y_min_dom = min(yvec_um);  y_max_dom = max(yvec_um);
z_min_dom = min(z_roi_vec); z_max_dom = max(z_roi_vec);

for k = 1:num_hotspots

    cur_i   = top_i(k);
    cur_j   = top_j(k);
    cur_k   = top_k(k);
    cur_val = top_vals(k);

    x_coord = xvec_um(cur_i);
    peak_y  = yvec_um(cur_j);
    peak_z  = zvec_um(cur_k);

    vm_slice = squeeze(vm_MPa_roi(cur_i, :, :))';
    vm_slice_roi = vm_slice(z_roi_idx, :);

    [YY, ZZ] = meshgrid(yvec_um, z_roi_vec);

    Wy = 2*zoom_halfwidth_y_um;
    Wz = 2*zoom_halfwidth_z_um;

    if Wy >= (y_max_dom - y_min_dom)
        y_lo = y_min_dom; y_hi = y_max_dom;
    else
        y_lo = peak_y - zoom_halfwidth_y_um;
        y_hi = peak_y + zoom_halfwidth_y_um;

        if y_lo < y_min_dom
            shift = y_min_dom - y_lo;
            y_lo = y_lo + shift;
            y_hi = y_hi + shift;
        end
        if y_hi > y_max_dom
            shift = y_hi - y_max_dom;
            y_lo = y_lo - shift;
            y_hi = y_hi - shift;
        end

        y_lo = max(y_lo, y_min_dom);
        y_hi = min(y_hi, y_max_dom);
    end

    if Wz >= (z_max_dom - z_min_dom)
        z_lo = z_min_dom; z_hi = z_max_dom;
    else
        z_lo = peak_z - zoom_halfwidth_z_um;
        z_hi = peak_z + zoom_halfwidth_z_um;

        if z_lo < z_min_dom
            shift = z_min_dom - z_lo;
            z_lo = z_lo + shift;
            z_hi = z_hi + shift;
        end
        if z_hi > z_max_dom
            shift = z_hi - z_max_dom;
            z_lo = z_lo - shift;
            z_hi = z_hi - shift;
        end

        z_lo = max(z_lo, z_min_dom);
        z_hi = min(z_hi, z_max_dom);
    end


    figure; clf; hold on;   
    contourf(YY, ZZ, vm_slice_roi, 600, 'LineColor','none');
    if use_turbo, colormap(turbo); else, colormap(parula); end

    cb = colorbar;
    cb.Label.String   = 'von Mises Stress (MPa)';
    cb.Label.FontSize = font_sz;
    cb.Location       = 'eastoutside';
    clim([0, str_lim]);

    plot(peak_y, peak_z, 'kp', 'MarkerSize', 20, 'MarkerFaceColor', ...
        'magenta', 'LineWidth', 1.8);
    text(peak_y, peak_z, sprintf('  %.1f MPa', cur_val), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.80], 'EdgeColor', 'k', ...
        'VerticalAlignment', 'bottom');

    xlabel('X_2 (µm)');
    ylabel('X_3 (µm)');
    title(sprintf('X_1=%.1f µm | (X_2,X_3)=(%.1f, %.1f) µm', ...
        x_coord, peak_y, peak_z)); 
    set(gca, 'FontSize', font_sz);

    axis equal tight;
    xlim([y_lo y_hi]);
    ylim([z_lo z_hi]);
    xticks(linspace(y_lo, y_hi, zoom_ticks));
    yticks(linspace(z_lo, z_hi, zoom_ticks));
    grid on; box on;

    hold off;
end


%--------------------------------------------------------------------------
% HOTSPOT XY-SLICE PLOTS
zoom_halfwidth_x_um    = zoom_fullwidth*0.5;
zoom_halfwidth_y_um_xy = zoom_fullwidth*0.5;
zoom_ticks_xy          = num_ticks;


x_min_dom = min(xvec_um);
x_max_dom = max(xvec_um);

y_min_dom = min(yvec_um);
y_max_dom = max(yvec_um);

for k = 1:num_hotspots

    cur_i   = top_i(k);
    cur_j   = top_j(k);
    cur_k   = top_k(k);
    cur_val = top_vals(k);

    peak_x  = xvec_um(cur_i);
    peak_y  = yvec_um(cur_j);
    z_coord = zvec_um(cur_k);

    vm_slice_xy = squeeze(vm_MPa_roi(:, :, cur_k))';
    [XX, YY] = meshgrid(xvec_um, yvec_um);

    Wx = 2 * zoom_halfwidth_x_um;
    Wy = 2 * zoom_halfwidth_y_um_xy;

    if Wx >= (x_max_dom - x_min_dom)
        x_lo = x_min_dom;
        x_hi = x_max_dom;
    else
        x_lo = peak_x - zoom_halfwidth_x_um;
        x_hi = peak_x + zoom_halfwidth_x_um;

        if x_lo < x_min_dom
            shift = x_min_dom - x_lo;
            x_lo = x_lo + shift;
            x_hi = x_hi + shift;
        end

        if x_hi > x_max_dom
            shift = x_hi - x_max_dom;
            x_lo = x_lo - shift;
            x_hi = x_hi - shift;
        end

        x_lo = max(x_lo, x_min_dom);
        x_hi = min(x_hi, x_max_dom);
    end

    if Wy >= (y_max_dom - y_min_dom)
        y_lo = y_min_dom;
        y_hi = y_max_dom;
    else
        y_lo = peak_y - zoom_halfwidth_y_um_xy;
        y_hi = peak_y + zoom_halfwidth_y_um_xy;

        if y_lo < y_min_dom
            shift = y_min_dom - y_lo;
            y_lo = y_lo + shift;
            y_hi = y_hi + shift;
        end

        if y_hi > y_max_dom
            shift = y_hi - y_max_dom;
            y_lo = y_lo - shift;
            y_hi = y_hi - shift;
        end

        y_lo = max(y_lo, y_min_dom);
        y_hi = min(y_hi, y_max_dom);
    end

    figure; clf; hold on;

    contourf(XX, YY, vm_slice_xy, 600, 'LineColor', 'none');

    if use_turbo, colormap(turbo); else, colormap(parula); end

    cb = colorbar;
    cb.Label.String   = 'von Mises Stress (MPa)';
    cb.Label.FontSize = font_sz;
    cb.Location       = 'eastoutside';
    clim([0, str_lim]);

    plot(peak_x, peak_y, 'kp', ...
        'MarkerSize', 20, ...
        'MarkerFaceColor', 'magenta', ...
        'LineWidth', 1.8);

    text(peak_x, peak_y, sprintf('  %.1f MPa', cur_val), ...
        'FontSize', 13, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.80], ...
        'EdgeColor', 'k', ...
        'VerticalAlignment', 'bottom');

    xlabel('X_1 (µm)');
    ylabel('X_2 (µm)');

    title(sprintf('X_3=%.1f µm | (X_1,X_2)=(%.1f, %.1f) µm', ...
        z_coord, peak_x, peak_y)); 

    set(gca, 'FontSize', font_sz);

    axis equal tight;
    xlim([x_lo, x_hi]);
    ylim([y_lo, y_hi]);

    xticks(linspace(x_lo, x_hi, zoom_ticks_xy));
    yticks(linspace(y_lo, y_hi, zoom_ticks_xy));

    grid on; box on;

    hold off;
end






% -------------------------------------------------------------------------
% 3D plot
figure; clf; hold on;

mask_solid_roi = mask_solid & mask_roi;
exposed = false(nx, ny, nz);
exposed(1:end-1,:,:) = exposed(1:end-1,:,:) | ~mask_solid_roi(2:end,:,:);
exposed(2:end,:,:)   = exposed(2:end,:,:)   | ~mask_solid_roi(1:end-1,:,:);
exposed(:,1:end-1,:) = exposed(:,1:end-1,:) | ~mask_solid_roi(:,2:end,:);
exposed(:,2:end,:)   = exposed(:,2:end,:)   | ~mask_solid_roi(:,1:end-1,:);
exposed(:,:,1:end-1) = exposed(:,:,1:end-1) | ~mask_solid_roi(:,:,2:end);
exposed(:,:,2:end)   = exposed(:,:,2:end)   | ~mask_solid_roi(:,:,1:end-1);

z_restrict = false(nx, ny, nz);
z_restrict(:,:,iz_keep_start:iz_keep_end) = true;
surface_mask_roi = mask_solid_roi & exposed & z_restrict;

surf_idx     = find(surface_mask_roi);
[ii, jj, kk] = ind2sub([nx, ny, nz], surf_idx);

Xs = xvec_um(ii);
Ys = yvec_um(jj);
Zs = zvec_um(kk);
Cs = vm_MPa_roi(surf_idx);

x_min = min(xvec)*1e3; x_max = max(xvec)*1e3;
y_min = min(yvec)*1e3; y_max = max(yvec)*1e3;
z_min = min(zvec)*1e3; z_max = max(zvec)*1e3;

x_tick = linspace(x_min, x_max, num_ticks);
y_tick = linspace(y_min, y_max, num_ticks);
z_tick = linspace(z_min, z_max, num_ticks);

hSurf = scatter3(Xs, Ys, Zs, 12, Cs, 'filled');
hSurf.MarkerFaceAlpha = 0.55;
hSurf.MarkerEdgeAlpha = 0.20;

if use_turbo, colormap(turbo); else, colormap(parula); end
cb = colorbar;
cb.Label.String   = 'von Mises Stress (MPa)';
cb.Label.FontSize = font_sz;
cb.Location       = 'eastoutside';
clim([0, str_lim]);


if exist('num_hotspots','var') && num_hotspots > 0

    lift_um = 0.015 * max([x_max - x_min, y_max - y_min, z_max - z_min]);

    for kk_hot = 1:num_hotspots

        hx0 = h_x(kk_hot);
        hy0 = h_y(kk_hot);
        hz0 = h_z(kk_hot);

        plot3(hx0, hy0, hz0, ...
            'p', ...
            'MarkerSize', 24, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'k', ...
            'LineWidth', 2.5);

        plot3(hx0, hy0, hz0, ...
            'p', ...
            'MarkerSize', 18, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'y', ...
            'LineWidth', 2.0);

        plot3(hx0, hy0, hz0, ...
            'o', ...
            'MarkerSize', 8, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'm', ...
            'LineWidth', 1.2);

        if exist('top_vals','var')
            label_str = sprintf('  #%d | %.1f MPa', kk_hot, top_vals(kk_hot));
        else
            label_str = sprintf('  #%d', kk_hot);
        end

        text(hx0, hy0, hz0 + 1.35*lift_um, label_str, ...
            'FontSize', max(10, font_sz-3), ...
            'FontWeight', 'bold', ...
            'Color', 'k', ...
            'BackgroundColor', [1 1 1 0.80], ...
            'EdgeColor', 'k', ...
            'Margin', 3, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom');
    end
end

axis equal tight;
grid on; box on;
view([-42.71, 22.97]);

hxlab = xlabel('X_1 (μm)');
set(hxlab, 'Rotation', 21);
hylab = ylabel('X_2 (μm)');
set(hylab, 'Rotation', -21);
zlabel('X_3 (μm)');
set(gca, 'FontSize', font_sz, 'LineWidth', 1.2);

xlim([x_min, x_max]);
ylim([y_min, y_max]);
zlim([z_min, z_max + 0.0*lift_um]);

xticks(x_tick);
yticks(y_tick);
zticks(z_tick);
set(gca, 'SortMethod', 'childorder');
hold off;


%% ---------------------- LOCAL FUNCTIONS ---------------------------------
function xnorm = minmax_safe(x)
%MINMAX_SAFE Safely normalize an array to the range [0, 1].
%
% This helper performs standard min-max normalization:
%
%   xnorm = (x - min(x)) / (max(x) - min(x))
%
% but handles the degenerate case where all entries are identical or nearly
% identical.
%
% Why this is needed:
%   In hotspot severity scoring, each score component is normalized before
%   weighting. If all candidate hotspots have the same value for a component,
%   standard min-max normalization would involve division by zero and would
%   produce NaN values or an all-zero contribution.
%
% Behavior:
%   - If x is empty:
%       return x unchanged.
%
%   - If all values are effectively equal:
%       return an array of ones with the same size as x.
%
%   - Otherwise:
%       return standard min-max normalized values in [0, 1].
%
% Interpretation for severity ranking:
%   Returning ones for equal-valued inputs treats all candidates as equally
%   severe with respect to that score component. This makes the component
%   neutral in relative ranking rather than artificially suppressing all scores.
%
% Input:
%   x
%     Numeric array of raw score values.
%
% Output:
%   xnorm
%     Array of normalized score values with the same size as x.

x = double(x);

% Empty input: nothing to normalize.
if isempty(x)
    xnorm = x;
    return;
end

% Compute minimum and maximum score values.
xmin = min(x);
xmax = max(x);

% If the range is effectively zero, all candidates are tied for this metric.
% Use a tolerance scaled by the magnitude of x to avoid numerical noise.
if abs(xmax - xmin) < eps(max(abs(x),[],'all') + 1)
    xnorm = ones(size(x));
else
    xnorm = (x - xmin) ./ (xmax - xmin);
end

end