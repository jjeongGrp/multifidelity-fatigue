%% CPFE-Based Fatigue Life Prediction Using Hotspot Equivalent Strain
% -------------------------------------------------------------------------
% This script estimates fatigue life by combining crystal-plasticity finite
% element, CPFE, hotspot strain results with a strain-life, Coffin-Manson-
% Basquin, CMB, fatigue model.
%
% The script reads equivalent strain fields exported from CPFE simulations,
% identifies a percentile-based hotspot strain inside the common simulation
% domain, computes a strain concentration factor relative to an LPBF baseline
% strain-life curve, and generates a shifted fatigue-life prediction curve.
%
% The workflow is intended for comparing:
%   - an LPBF SS316L strain-life baseline,
%   - a wrought SS316L reference baseline,
%   - experimental fatigue data, and
%   - a CPFE-informed knockdown curve based on hotspot strain amplification.
%
% Main workflow:
%   1. Read CPFE strain-field CSV files for positive and negative loading states.
%   2. Verify that required columns are present in each CSV file.
%   3. Determine a common cubic region of interest, ROI, from the coordinate
%      bounds of all input CSV files.
%   4. Extract equivalent strain values inside the ROI.
%   5. Compute a percentile hotspot equivalent strain for each loading state.
%   6. Average the positive and negative hotspot strains to define a reference
%      hotspot strain amplitude.
%   7. Evaluate the LPBF CMB baseline strain-life curve.
%   8. Compute a CPFE strain concentration factor:
%
%          K_epsilon = epsilon_hotspot / epsilon_nominal
%
%   9. Shift the LPBF baseline curve by K_epsilon to obtain a CPFE-informed
%      nominal strain-life knockdown curve.
%  10. Estimate the predicted life for the reference CPFE hotspot strain.
%  11. Plot LPBF baseline, wrought baseline, experimental data, and CPFE
%      knockdown prediction.
%
% Inputs:
%   - data_export_emax_plus_full_set6_1.csv
%   - data_export_emax_minus_full_set6_1.csv
%
%     CSV files exported from CPFE postprocessing. Each file must contain:
%       * Eqv_strain : equivalent strain field
%       * Points_0   : X coordinate
%       * Points_1   : Y coordinate
%       * Points_2   : Z coordinate
%
% Outputs:
%   - Figure comparing:
%       * ANL wrought SS316L baseline
%       * LPBF SS316L CMB baseline
%       * experimental fatigue data
%       * CPFE-informed knockdown strain-life curve
%
% Key assumptions:
%   - The equivalent strain percentile, p, represents the fatigue-relevant
%     hotspot strain measure.
%   - The CPFE strain concentration factor K_epsilon is constant over the
%     fatigue-life range considered.
%   - The LPBF baseline CMB model provides the nominal strain-life relation
%     before CPFE hotspot correction.
%   - The positive and negative CPFE loading-state hotspot strains are averaged
%     to estimate the reference strain amplitude.
%
% Units:
%   - Strain is dimensionless, mm/mm.
%   - Elastic modulus E and fatigue strength coefficient sigma_f_prime are in MPa.
%   - Cycles to failure are reported as N_f.
%   - Coordinates are assumed to be in the units exported by the CPFE/CSV files.
%
% Notes:
%   - The script uses robust delimiter handling for CSV import and accepts either
%     tab- or comma-delimited files.
%   - The ROI is automatically defined from the coordinate bounds across all CSV
%     files.
%   - The selected hotspot percentile p can be adjusted depending on how
%     conservative the hotspot definition should be.
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
clear; clc;

%% ---------------------- USER INPUT SETTINGS -----------------------------
% CPFE strain-field CSV files to analyze.
%
% Expected interpretation:
%   - The first file corresponds to the positive maximum loading state.
%   - The second file corresponds to the negative maximum loading state.
%
% The script computes a hotspot strain percentile from each file and then
% averages the two hotspot values to estimate the strain amplitude used for
% fatigue-life prediction.
csvFiles = { ...
    'data_export_emax_plus_full_set6_1.csv', ...
    'data_export_emax_minus_full_set6_1.csv'};

% Percentile used to define the hotspot equivalent strain.
%
% Example:
%   p = 95 means the hotspot strain is taken as the 95th percentile of
%   equivalent strain values inside the ROI.
%
% A higher percentile gives a more conservative hotspot measure but is still
% less sensitive to isolated numerical spikes than the absolute maximum.
p = 95.0;

% Number of CPFE strain-field files.
nFiles = numel(csvFiles);

% Required CSV column headers.
%
% Eqv_strain:
%   Equivalent strain value at each exported point/cell/node.
%
% Points_0, Points_1, Points_2:
%   Spatial coordinates of each exported point.
%
% The script preserves original variable names during import, so these names
% must match the CSV headers exactly.
requiredVars = { ...
    'Eqv_strain', ...
    'Points_0', ...
    'Points_1', ...
    'Points_2'};

%% ---------------------- DETERMINE COMMON CUBIC ROI ----------------------
% Determine a common cubic region of interest, ROI, from the coordinate bounds
% of all input CSV files.
%
% The script first finds the global minimum and maximum coordinate values over
% all files. These bounds are then used to define a cube centered at the global
% bounding-box center with side length equal to the largest bounding-box extent.
%
% This ensures that both positive and negative loading-state strain fields are
% evaluated over a consistent spatial region.
xmin_all =  inf; xmax_all = -inf;
ymin_all =  inf; ymax_all = -inf;
zmin_all =  inf; zmax_all = -inf;

for ftmp = 1:nFiles

    csvFile_tmp = csvFiles{ftmp};

    success = false;
    Ttmp = table();

    % First attempt: let MATLAB detect delimiter and import options
    % automatically while preserving original CSV header names.
    try
        opts = detectImportOptions(csvFile_tmp, ...
            'FileType', 'text', ...
            'VariableNamingRule', 'preserve');

        Ttest = readtable(csvFile_tmp, opts);

        if all(ismember(requiredVars, Ttest.Properties.VariableNames))
            Ttmp = Ttest;
            success = true;
        end
    catch
        % If automatic import fails, fall back to explicit delimiter attempts.
    end

    % Second attempt: explicitly try tab- and comma-delimited formats.
    % This makes the script robust to CSV/TSV exports from different tools.
    if ~success
        delimiters = {'\t', ','};

        for idel = 1:numel(delimiters)
            try
                opts = detectImportOptions(csvFile_tmp, ...
                    'FileType', 'text', ...
                    'Delimiter', delimiters{idel}, ...
                    'VariableNamingRule', 'preserve');

                Ttest = readtable(csvFile_tmp, opts);

                if all(ismember(requiredVars, Ttest.Properties.VariableNames))
                    Ttmp = Ttest;
                    success = true;
                    break;
                end
            catch
                % Try the next delimiter.
            end
        end
    end

    if ~success
        error('Could not read required columns from file: %s', csvFile_tmp);
    end

    % Extract coordinate columns.
    % Values may be imported as numeric arrays or strings depending on the CSV
    % format, so each column is converted explicitly to double.
    x_tmp = Ttmp.('Points_0');
    y_tmp = Ttmp.('Points_1');
    z_tmp = Ttmp.('Points_2');

    % Convert coordinates to double precision.
    if isnumeric(x_tmp)
        x_tmp = double(x_tmp);
    else
        x_tmp = str2double(string(x_tmp));
    end

    if isnumeric(y_tmp)
        y_tmp = double(y_tmp);
    else
        y_tmp = str2double(string(y_tmp));
    end

    if isnumeric(z_tmp)
        z_tmp = double(z_tmp);
    else
        z_tmp = str2double(string(z_tmp));
    end

    % Force column-vector format.
    x_tmp = x_tmp(:);
    y_tmp = y_tmp(:);
    z_tmp = z_tmp(:);

    % Keep only rows with finite coordinate values.
    good_tmp = isfinite(x_tmp) & isfinite(y_tmp) & isfinite(z_tmp);

    x_tmp = x_tmp(good_tmp);
    y_tmp = y_tmp(good_tmp);
    z_tmp = z_tmp(good_tmp);

    if isempty(x_tmp)
        error('No valid coordinate data found in %s.', csvFile_tmp);
    end

    % Update global coordinate bounds over all files.
    xmin_all = min(xmin_all, min(x_tmp));
    xmax_all = max(xmax_all, max(x_tmp));

    ymin_all = min(ymin_all, min(y_tmp));
    ymax_all = max(ymax_all, max(y_tmp));

    zmin_all = min(zmin_all, min(z_tmp));
    zmax_all = max(zmax_all, max(z_tmp));

end

%% ---------------------- DEFINE CUBIC ROI --------------------------------
% Define a cubic region of interest from the global coordinate bounds found
% across all CPFE CSV files.
%
% The initial coordinate extents are:
%   Lx = xmax_all - xmin_all
%   Ly = ymax_all - ymin_all
%   Lz = zmax_all - zmin_all
%
% The ROI is chosen as a cube centered at the bounding-box center:
%
%   (x0, y0, z0)
%
% with side length:
%
%   L = max([Lx, Ly, Lz])
%
% and half-width:
%
%   h = L/2
%
% This gives a common cubic ROI large enough to contain the coordinate extents
% of all input files.
Lx = xmax_all - xmin_all;
Ly = ymax_all - ymin_all;
Lz = zmax_all - zmin_all;

% Center of the global coordinate bounding box.
x0 = 0.5 * (xmin_all + xmax_all);
y0 = 0.5 * (ymin_all + ymax_all);
z0 = 0.5 * (zmin_all + zmax_all);

% Cubic ROI side length and half-width.
L = max([Lx, Ly, Lz]);
h = L / 2;

% Numerical tolerance for ROI membership checks.
% The tolerance scales with the ROI size to avoid excluding points due to
% floating-point roundoff near cube boundaries.
tolROI = 1e-10 * max(1, L);

%% ---------------------- COMPUTE HOTSPOT STRAIN VALUES -------------------
% Compute one percentile hotspot equivalent strain value from each CPFE CSV
% file.
%
% For each file:
%   1. Read the table robustly using automatic delimiter detection, with
%      tab/comma fallbacks.
%   2. Extract equivalent strain and coordinate columns.
%   3. Remove invalid rows.
%   4. Apply the common cubic ROI mask.
%   5. Compute the p-th percentile of Eqv_strain inside the ROI.
%
% Output:
%   eps_hot_all(f) = p-th percentile equivalent strain for csvFiles{f}
eps_hot_all = zeros(1, nFiles);

for f = 1:nFiles

    csvFile = csvFiles{f};

    success = false;
    T = table();

    % First attempt: automatic delimiter/import-option detection.
    try
        opts = detectImportOptions(csvFile, ...
            'FileType', 'text', ...
            'VariableNamingRule', 'preserve');

        Ttest = readtable(csvFile, opts);

        % Accept import only if required columns are present.
        if all(ismember(requiredVars, Ttest.Properties.VariableNames))
            T = Ttest;
            success = true;
        end
    catch
        % Fall back to explicit delimiter attempts below.
    end

    % Second attempt: explicitly try tab and comma delimiters.
    if ~success
        delimiters = {'\t', ','};

        for idel = 1:numel(delimiters)
            try
                opts = detectImportOptions(csvFile, ...
                    'FileType', 'text', ...
                    'Delimiter', delimiters{idel}, ...
                    'VariableNamingRule', 'preserve');

                Ttest = readtable(csvFile, opts);

                if all(ismember(requiredVars, Ttest.Properties.VariableNames))
                    T = Ttest;
                    success = true;
                    break;
                end
            catch
                % Try the next delimiter.
            end
        end
    end

    if ~success
        error('Could not read required columns from file: %s', csvFile);
    end

    % Extract equivalent strain and coordinates.
    % Columns may be imported as numeric arrays or text, so conversion to
    % double is handled explicitly below.
    eps_raw = T.('Eqv_strain');
    x_raw   = T.('Points_0');
    y_raw   = T.('Points_1');
    z_raw   = T.('Points_2');

    % Convert strain and coordinate columns to double precision.
    if isnumeric(eps_raw)
        eps_raw = double(eps_raw);
    else
        eps_raw = str2double(string(eps_raw));
    end

    if isnumeric(x_raw)
        x_raw = double(x_raw);
    else
        x_raw = str2double(string(x_raw));
    end

    if isnumeric(y_raw)
        y_raw = double(y_raw);
    else
        y_raw = str2double(string(y_raw));
    end

    if isnumeric(z_raw)
        z_raw = double(z_raw);
    else
        z_raw = str2double(string(z_raw));
    end

    % Force column-vector format.
    eps_raw = eps_raw(:);
    x_raw   = x_raw(:);
    y_raw   = y_raw(:);
    z_raw   = z_raw(:);

    % Remove rows with invalid strain or coordinate values.
    good = isfinite(eps_raw) & ...
           isfinite(x_raw)   & ...
           isfinite(y_raw)   & ...
           isfinite(z_raw);

    eps = eps_raw(good);
    x   = x_raw(good);
    y   = y_raw(good);
    z   = z_raw(good);

    if isempty(eps)
        error('No valid strain/coordinate rows found in %s.', csvFile);
    end

    % Apply the common cubic ROI mask.
    %
    % A point is inside the ROI if it lies within h of the ROI center in all
    % three coordinate directions.
    ROI = abs(x - x0) <= h + tolROI & ...
          abs(y - y0) <= h + tolROI & ...
          abs(z - z0) <= h + tolROI;

    if ~any(ROI)
        error('ROI empty for %s.', csvFile);
    end

    % Compute the percentile hotspot equivalent strain inside the ROI.
    eps_roi = eps(ROI);

    eps_hot_all(f) = prctile(eps_roi, p);

end

%% ---------------------- STRAIN-LIFE MODEL SETUP -------------------------
% Define strain-life curves used for fatigue-life prediction and comparison.
%
% This section sets up:
%   1. An LPBF SS316L baseline strain-life curve using the
%      Coffin-Manson-Basquin, CMB, relation.
%   2. A single experimental fatigue data point from DARPA Set-6 LOF.
%   3. A wrought SS316L baseline curve from ANL for reference.
%
% The CPFE-informed knockdown curve is computed in the next section by applying
% a hotspot strain concentration factor to the LPBF baseline.
%% ---------------------- MATERIAL PROPERTY --------------------------------
% Elastic modulus used in the elastic Basquin term of the CMB equation.
%
% Units:
%   E is in MPa.
E       = 225.4507e3;

%% ---------------------- LPBF CMB PARAMETERS ------------------------------
% Coffin-Manson-Basquin, CMB, parameters for the LPBF SS316L baseline.
%
% CMB strain-amplitude relation:
%
%   epsilon_a = (sigma_f' / E) * (2N_f)^b + epsilon_f' * (2N_f)^c
%
% where:
%   epsilon_a     = strain amplitude
%   N_f           = cycles to failure
%   2N_f          = reversals to failure
%   sigma_f'      = fatigue strength coefficient, MPa
%   b             = fatigue strength exponent
%   epsilon_f'    = fatigue ductility coefficient
%   c             = fatigue ductility exponent
%
% These parameters are fitted from LPBF experimental fatigue data.
sigma_f_prime = 690.2592;
b_exp         = -0.0200;
eps_f_prime   = 0.4423;
c_exp         = -0.5321;

%% ---------------------- LPBF BASELINE CURVE ------------------------------
% Define the CMB strain-amplitude function.
%
% Parameter vector p:
%   p(1) = sigma_f'
%   p(2) = b
%   p(3) = epsilon_f'
%   p(4) = c
cmb_epsa = @(p, Nf) (p(1)/E).*(2.*Nf).^p(2) + p(3).*(2.*Nf).^p(4);

% LPBF baseline parameter vector.
p_lpbf   = [sigma_f_prime, b_exp, eps_f_prime, c_exp];

% Cycles-to-failure range used for plotting and interpolation.
Nf_plot       = logspace(2, 7, 300);

% Smooth LPBF baseline strain-amplitude curve.
epsa_baseline = cmb_epsa(p_lpbf, Nf_plot);

%% ---------------------- EXPERIMENTAL DATA POINT --------------------------
% Experimental fatigue data used for comparison.
%
% DARPA Set-6 LOF:
%   Nf_set6_exp : measured cycles to failure
%   A_set6_exp  : applied nominal strain amplitude, mm/mm
Nf_set6_exp       = [19087];
A_set6_exp        = [0.003];

% Horizontal error bar for experimental life.
% Here, ±30% life scatter is shown for visualization.
error_bar_percent = 0.30;
error_bar_cycles  = error_bar_percent * Nf_set6_exp;

%% ---------------------- ANL WROUGHT SS316L BASELINE ----------------------
% Wrought SS316L reference baseline from ANL.
%
% The expression returns strain amplitude in percent, so it is divided by 100
% to convert to strain, mm/mm.
ANL_Nf_plot  = logspace(2, 7, 100);
ANL_baseline = (23.0 .* ANL_Nf_plot.^(-0.457) + 0.11) ./ 100;

%% ---------------------- CPFE-INFORMED LIFE SHIFT -------------------------
% Use the CPFE hotspot strain to shift the LPBF baseline strain-life curve.
%
% Concept:
%   The LPBF CMB curve gives nominal strain amplitude as a function of fatigue
%   life:
%
%       epsilon_nominal = f(N_f)
%
%   The CPFE simulation provides a hotspot equivalent strain amplitude for a
%   reference nominal loading condition.
%
%   The ratio between hotspot strain and nominal baseline strain defines a
%   strain concentration factor:
%
%       K_epsilon = epsilon_hotspot / epsilon_nominal
%
%   Assuming K_epsilon is approximately constant over the life range, the
%   nominal strain-life curve is shifted downward:
%
%       epsilon_nominal_CPFE(N_f) = epsilon_baseline(N_f) / K_epsilon
%
%   This gives a CPFE-informed knockdown curve for nominal strain amplitude.

%% ---- Select reference nominal strain from LPBF baseline -----------------
% Choose a reference nominal strain amplitude corresponding to the experimental
% loading level.
%
% eps_nom_target:
%   Target nominal strain amplitude, here taken from the experimental Set-6
%   loading amplitude.
%
% eps_nom_ref:
%   Closest strain amplitude on the discretized LPBF baseline curve.
%
% Nf_ref:
%   Fatigue life on the LPBF baseline corresponding to eps_nom_ref.
eps_nom_target = A_set6_exp(1);
[~, i_ref]     = min(abs(epsa_baseline - eps_nom_target));
eps_nom_ref    = epsa_baseline(i_ref);
Nf_ref         = Nf_plot(i_ref);

%% ---- Compute reference CPFE hotspot strain amplitude --------------------
% Extract percentile hotspot equivalent strains from the two CPFE loading
% states.
%
% Assumption:
%   csvFiles{1} corresponds to the positive maximum loading state.
%   csvFiles{2} corresponds to the negative maximum loading state.
%
% The reference hotspot strain amplitude is estimated by averaging the two
% percentile hotspot values.
eps_a_hot_ref_plus = eps_hot_all(1);
eps_a_hot_ref_minus = eps_hot_all(2);
eps_a_hot_ref = 0.5*(eps_a_hot_ref_plus+eps_a_hot_ref_minus);

%% ---- Compute strain concentration factor --------------------------------
% Compute CPFE strain concentration factor.
%
% K_epsilon > 1 means the local hotspot strain is higher than the nominal
% baseline strain amplitude.
%
% The lower bound of 1.0 prevents the CPFE correction from increasing life if
% the computed hotspot strain is less than the nominal reference strain.
Keps = eps_a_hot_ref / eps_nom_ref;
Keps = max(Keps, 1.0);

%% ---- Generate CPFE knockdown strain-life curve --------------------------
% Apply the CPFE strain concentration factor to the full LPBF baseline curve.
%
% For a given fatigue life N_f, the allowable nominal strain amplitude is
% reduced by K_epsilon.
Nf_cp_curve = Nf_plot;
epsa_cp_nom = epsa_baseline ./ Keps;

%% ---- Predict life for the reference CPFE case ---------------------------
% Estimate the fatigue life corresponding to the CPFE hotspot strain.
%
% The baseline curve gives:
%
%   epsilon_a = f(N_f)
%
% To estimate life from a strain amplitude, invert this relation by
% interpolation. Because epsa_baseline decreases with N_f, the strain array is
% sorted in ascending order before calling interp1.
%
% Nf_pred_ref:
%   Predicted cycles to failure corresponding to eps_a_hot_ref based on the
%   LPBF baseline curve.
A = epsa_baseline(:);
N = Nf_plot(:);

% Sort strain amplitudes in ascending order for interpolation.
[A_inc, idx] = sort(A, 'ascend');  
N_inc = N(idx);

% Remove duplicate strain values, if any, to satisfy interp1 requirements.
[A_inc, ia]  = unique(A_inc, 'stable');  
N_inc = N_inc(ia);

% Interpolate predicted life at the CPFE hotspot strain.
Nf_pred_ref = interp1(A_inc, N_inc, eps_a_hot_ref, 'pchip', 'extrap');

%% ---------------------- PLOT STRAIN-LIFE COMPARISON ----------------------
% Plot the fatigue strain-life curves and reference data:
%
%   1. ANL wrought SS316L reference baseline
%   2. LPBF SS316L CMB baseline
%   3. Experimental Set-6 LOF fatigue data point
%   4. CPFE-informed knockdown curve
%   5. CPFE-predicted life marker for the reference loading condition
figure; clf; hold on;

% Plot wrought SS316L reference baseline.
loglog(ANL_Nf_plot, ANL_baseline, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, ...
    'DisplayName', 'ANL Baseline (Wrought SS316L)');

% Plot LPBF CMB baseline fitted from LPBF fatigue data.
loglog(Nf_plot, epsa_baseline,    'k-',  'LineWidth', 2, ...
    'DisplayName','LPBF Baseline');

% Plot experimental fatigue data.
%
% Horizontal error bars represent ±30% scatter in fatigue life.
p_exp = errorbar(Nf_set6_exp, A_set6_exp, [], [], ...
    error_bar_cycles, error_bar_cycles, ...
    's', 'MarkerSize', 11, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2, 'Color', 'k', 'CapSize', 8, ...
    'DisplayName', 'Exp: Set6-1 (±30%)');

% Annotate the experimental fatigue life.
text(Nf_set6_exp, A_set6_exp * 1.10, ...
    sprintf('Exp N_f \n%.0f cycles', Nf_set6_exp), ...
    'Color', 'k', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');

% Plot the CPFE-informed knockdown curve.
%
% This curve is obtained by dividing the LPBF baseline strain amplitude by
% K_epsilon. It represents the nominal strain amplitude that would produce the
% same local hotspot strain demand.
loglog(Nf_cp_curve, epsa_cp_nom, 'b-', 'LineWidth', 2, ...
    'DisplayName', sprintf('CPFE knockdown (K_\\epsilon=%.2f)', Keps));

% Plot the CPFE-predicted life marker.
%
% The x-coordinate, Nf_pred_ref, is obtained by evaluating the LPBF baseline
% life corresponding to the CPFE hotspot strain.
%
% The y-coordinate is the nominal reference strain amplitude, eps_nom_ref, so
% the marker lies on the nominal strain-life plot.
loglog(Nf_pred_ref, eps_nom_ref, 'bx', 'MarkerFaceColor','b', ...
    'MarkerSize', 12,'LineWidth', 2, 'HandleVisibility','off');

% Annotate the CPFE-predicted life.
text(Nf_pred_ref, eps_nom_ref * 1.15, ...
    sprintf('CP N_f \n%d cycles', round(Nf_pred_ref)), ...
    'Color', 'b', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');


grid on; grid minor;
xlim([1e2, 1e7]);
ylim([0, 0.016])
set(gca,'XScale','log','YScale','linear');
xlabel('Cycles to Failure, N_f', 'FontSize', 13, 'FontWeight','bold');
ylabel('Strain Amplitude, \epsilon_a  (mm/mm)', ...
    'FontSize', 13, 'FontWeight', 'bold');
legend('Location','northeast');

hold off
set(gca,'FontSize',15);
