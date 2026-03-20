% =========================================================================
%  FAVAR – Monetary Policy Transmission in the Euro Area
%  Author: Lorenzo Pisa
% =========================================================================

clear all; close all;
addpath(genpath(pwd))

% -------------------------------------------------------------------------
%  SETTINGS
% -------------------------------------------------------------------------
confidence       = [5 95];   % bootstrap confidence band (90%)
n_bootstrap      = 1000;
H                = 60;       % IRF horizon (months)
min_var_explained = 0.60;    % factor selection threshold

forcing  = 1;   % 1 = force common lag across all models
forced_p = 7;   % ignored if forcing = 0 (AIC used instead)
p_max    = 12;

tic

% =========================================================================
%  1. DATA LOADING
% =========================================================================

varNames = string(readcell("data/EAdata2.xlsx", "Range", "B1:DO1"));
X_raw    = readmatrix("data/EAdata2.xlsx",       "Range", "B3:DO313");
t        = datetime(2000,1,1) : calmonths(1) : datetime(2025,11,1);

% Standardize
X_std          = zscore(X_raw);
varNames_baseline = varNames;
X_std_baseline    = X_std;

% Standard deviations for IRF rescaling
sigma          = std(X_raw, 'omitnan');
sigma_baseline = sigma;

% =========================================================================
%  2. EXPLICIT VARIABLES (Y block: GDP, HICP, Policy Rate)
% =========================================================================

% Column indices in the original dataset
idx_gdp  = 1;
idx_hicp = 98;
idx_rate = 75;

Y       = X_std(:, [idx_gdp, idx_hicp, idx_rate]);
Y       = Y - mean(Y, 'omitnan');
sigma_Y = sigma([idx_gdp, idx_hicp, idx_rate]);

% Remove explicit variables from X (they enter Z directly)
idx_remove          = [idx_gdp, idx_rate, idx_hicp];
idx_remove_baseline = idx_rate;

varNames(idx_remove)           = [];
varNames_baseline(idx_remove_baseline) = [];
X_std(:, idx_remove)           = [];
X_std_baseline(:, idx_remove_baseline) = [];
sigma(idx_remove)              = [];
sigma_baseline(idx_remove_baseline) = [];

% =========================================================================
%  3. FAST / SLOW VARIABLE CLASSIFICATION
% =========================================================================

fast_keywords = ["ERUS","REER","IRT","SHIX","BCI","CCI","SENTIX", ...
                 "TASS","TLB","M1_","M2_","CURR"];

is_fast          = false(1, size(varNames, 2));
is_fast_baseline = false(1, size(varNames_baseline, 2));

for j = 1:length(fast_keywords)
    is_fast          = is_fast          | contains(varNames,          fast_keywords(j));
    is_fast_baseline = is_fast_baseline | contains(varNames_baseline, fast_keywords(j));
end

X_slow          = X_std(:, ~is_fast);
X_slow_baseline = X_std_baseline(:, ~is_fast_baseline);

% =========================================================================
%  4. SHADOW RATE
% =========================================================================

shadowrate = cell2mat(struct2cell(load("data/shadowrate_spliced.mat")));
shadowrate = zscore(shadowrate);

% =========================================================================
%  5. FACTOR SELECTION
% =========================================================================

r_opt          = optimal_r(X_std,          min_var_explained, "FAVAR");
r_opt_baseline = optimal_r(X_std_baseline, min_var_explained, "FAVAR Baseline");

policy_rate = Y(:, 3);

% Number of latent factors (explicit variables are subtracted)
r          = r_opt - size(Y, 2);
r_baseline = r_opt_baseline - 1;

% =========================================================================
%  6. FACTOR EXTRACTION – BERNANKE CLEANING
% =========================================================================

Fhat                      = bernanke_cleaning(r,          X_std,          X_slow,          policy_rate);
Fhat_baseline             = bernanke_cleaning(r_baseline, X_std_baseline, X_slow_baseline, policy_rate);
Fhat_shadow               = bernanke_cleaning(r,          X_std,          X_slow,          shadowrate);
Fhat_baseline_shadow      = bernanke_cleaning(r_baseline, X_std_baseline, X_slow_baseline, shadowrate);

% =========================================================================
%  7. STATE VECTORS
% =========================================================================

Z                    = [Fhat,               Y];
Z_shadow             = [Fhat_shadow,        Y(:,[1,2]), shadowrate];
Z_baseline           = [Fhat_baseline,      policy_rate];
Z_baseline_shadow    = [Fhat_baseline_shadow, shadowrate];

% =========================================================================
%  8. LAG SELECTION
% =========================================================================

[optimal_p, ~, res_optimalp] = VarOLSbestP(Z, p_max);

if forcing && forced_p >= 1
    optimal_p = forced_p;
    [~, res_optimalp] = VarOLS(Z, optimal_p);
end

% Force same lag across all models (or re-select per model)
if forcing
    optimal_p_shadow          = optimal_p;
    optimal_p_baseline        = optimal_p;
    optimal_p_baseline_shadow = optimal_p;
    optimal_p_var             = optimal_p;
else
    [optimal_p_baseline,        ~, ~] = VarOLSbestP(Z_baseline,        p_max);
    [optimal_p_baseline_shadow, ~, ~] = VarOLSbestP(Z_baseline_shadow, p_max);
    [optimal_p_shadow,          ~, ~] = VarOLSbestP(Z_shadow,          p_max);
    [optimal_p_var,             ~, ~] = VarOLSbestP(Y,                 p_max);
end

% Residual autocorrelation check
autocorr_labels = [string(1:r), "GDP", "HICP", "Policy Rate"];
figure
for i = 1:size(Z, 2)
    subplot(3, 3, i)
    autocorr(res_optimalp(:, i))
    title(autocorr_labels(i))
end

% =========================================================================
%  9. IDENTIFICATION – CHOLESKY
% =========================================================================

constant   = 1;
cum_index  = [];
display_on = 0;

[favar_irf,          ~, favar_irf_boot,          ~, ~, ~] = CholeskyIdentification(Z,                 optimal_p,                H, constant, n_bootstrap, cum_index, [string(1:r), "GDP","HICP","Policy Rate"],    display_on);
[var_irf,            ~, var_irf_boot,            ~, ~, ~] = CholeskyIdentification(Y,                 optimal_p_var,            H, constant, n_bootstrap, cum_index, ["GDP","HICP","Policy Rate"],                  display_on);
[favar_baseline_irf, ~, favar_baseline_irf_boot, ~, ~, ~] = CholeskyIdentification(Z_baseline,        optimal_p_baseline,       H, constant, n_bootstrap, cum_index, [string(1:r_baseline), "Policy Rate"],         display_on);
[favar_shadow_irf,   ~, favar_shadow_irf_boot,   ~, ~, ~] = CholeskyIdentification(Z_shadow,          optimal_p_shadow,         H, constant, n_bootstrap, cum_index, [string(1:r), "GDP","HICP","Shadow Rate"],     display_on);
[favar_bl_sh_irf,    ~, favar_bl_sh_irf_boot,    ~, ~, ~] = CholeskyIdentification(Z_baseline_shadow, optimal_p_baseline_shadow, H, constant, n_bootstrap, cum_index, [string(1:r_baseline), "Shadow Rate"],        display_on);

var_labels = ["GDP", "HICP", "Policy Rate"];

% =========================================================================
%  10. IRF LOADING (factor IRFs → observed variables)
% =========================================================================

idx_slow          = find(~is_fast);
idx_slow_baseline = find(~is_fast_baseline);

[IRF_favar,       IRF_favar_boot,       ~] = irf_loading(X_std,          Z,                sigma,          favar_irf,       favar_irf_boot,       idx_slow,          r + 3);
[IRF_shadow,      IRF_shadow_boot,      ~] = irf_loading(X_std,          Z_shadow,         sigma,          favar_shadow_irf, favar_shadow_irf_boot, idx_slow,          r + 1);
[IRF_baseline,    IRF_baseline_boot,    ~] = irf_loading(X_std_baseline, Z_baseline,       sigma_baseline, favar_baseline_irf, favar_baseline_irf_boot, idx_slow_baseline, r_baseline + 1);
[IRF_bl_sh,       IRF_bl_sh_boot,       ~] = irf_loading(X_std_baseline, Z_baseline_shadow, sigma_baseline, favar_bl_sh_irf, favar_bl_sh_irf_boot, idx_slow_baseline, r_baseline + 1);

% =========================================================================
%  11. CUMULATIVE IRFs
% =========================================================================

[cum_var_irf,      cum_var_irf_boot]      = cumulative_irf(var_irf,          var_irf_boot,          [1, 2]);
[cum_favar_irf,    cum_favar_irf_boot]    = cumulative_irf(favar_irf,        favar_irf_boot,        [r+1, r+2]);
[cum_shadow_irf,   cum_shadow_irf_boot]   = cumulative_irf(favar_shadow_irf, favar_shadow_irf_boot, [r+1, r+2]);

% Baseline: extract GDP + HICP from loadings, rate from VAR block
baseline_irf     = [IRF_baseline([1,97],    r_baseline+1, :);   favar_baseline_irf(r_baseline+1,  r_baseline+1, :)];
baseline_irf_boot= [IRF_baseline_boot([1,97], r_baseline+1, :, :); favar_baseline_irf_boot(r_baseline+1, r_baseline+1, :, :)];
[cum_bl_irf,   cum_bl_irf_boot]   = cumulative_irf(baseline_irf,      baseline_irf_boot,      [1, 2]);

bl_sh_irf      = [IRF_bl_sh([1,97],      r_baseline+1, :);   favar_bl_sh_irf(r_baseline+1,    r_baseline+1, :)];
bl_sh_irf_boot = [IRF_bl_sh_boot([1,97], r_baseline+1, :, :); favar_bl_sh_irf_boot(r_baseline+1, r_baseline+1, :, :)];
[cum_bl_sh_irf, cum_bl_sh_irf_boot] = cumulative_irf(bl_sh_irf, bl_sh_irf_boot, [1, 2]);

% =========================================================================
%  12. MODEL COMPARISON PLOT (VAR vs all FAVAR variants)
% =========================================================================

% Collect policy-rate shock IRFs for GDP, HICP, Rate across models
irf_s      = {cum_var_irf(:,3,:), ...
              cum_favar_irf(r+1:r+3, r+3, :), ...
              cum_bl_irf, ...
              cum_shadow_irf(r+1:r+3, r+3, :), ...
              cum_bl_sh_irf};

irf_boot_s = {cum_var_irf_boot(:,3,:,:), ...
              cum_favar_irf_boot(r+1:r+3, r+3, :, :), ...
              cum_bl_irf_boot, ...
              cum_shadow_irf_boot(r+1:r+3, r+3, :, :), ...
              cum_bl_sh_irf_boot};

% Rescale to original units (×sigma) and percentage points (×100 for GDP/HICP)
scale_idx_full     = {1:3, 1:3, 1:2, 1:3, 1:2};   % rows to rescale with sigma_Y
scale_idx_x100     = {[1,2], [1,2], [1,2], [1,2], [1,2]};
sigma_Y_map        = {sigma_Y', sigma_Y', sigma_Y(3)', sigma_Y', sigma_Y(3)'};

for m = 1:5
    rows = scale_idx_full{m};
    irf_s{m}(rows, 1, :)          = irf_s{m}(rows, 1, :)          .* sigma_Y_map{m}(rows);
    irf_boot_s{m}(rows, 1, :, :)  = irf_boot_s{m}(rows, 1, :, :)  .* sigma_Y_map{m}(rows);
    r100 = scale_idx_x100{m};
    irf_s{m}(r100, 1, :)         = irf_s{m}(r100, 1, :)         * 100;
    irf_boot_s{m}(r100, 1, :, :) = irf_boot_s{m}(r100, 1, :, :) * 100;
end

model_legends = ["VAR", "FAVAR", "FAVAR Baseline", "FAVAR Shadow", "FAVAR Baseline Shadow"];
compare_irf(irf_s, irf_boot_s, var_labels, "Policy Rate Shock", model_legends, 0, confidence);

% =========================================================================
%  13. ROBUSTNESS: IRFs ACROSS DIFFERENT NUMBER OF FACTORS
% =========================================================================

irf_s_r      = {};
irf_boot_s_r = {};

for i = 1:r
    F_r   = bernanke_cleaning(i, X_std, X_slow, policy_rate);
    Z_r   = [F_r, Y];
    [irf_r, ~, irf_boot_r, ~, ~, ~] = CholeskyIdentification(Z_r, optimal_p, H, constant, n_bootstrap, [], "", 0);
    [cum_r, cum_boot_r] = cumulative_irf(irf_r, irf_boot_r, [i+1, i+2]);

    irf_s_r{i}      = cum_r(i+1:i+3,    i+3, :)    .* sigma_Y';
    irf_boot_s_r{i} = cum_boot_r(i+1:i+3, i+3, :, :) .* sigma_Y';
    irf_s_r{i}([1,2], :, :)         = irf_s_r{i}([1,2], :, :)         * 100;
    irf_boot_s_r{i}([1,2], :, :, :) = irf_boot_s_r{i}([1,2], :, :, :) * 100;
end

compare_irf(irf_s_r, irf_boot_s_r, var_labels, ...
    "FAVAR – Robustness: varying number of factors", ...
    [string(1+3 : r+3-1), "Optimal r"], 0, confidence);

% =========================================================================
%  14. MAIN IRF DISPLAY (loaded variables)
% =========================================================================

target_var = ["LTIRT_EACC","M1_EACC","HFCE_EA","M2_EACC","CONSD_EA", ...
              "ERUS_EA","CONSND_EA","UNETOT_EA","HICPIN_EA","SHIX_EA","TEMP_EA","CCI_EA"];
idx_irfs   = get_indices(varNames, target_var);

% Labels: rate shock + explicit vars + loaded vars (strip suffixes)
labels_disp = ["IRT3M: shock", "REAL GDP", "HICP", target_var];
labels_disp = replace(labels_disp, "_EACC", "");
labels_disp = replace(labels_disp, "_EA",   "");

% Stack explicit and loaded IRFs
disp_irf      = [favar_irf(     [r+3,r+1,r+2], :, :)    .* sigma_Y([3,1,2])'; IRF_favar(    idx_irfs, :, :)];
disp_irf_boot = [favar_irf_boot([r+3,r+1,r+2], :, :, :) .* sigma_Y([3,1,2])'; IRF_favar_boot(idx_irfs, :, :, :)];

% Cumulate level variables and convert to percentage points
cum_idx = [2,3,5,6,7,8,9,10,12,13,14];
labels_disp(cum_idx) = labels_disp(cum_idx) + ": level";
[cum_disp_irf, cum_disp_irf_boot] = cumulative_irf(disp_irf, disp_irf_boot, cum_idx);
cum_disp_irf(     cum_idx, :, :)       = cum_disp_irf(     cum_idx, :, :)       * 100;
cum_disp_irf_boot(cum_idx, :, :, :)    = cum_disp_irf_boot(cum_idx, :, :, :)    * 100;

% Plot
x = (1:H)';
figure
for i = 1:length(labels_disp)
    subplot(4, 4, i)
    band  = squeeze(prctile(cum_disp_irf_boot(i, r+3, :, :), confidence, 4));
    lower = band(:, 1);
    upper = band(:, 2);
    fill([x; flipud(x)], [upper; flipud(lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    hold on
    plot(1:H, squeeze(cum_disp_irf(i, r+3, :)), 'b')
    yline(0, 'r--', 'HandleVisibility', 'off')
    axis tight; hold off
    title(labels_disp(i))
end

% =========================================================================
%  15. DIAGNOSTICS: FEVD AND R²
% =========================================================================

FEVD_result          = calc_favar_fevd(disp_irf);
R2_vector            = calc_favar_r2(X_std, Z);

Variance_Decomp      = FEVD_result(:, r+3, 60);
R2                   = [1, 1, 1, R2_vector(idx_irfs)];

% Summary table
fig = uifigure;
uitable(fig, ...
    "Data",       [labels_disp', Variance_Decomp, R2'], ...
    "ColumnName", ["Variable", "FEVD at h=60", "R²"]);

toc