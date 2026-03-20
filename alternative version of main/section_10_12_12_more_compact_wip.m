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