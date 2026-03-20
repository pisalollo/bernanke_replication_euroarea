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


%% ============================================================
%  EXPLICIT VARIABLES (Y block)
% ============================================================

% Explicit macro variables:
% 1 = GDP
% 98 = HICP overall
% 75 = 3M policy rate
Y = X_std(:,[1,98,75]);

% Remove sample mean
Y = Y - mean(Y,'omitnan');

% Standard deviations of explicit variables
sigma_Y = sigma([1,98,75]);

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

% Mark fast variables based on keywords
for j = 1:length(fast_keywords)
    is_fast          = is_fast          | contains(varNames,          fast_keywords(j));
    is_fast_baseline = is_fast_baseline | contains(varNames_baseline, fast_keywords(j));
end

% Slow-moving variables only (used for Bernanke restriction)
X_slow = X_std(:,~is_fast);
X_slow_baseline = X_std_baseline(:,~is_fast_baseline);

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

% Factors with policy rate
Fhat = bernanke_cleaning(r,X_std,X_slow,policy_rate);
Fhat_baseline = bernanke_cleaning(r_baseline,X_std_baseline,X_slow_baseline,policy_rate);
% Factors with shadow rate
Fhat_shadow = bernanke_cleaning(r,X_std,X_slow,shadowrate);
Fhat_baseline_shadow = bernanke_cleaning(r_baseline,X_std_baseline,X_slow_baseline,shadowrate);

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

[optimal_p, ~, res_optimalp, ~, ~]= VarOLSbestP(Z,p_max);

if forcing && forced_p >= 1
    optimal_p = forced_p;
    [~, res_optimalp] = VarOLS(Z, optimal_p);
end

% Force same lag across all models (or re-select per model)
if forcing
    optimal_p_shadowrate=optimal_p;
    optimal_p_baseline=optimal_p;
    optimal_p_baseline_shadowrate=optimal_p;
    optimal_p_var=optimal_p;
else %robustness test on different p lags for various models
    [optimal_p_baseline, ~, ~, ~, ~]= VarOLSbestP(Z_baseline,p_max);
    [optimal_p_baseline_shadowrate, ~, ~, ~, ~]= VarOLSbestP(Z_baseline_shadow,p_max);
    [optimal_p_shadowrate, ~, ~, ~, ~]= VarOLSbestP(Z_shadow,p_max);
    [optimal_p_var, ~, ~, ~, ~]= VarOLSbestP(Y,p_max);
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


constant = 1;
cum_index = [];
display_on = 0;

% Main FAVAR
favar_labels=[1:r, "GDP","HICP","Policy Rate"];
[favar_irf,~,favar_irf_boot,~,~,~] = CholeskyIdentification(Z,optimal_p,H,constant,n_bootstrap,cum_index,favar_labels,display_on);

% Standard VAR
var_labels=["GDP","HICP","Policy Rate"];
[var_irf,~,var_irf_boot,~,~,~] = CholeskyIdentification(Y,optimal_p_var,H,constant,n_bootstrap,cum_index,var_labels,display_on);

% Baseline FAVAR
favar_baseline_labels=[1:r_baseline,"Policy Rate"];
[favar_baseline_irf,~,favar_baseline_irf_boot,~,~,~] = CholeskyIdentification(Z_baseline,optimal_p_baseline,H,constant,n_bootstrap,cum_index,favar_baseline_labels,display_on);

% Shadow rate Baseline
favar_baseline_shadow_labels=[1:r_baseline,"Shadow Rate"];
[favar_bs_sh_irf,~,favar_bl_sh_irf_boot,~,~,~] = CholeskyIdentification(Z_baseline_shadow,optimal_p_baseline_shadowrate,H,constant,n_bootstrap,cum_index,favar_baseline_shadow_labels,display_on);

% Shadow rate (Main variant)
favar_shadow_labels=[1:r,"GDP","HICP","Shadow Rate"];
[favar_shadow_irf,~,favar_shadow_irf_boot,~,~,~] = CholeskyIdentification(Z_shadow,optimal_p_shadowrate,H,constant,n_bootstrap,cum_index,favar_shadow_labels,display_on);

%% ============================================================
%  LOADINGS: MAP FACTOR IRFs TO OBSERVED VARIABLES
%% ============================================================

index_favar_slow = find(~is_fast);
[IRF_loaded_favar, IRF_loaded_boot_favar,~]=irf_loading(X_std,Z,sigma,favar_irf,favar_irf_boot,index_favar_slow,r+3);
[IRF_loaded_favar_shadowrate, IRF_loaded_boot_favar_shadowrate,~]=irf_loading(X_std,Z_shadow,sigma,favar_shadow_irf,favar_shadow_irf_boot,index_favar_slow,r+1);

index_baseline_slow = find(~is_fast_baseline);
[IRF_loaded_favar_baseline, IRF_loaded_boot_favar_baseline,~]=irf_loading(X_std_baseline,Z_baseline,sigma_baseline,favar_baseline_irf,favar_baseline_irf_boot,index_baseline_slow,r_baseline+1);
[IRF_loaded_favar_baseline_shadowrate, IRF_loaded_boot_favar_baseline_shadowrate,~]=irf_loading(X_std_baseline,Z_baseline_shadow,sigma_baseline,favar_bs_sh_irf,favar_bl_sh_irf_boot,index_baseline_slow,r_baseline+1);

%% ============================================================
%  CUMULATIVE IMPULSE RESPONSES
%% ============================================================

%var con gdp, hicp, policy rate
cum_index=[1,2];
[cumulative_var_irf, cumulative_var_irf_boot]=cumulative_irf(var_irf,var_irf_boot,cum_index);

%favar con r ottimo con gdp, hicp, policy rate
cum_index=[r+1,r+2];
[cumulative_favar_irf, cumulative_favar_irf_boot]=cumulative_irf(favar_irf,favar_irf_boot,cum_index);

%favar con r ottimo con gdp, hicp, shadow rate
cum_index=[r+1,r+2];
[cumulative_favar_shadow_irf, cumulative_favar_shadow_irf_boot]=cumulative_irf(favar_shadow_irf,favar_shadow_irf_boot,cum_index);

%favar baseline
cum_index=[1,2];
baselineirf=[IRF_loaded_favar_baseline([1,97],r_baseline+1,:);favar_baseline_irf(r_baseline+1,r_baseline+1,:)];
baselineirfboot=[IRF_loaded_boot_favar_baseline([1,97],r_baseline+1,:,:);favar_baseline_irf_boot(r_baseline+1,r_baseline+1,:,:)];
[cumulative_favar_baseline_irf, cumulative_favar_baseline_irf_boot]=cumulative_irf(baselineirf,baselineirfboot,cum_index);

%favar baseline shadowrate
cum_index=[1,2];
baseline_shadowirf=[IRF_loaded_favar_baseline_shadowrate([1,97],r_baseline+1,:);favar_bs_sh_irf(r_baseline+1,r_baseline+1,:)];
baseline_shadowirfboot=[IRF_loaded_boot_favar_baseline_shadowrate([1,97],r_baseline+1,:,:);favar_bl_sh_irf_boot(r_baseline+1,r_baseline+1,:,:)];
[cumulative_favar_baseline_shadow_irf, cumulative_baseline_favar_shadow_irf_boot]=cumulative_irf(baseline_shadowirf,baseline_shadowirfboot,cum_index);


% preapring cells to compare irf models and display
%irf_s={cumulative_favar_irf([r+1:r+3],r+3,:),cumulative_var_irf(:,3,:),cumulative_favar_baseline_irf,cumulative_favar_baseline_shadow_irf,cumulative_favar_shadow_irf([r+1:r+3],r+3,:)};
%irf_boot_s={cumulative_favar_irf_boot([r+1:r+3],r+3,:,:),cumulative_var_irf_boot(:,3,:,:),cumulative_favar_baseline_irf_boot,cumulative_baseline_favar_shadow_irf_boot,cumulative_favar_shadow_irf_boot([r+1:r+3],r+3,:,:)};


irf_s={cumulative_var_irf(:,3,:),cumulative_favar_irf([r+1:r+3],r+3,:),cumulative_favar_baseline_irf,cumulative_favar_shadow_irf([r+1:r+3],r+3,:),cumulative_favar_baseline_shadow_irf};
irf_boot_s={cumulative_var_irf_boot(:,3,:,:),cumulative_favar_irf_boot([r+1:r+3],r+3,:,:),cumulative_favar_baseline_irf_boot,cumulative_favar_shadow_irf_boot([r+1:r+3],r+3,:,:),cumulative_baseline_favar_shadow_irf_boot};


irf_s{1}(:,1,:)=irf_s{1}(:,1,:).*sigma_Y';
irf_boot_s{1}(:,1,:,:)=irf_boot_s{1}(:,1,:,:).*sigma_Y';
irf_s{1}([1,2],1,:)=irf_s{1}([1,2],1,:)*100;
irf_boot_s{1}([1,2],1,:,:)=irf_boot_s{1}([1,2],1,:,:)*100;

irf_s{2}(:,1,:)=irf_s{2}(:,1,:).*sigma_Y';
irf_boot_s{2}(:,1,:,:)=irf_boot_s{2}(:,1,:,:).*sigma_Y';
irf_s{2}([1,2],1,:)=irf_s{2}([1,2],1,:)*100;
irf_boot_s{2}([1,2],1,:,:)=irf_boot_s{2}([1,2],1,:,:)*100;

irf_s{3}(3,1,:)=irf_s{3}(3,1,:).*sigma_Y(3)';
irf_boot_s{3}(3,1,:,:)=irf_boot_s{3}(3,1,:,:).*sigma_Y(3)';
irf_s{3}([1,2],1,:)=irf_s{3}([1,2],1,:)*100;
irf_boot_s{3}([1,2],1,:,:)=irf_boot_s{3}([1,2],1,:,:)*100;

irf_s{4}(:,1,:)=irf_s{4}(:,1,:).*sigma_Y';
irf_boot_s{4}(:,1,:,:)=irf_boot_s{4}(:,1,:,:).*sigma_Y';
irf_s{4}([1,2],1,:)=irf_s{4}([1,2],1,:)*100;
irf_boot_s{4}([1,2],1,:,:)=irf_boot_s{4}([1,2],1,:,:)*100;

irf_s{5}(3,1,:)=irf_s{5}(3,1,:).*sigma_Y(3)';
irf_boot_s{5}(3,1,:,:)=irf_boot_s{5}(3,1,:,:).*sigma_Y(3)';
irf_s{5}([1,2],1,:)=irf_s{5}([1,2],1,:)*100;
irf_boot_s{5}([1,2],1,:,:)=irf_boot_s{5}([1,2],1,:,:)*100;

%dispplay comparison
compare_irf(irf_s,irf_boot_s,var_labels,"Policy Rate shock VAR vs FAVAR (baseline) vs FAVAR(ibrido) vs FAVAR (shadow) vs FAVAR(ibrido e shadow): livelli", ...
    ["VAR","FAVAR[r+GDP,HICP,PolicyRate","Favar Baseline","Favar Shadow", "Favar baseline shadow"],0,confidence);

%% ============================================================
%  ROBUSTNESS: DIFFERENT NUMBER OF FACTORS
%% ============================================================

irf_s_r={};
irf_s_boot_r={};

for i=1:r
    [F_r]= bernanke_cleaning(i,X_std,X_slow,policy_rate);
    Z_r = [F_r Y];

    [favar_irf_r,~,favar_irf_boot_r,~,~,~] = CholeskyIdentification(Z_r,optimal_p,H,constant,n_bootstrap,[],"",0);
    [cumulative_favar_irf_r, cumulative_favar_irf_boot_r]=cumulative_irf(favar_irf_r,favar_irf_boot_r,[i+1,i+2]);
    
    irf_s_r{i}=cumulative_favar_irf_r(i+1:i+3,i+3,:).*sigma_Y';
    irf_s_boot_r{i}=cumulative_favar_irf_boot_r(i+1:i+3,i+3,:,:).*sigma_Y';
    
    irf_s_r{i}([1,2],:,:)=irf_s_r{i}([1,2],:,:)*100;
    irf_s_boot_r{i}([1,2],:,:,:)=irf_s_boot_r{i}([1,2],:,:,:)*100;
end

%displaying
compare_irf(irf_s_r,irf_s_boot_r,var_labels,"Favar comparison - various r: from 1 factor + 3 (gdp,hicp,rate) to the optimal computed ( minus the varaibile made explicit)- livelli",[1+3:r+3-1,"optimal R"],0,confidence);

%% ============================================================
%  DISPLAYING LOADED IRF MAIN FAVAR
%% ============================================================

%FFR - IRT3M 
%IP - GDP - cum - 2
%CPI -HICP - cum - 3
%3m TREASURY BILLS - NaN
%5y TREASURY BONDS - LTIRT_EACC - 4
%MONETARY BASE - M1_EACC - cum - 5
%PERSONAL CONSUMPTION - HFCE_EA - cum - 6
%AVG HOURLY EARNINGS -
%M2 - M2_EACC - cum - 7
%DURABLE CONS - CONSD_EA - cum - 8
%HOUSING STARTS - NaN
%EXCHANGE RATE YEN - ERUS_EA - cum - 9
%NONDURABLE CONS - CONSND_EA - cum - 10
%UNEMPLOYMENT - UNETOT_EA - 11
%NEW ORDERS - NaN
%COMMODITY PRICE INDEX - HICPIN_EA - cum - 12
%DIVIDENDS - SHIX_EA - cum - 13
%CAPACITY UTIL RATE - NaN 
%EMPLOYMENT - TEMP_EA - cum - 14
%CONSUMER EXPECTATIONS - CCI_EA - 15

target_var=["LTIRT_EACC","M1_EACC","HFCE_EA","M2_EACC","CONSD_EA","ERUS_EA","CONSND_EA","UNETOT_EA","HICPIN_EA","SHIX_EA","TEMP_EA","CCI_EA"];
[idx_irfs] = get_indices(varNames,target_var);

labels_disp=["IRT3M: shock","REAL GDP","HICP",target_var];
labels_disp=replace(labels_disp,"_EACC","");
labels_disp=replace(labels_disp,"_EA","");



disp_irf=[favar_irf([r+3,r+1,r+2],:,:).*sigma_Y([3,1,2])';IRF_loaded_favar(idx_irfs,:,:)];
disp_irf_boot=[favar_irf_boot([r+3,r+1,r+2],:,:,:).*sigma_Y([3,1,2])';IRF_loaded_boot_favar(idx_irfs,:,:,:)];


cum_disp_idx=[2,3,5,6,7,8,9,10,12,13,14];
labels_disp(cum_disp_idx)=labels_disp(cum_disp_idx)+": level";
[cum_disp_irf, cum_disp_irf_boot]=cumulative_irf(disp_irf,disp_irf_boot,cum_disp_idx);

cum_disp_x100_idx=cum_disp_idx;
cum_disp_irf(cum_disp_x100_idx,:,:)=cum_disp_irf(cum_disp_x100_idx,:,:)*100;
cum_disp_irf_boot(cum_disp_x100_idx,:,:,:)=cum_disp_irf_boot(cum_disp_x100_idx,:,:,:)*100;


H=size(disp_irf,3);
x = 1:H;
x=x(:);

figure()
for i=1:max(size(labels_disp))
    subplot(4,4,i),
    band = squeeze(prctile(cum_disp_irf_boot(i,r+3,:,:), confidence, 4)); 
    lower = squeeze(band(:,1));
    upper = squeeze(band(:,2));
    lower=lower(:);
    upper=upper(:);
    x_area = [x; flipud(x)];          
    y_area = [upper; flipud(lower)];
    plot(1:H,squeeze(cum_disp_irf(i,r+3,:)),'b')
    hold on
    Patch = fill(x_area, y_area,'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    yline(0,'r--','HandleVisibility','off')
    axis tight
    hold off
    
    title(labels_disp(i))

end

%% ============================================================
%  FEVD AND R²
%% ============================================================

FEVD_result = calc_favar_fevd(disp_irf);
[R2_vector] = calc_favar_r2(X_std,Z);

Variance_Decomposition=FEVD_result(:,r+3,60);

R2=[1,1,1,R2_vector(idx_irfs)];

for i=1:size(labels_disp,2)
    print_terminal(i,1)=labels_disp(i)+":  "+Variance_Decomposition(i)+"     "+R2(i);
    if i<=3
        print_terminal(i,1)=print_terminal(i,1)+"*";

    end
end

figure();
fig = uifigure;
uit=uitable(fig,"Data",[labels_disp',Variance_Decomposition,R2'],"ColumnName",["Variable","Variande Decompoistion: " + 60,"R2"]);

toc
return


%% ============================================================
%  HISTORICAL DECOMPOSITION(WIP)
%% ============================================================
% 
% [para_hd, res_hd] = VarOLS(Z, optimal_p);
% [HD_favar, baseline_favar] = historical_decomposition(Z, para_hd, res_hd, optimal_p);
% 
% shock_idx = r + 3;   % policy rate shock (ultimo in Z)
% t_hd = t(optimal_p+1:end);
% 
% % Rescale: standardized -> unità originali
% HD_gdp  = squeeze(HD_favar(r+1, shock_idx, optimal_p+1:end)) * sigma_Y(1) * 100;
% HD_hicp = squeeze(HD_favar(r+2, shock_idx, optimal_p+1:end)) * sigma_Y(2) * 100;
% HD_rate = squeeze(HD_favar(r+3, shock_idx, optimal_p+1:end)) * sigma_Y(3);
% 
% % Plot
% figure('Name','HD – Contributo shock di politica monetaria')
% titles_hd = ["GDP (p.p.)", "HICP (p.p.)", "Policy Rate (std)"];
% hd_series = {HD_gdp, HD_hicp, HD_rate};
% colors = {[0.2 0.4 0.8],[0.8 0.2 0.2],[0.2 0.7 0.3]};
% 
% for i = 1:3
%     subplot(3,1,i)
%     bar(t_hd, hd_series{i}, 'FaceColor', colors{i}, 'EdgeColor','none')
%     yline(0,'k-','LineWidth',1)
%     title(titles_hd(i))
%     axis tight
% end
% sgtitle('Historical Decomposition – Shock di Politica Monetaria')
% 
% %% ============================================================
% %  HD – STACKED BAR: TUTTI GLI SHOCK (WIP)
% %% ============================================================
% 
% shock_labels = [arrayfun(@(x) "Factor "+x, 1:r, 'UniformOutput',false), ...
%                 "GDP shock", "HICP shock", "Rate shock"];
% shock_labels = string(shock_labels);
% 
% vars_to_plot = struct();
% vars_to_plot(1).name   = "GDP";
% vars_to_plot(1).idx    = r+1;
% vars_to_plot(1).scale  = sigma_Y(1)*100;
% vars_to_plot(1).ylabel = "p.p. contribution";
% 
% vars_to_plot(2).name   = "HICP";
% vars_to_plot(2).idx    = r+2;
% vars_to_plot(2).scale  = sigma_Y(2)*100;
% vars_to_plot(2).ylabel = "p.p. contribution";
% 
% vars_to_plot(3).name   = "Policy Rate";
% vars_to_plot(3).idx    = r+3;
% vars_to_plot(3).scale  = sigma_Y(3);
% vars_to_plot(3).ylabel = "std. dev. contribution";
% 
% % Colori: fattori in grigio, ultimi 3 in colori
% colors_shocks = [repmat([0.75 0.75 0.75], r, 1);
%                  0.2 0.5 0.8;
%                  0.8 0.3 0.2;
%                  0.2 0.7 0.3];
% 
% figure('Name','HD – All Shocks Stacked','Position',[100 100 1200 900])
% 
% for v = 1:3
%     var_idx = vars_to_plot(v).idx;
%     scale   = vars_to_plot(v).scale;
% 
%     % [T_hd x n_shocks]
%     HD_all = squeeze(HD_favar(var_idx, :, optimal_p+1:end))' * scale;
% 
%     % Deviazione osservata dalla baseline
%     observed  = Z(optimal_p+1:end, var_idx) * scale;
%     baseline_v = baseline_favar(var_idx, optimal_p+1:end)' * scale;
%     deviation  = observed - baseline_v;
% 
%     % ---- Separa positivi e negativi ----
%     HD_pos = max(HD_all, 0);   % [T_hd x n_shocks]
%     HD_neg = min(HD_all, 0);
% 
%     subplot(3,1,v)
%     hold on
% 
%     % Un unico bar() per i positivi con matrice → MATLAB fa stacking corretto
%     b_pos = bar(t_hd, HD_pos, 'stacked', 'EdgeColor','none', 'BarWidth',1);
%     b_neg = bar(t_hd, HD_neg, 'stacked', 'EdgeColor','none', 'BarWidth',1);
% 
%     % Applica colori
%     for j = 1:length(b_pos)
%         b_pos(j).FaceColor = colors_shocks(j,:);
%         b_neg(j).FaceColor = colors_shocks(j,:);
%         % Nascondi b_neg dalla legenda (duplicato)
%         b_neg(j).HandleVisibility = 'off';
%     end
% 
%     % Linea osservato
%     plot(t_hd, deviation, 'k-', 'LineWidth', 1.8, 'DisplayName','Observed deviation')
%     yline(0,'k--','LineWidth',0.8,'HandleVisibility','off')
% 
%     axis tight
%     title(vars_to_plot(v).name)
%     ylabel(vars_to_plot(v).ylabel)
% 
%     if v == 1
%         legend([b_pos, plot(NaN,NaN,'k-','LineWidth',1.5)], ...
%                [shock_labels, "Observed deviation"], ...
%                'Location','northeastoutside','FontSize',7)
%     end
% 
%     hold off
% end
% 
% sgtitle('Historical Decomposition – All Structural Shocks')