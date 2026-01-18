% Lorenzo Pisa Project
clear all;
close all;
varNames=string(readcell("EAdata2.xlsx","Range","B1:DO1"));
varNames_baseline=varNames;
EAdata=readmatrix("EAdata2.xlsx","Range","B3:DO313");
t = datetime(2000,1,1):calmonths(1):datetime(2025,11,1);


X_raw=EAdata;
X_std = zscore(X_raw); % Standardize the raw data
X_std_baseline=X_std;

Y=X_std(:,[1,98,75]);
Y = Y - mean(Y,'omitnan');

for i=1:size(X_std,2)
    sigma(i)=std(X_raw(:,i), 'omitnan');
    sigma_baseline(i)=sigma(i);
end

sigma_Y=sigma([1,98,75]);

%1: gdp, 98:hicpOV, 75:IRT3M_EACC
indextoremove=[1,75,98];
indextoremove_baseline=[75];

varNames(indextoremove)=[];
varNames_baseline(indextoremove_baseline)=[];

X_std(:,indextoremove)=[];
X_std_baseline(:,indextoremove_baseline)=[];

sigma(indextoremove)=[];
sigma_baseline(indextoremove_baseline)=[];

fast_keywords = string(["ERUS", "REER", "IRT", "SHIX", "BCI", "CCI", "SENTIX", "TASS", "TLB","M1_","M2_","CURR"]);

is_fast = false(1, size(varNames,2));
is_fast_baseline=false(1, size(varNames_baseline,2));

for j = 1:length(fast_keywords)
    is_fast = is_fast | contains(varNames, fast_keywords(j));
    is_fast_baseline= is_fast_baseline | contains(varNames_baseline, fast_keywords(j));
end

shadowrate=cell2mat(struct2cell(load("shadowrate_spliced.mat")));
shadowrate=zscore(shadowrate);
shadowrate=shadowrate-mean(shadowrate,'omitnan');

% --- 2. Separazione ---
X_s = X_std(:, ~is_fast); % Variabili Slow
X_s_baseline=X_std_baseline(:,~is_fast_baseline);



%% chacking end finding best r (number of factors)
min_var_explained=0.60;

[r_opt]= optimal_r(X_std,min_var_explained,"Favar con Hicp, GDP, IRT3M");
[r_opt_baseline]=optimal_r(X_std_baseline,min_var_explained,"Baseline");

%faccio combaciare forzatamente il numero di fattori propositamente
%r_opt=max(r_opt_baseline,r_opt);
%r_opt_baseline=r_opt;

%% ============================================================
%  4. FACTOR EXTRACTION (Bernanke)
%% ============================================================

r = r_opt-size(Y,2);   % numero fattori latenti tolti quelli espliciti
r_baseline = r_opt_baseline-1;

policy_rate = Y(:,3);

[Fhat]= bernanke_rotation(r,X_std,X_s,policy_rate);
[Fhat_baseline]= bernanke_rotation(r_baseline,X_std_baseline,X_s_baseline,policy_rate);

[Fhat_baseline_shadowrate]= bernanke_rotation(r_baseline,X_std_baseline,X_s_baseline,shadowrate);
[Fhat_shadowrate]= bernanke_rotation(r,X_std,X_s,shadowrate);

%% ============================================================
%  5. BUILD FAVAR STATE VECTOR
%% ============================================================

Z = [Fhat Y];     % (T x (r+k)) %favar con hicp e gdp espliciti
Z_shadowrate=[Fhat_shadowrate,Y(:,[1,2]),shadowrate];%shadow hicp e gdp espliciti

Z_baseline=[Fhat_baseline, policy_rate];%baseline
Z_baseline_shadowrate=[Fhat_baseline_shadowrate, shadowrate];%shadow baseline

%% ============================================================
%  6. VAR ESTIMATION
%% ============================================================
%stima ottimo var con massimo p_max
p_max=12;
[optimal_p, ~, res_optimalp, aic_values, ~]= VarOLSbestP(Z,p_max);
[optimal_p_shadowrate, ~, ~, aic_values, ~]= VarOLSbestP(Z_shadowrate,p_max);
[optimal_p_baseline, ~, ~, aic_values, ~]= VarOLSbestP(Z_baseline,p_max);
[optimal_p_baseline_shadowrate, ~, ~, aic_values, ~]= VarOLSbestP(Z_baseline_shadowrate,p_max);
[optimal_p_var, ~, ~, aic_values, ~]= VarOLSbestP(Y,p_max);

optimal_p=7;
optimal_p_shadowrate=optimal_p;
optimal_p_baseline=optimal_p;
optimal_p_baseline_shadowrate=optimal_p;
optimal_p_var=optimal_p;
%%
autocorr_labels=[1:r,"gdp","hicp","policy rate"];
figure
for i=1:size(Z,2)
    subplot(3,3,i)
    autocorr(res_optimalp(:,i))
    title(autocorr_labels(i))
end

%% ============================================================
%  7. IDENTIFICATION (CHOLESKY) VAR & FAVAR
%% ============================================================


p=optimal_p;
H=60;
constant=1;
n_bootstrap=100;
cum_index=[];
display_on=0;

%cholesky: FAVAR con fattori r ottimi e p ottimo del favar
favar_labels=[1:r, "GDP","HICP","Policy Rate"];
[favar_irf,~,favar_irf_boot,~,~,~] = CholeskyIdentification(Z,optimal_p,H,constant,n_bootstrap,cum_index,favar_labels,display_on);

%cholesky con var normale non Favar:VAR
var_labels=["GDP","HICP","Policy Rate"];
[var_irf,~,var_irf_boot,~,~,~] = CholeskyIdentification(Y,optimal_p_var,H,constant,n_bootstrap,cum_index,var_labels,display_on);

%cholesky con favar baseline
favar_baseline_labels=[1:r_baseline,"Policy Rate"];
[favar_baseline_irf,~,favar_baseline_irf_boot,~,~,~] = CholeskyIdentification(Z_baseline,optimal_p_baseline,H,constant,n_bootstrap,cum_index,favar_baseline_labels,display_on);

%cholesky con shadowrate
favar_baseline_shadow_labels=[1:r_baseline,"Shadow Rate"];
[favar_baseline_shadow_irf,~,favar_baseline_shadow_irf_boot,~,~,~] = CholeskyIdentification(Z_baseline_shadowrate,optimal_p_shadowrate,H,constant,n_bootstrap,cum_index,favar_baseline_shadow_labels,display_on);

%cholesky con shadowrate con hicp e gdp espliciti
favar_shadow_labels=[1:r,"GDP","HICP","Shadow Rate"];
[favar_shadow_irf,~,favar_shadow_irf_boot,~,~,~] = CholeskyIdentification(Z_shadowrate,optimal_p_baseline_shadowrate,H,constant,n_bootstrap,cum_index,favar_shadow_labels,display_on);


%% Loadings of favar con gdp,hicp,policy rate

index_favar_slow = find(~is_fast);
[IRF_loaded_favar, IRF_loaded_boot_favar,~]=irf_loading(X_std,Z,sigma,favar_irf,favar_irf_boot,index_favar_slow,r+3);
[IRF_loaded_favar_shadowrate, IRF_loaded_boot_favar_shadowrate,~]=irf_loading(X_std,Z_shadowrate,sigma,favar_shadow_irf,favar_shadow_irf_boot,index_favar_slow,r+1);

index_baseline_slow = find(~is_fast_baseline);
[IRF_loaded_favar_baseline, IRF_loaded_boot_favar_baseline,~]=irf_loading(X_std_baseline,Z_baseline,sigma_baseline,favar_baseline_irf,favar_baseline_irf_boot,index_baseline_slow,r_baseline+1);
[IRF_loaded_favar_baseline_shadowrate, IRF_loaded_boot_favar_baseline_shadowrate,~]=irf_loading(X_std_baseline,Z_baseline_shadowrate,sigma_baseline,favar_baseline_shadow_irf,favar_baseline_shadow_irf_boot,index_baseline_slow,r_baseline+1);



%% ============================================================
%%  cumulating irf to compare
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
baseline_shadowirf=[IRF_loaded_favar_baseline_shadowrate([1,97],r_baseline+1,:);favar_baseline_shadow_irf(r_baseline+1,r_baseline+1,:)];
baseline_shadowirfboot=[IRF_loaded_boot_favar_baseline_shadowrate([1,97],r_baseline+1,:,:);favar_baseline_shadow_irf_boot(r_baseline+1,r_baseline+1,:,:)];
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
    ["VAR","FAVAR[r+GDP,HICP,PolicyRate","Favar Baseline","Favar Shadow", "Favar baseline shadow"],0);


%% displaying r multiple comparison - r robustness
irf_s_r={};
irf_s_boot_r={};

for i=1:r
    [F_r]= bernanke_rotation(i,X_std,X_s,policy_rate);
    Z_r = [F_r Y];
    [favar_irf_r,~,favar_irf_boot_r,~,~,~] = CholeskyIdentification(Z_r,optimal_p,H,constant,n_bootstrap,[],"",0);
    [cumulative_favar_irf_r, cumulative_favar_irf_boot_r]=cumulative_irf(favar_irf_r,favar_irf_boot_r,[i+1,i+2]);
    irf_s_r{i}=cumulative_favar_irf_r(i+1:i+3,i+3,:).*sigma_Y';
    irf_s_boot_r{i}=cumulative_favar_irf_boot_r(i+1:i+3,i+3,:,:).*sigma_Y';
    
    irf_s_r{i}([1,2],:,:)=irf_s_r{i}([1,2],:,:)*100;
    irf_s_boot_r{i}([1,2],:,:,:)=irf_s_boot_r{i}([1,2],:,:,:)*100;
end



%displaying
compare_irf(irf_s_r,irf_s_boot_r,var_labels,"Favar comparison - various r: from 1 factor + 3 (gdp,hicp,rate) to the optimal computed (-3 +3 explicited)- livelli",[1+3:r+3-1,"optimal R"],0);

%% displaying some other irf

target_var=["LTIRT_EACC","M1_EACC","HFCE_EA","M2_EACC","CONSD_EA","ERUS_EA","CONSND_EA","UNETOT_EA","BCI_EA","PPINRG_EA","SHIX_EA","TEMP_EA","CCI_EA"];
[idx_irfs] = get_indices(varNames,target_var);

%favar_irf, favar_irf_boot r factors + gdp, hicp, policy rate
%IRF_loaded_favar, IRF_loaded_boot_favar

labels_disp=["policy rate","gdp","hicp",target_var];
disp_irf=[favar_irf([r+3,r+1,r+2],:,:).*sigma_Y([3,1,2])';IRF_loaded_favar(idx_irfs,:,:)];
disp_irf_boot=[favar_irf_boot([r+3,r+1,r+2],:,:,:).*sigma_Y([3,1,2])';IRF_loaded_boot_favar(idx_irfs,:,:,:)];


cum_disp_idx=[2,3,5,6,7,8,9,10,12,13,14,15];
%cum_disp_idx=[];
[cum_disp_irf, cum_disp_irf_boot]=cumulative_irf(disp_irf,disp_irf_boot,cum_disp_idx);

cum_disp_x100_idx=[2,3,5,6,7,8,9,10,12,13,14,15];
%cum_disp_x100_idx=[];
cum_disp_irf(cum_disp_x100_idx,:,:)=cum_disp_irf(cum_disp_x100_idx,:,:)*100;
cum_disp_irf_boot(cum_disp_x100_idx,:,:,:)=cum_disp_irf_boot(cum_disp_x100_idx,:,:,:)*100;


H=size(disp_irf,3);
x = 1:H;
x=x(:);

figure()
for i=1:max(size(labels_disp))
    subplot(4,4,i),
    band = squeeze(prctile(cum_disp_irf_boot(i,r+3,:,:), [16 84], 4)); 
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
return

%% fevd

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

print_terminal
%%
figure();
fig = uifigure;
uit=uitable(fig,"Data",[labels_disp',Variance_Decomposition,R2'],"ColumnName",["Variable","Variande Decompoistion: " + 60,"R2"]);

stileCentrato = uistyle('HorizontalAlignment', 'center');

% Applica lo stile a tutte le celle della tabella
addStyle(uit, stileCentrato, 'table');


%FFR - IRT3M
%IP - GDP 
%CPI -HICP 
%3m TREASURY BILLS - NaN
%5y TREASURY BONDS - LTIRT_EACC
%MONETARY BASE - M1_EACC
%PERSONAL CONSUMPTION - HFCE_EA
%AVG HOURLY EARNINGS -
%M2 - M2_EACC
%DURABLE CONS - CONSD_EA
%HOUSING STARTS
%EXCHANGE RATE YEN - ERUS_EA
%NONDURABLE CONS - CONSND_EA
%UNEMPLOYMENT - UNETOT_EA
%NEW ORDERS - BCI_EA (Business Climate Indicator
%COMMODITY PRICE INDEX - PPINRG_EA
%DIVIDENDS - SHIX_EA
%CAPACITY UTIL RATE
%EMPLOYMENT - TEMP_EA
%CONSUMER EXPECTATIONS - CCI_EA


return
% PARAMETRI
%%
window_size = 120; % 10 anni (se dati mensili)
[T_total, ~] = size(Fhat);
n_windows = T_total - window_size + 1;

% Preallocazione per salvare, ad esempio, la risposta del PIL allo shock Tasso
% Dimensioni: (Tempo x Orizzonte_IRF)
%Store_IRF_HICP = zeros(n_windows, H);
%Store_IRF_GDP = zeros(n_windows, H);
%
%Time_Axis = zeros(n_windows, 1);

% CICLO ROLLING
for i = 1:n_windows
    i
    % 1. Definisci gli indici temporali della finestra corrente
    idx_start = i;
    idx_end = i + window_size - 1;
    
    % 2. Seleziona i dati per questa finestra
    % Nota: I fattori F sono già estratti globalmente (metodo standard Bernanke)
    F_window = Fhat(idx_start:idx_end, :);
    Y_window = Y(idx_start:idx_end, :);
    
    % 3. Stima il VAR sulla finestra
    % (Qui chiami la tua funzione che calcola il VAR e le IRF)
    % [IRF_out] = calcola_var_irf(F_window, Y_window, ...);
    rolling_Z=[F_window Y_window];
    rolling_favar_labels=[1:r, "GDP","HICP","Policy Rate"];

    try
       [rolling_favar_irf,~,rolling_favar_irf_boot,~,~,~] = CholeskyIdentification(rolling_Z,optimal_p,H,constant,2,cum_index,rolling_favar_labels,0);
    catch
        rolling_favar_irf=zeros(r+3,r+3,H);
        size(rolling_Z)
        
    end
 
    
    
    % 4. Salva la IRF di interesse (es. GDP response to Rate Shock)
    irf_gdp_rate = rolling_favar_irf(r+1, r+3, :); 
    irf_hicp_rate = rolling_favar_irf(r+2, r+3, :); 
    
    irf_gdp_rate=cumsum(irf_gdp_rate,3);
    Store_IRF_GDP(i, :) = squeeze(irf_gdp_rate);

    irf_hicp_rate=cumsum(irf_hicp_rate,3);
    Store_IRF_HICP(i, :) = squeeze(irf_hicp_rate);
    
    % Salviamo la data centrale o finale della finestra per il grafico
    Time_Axis(i) = t(idx_end); 
end

% 5. PLOT 3D (Surface Plot) - Molto scenografico
figure;

surf(1:H, Time_Axis, Store_IRF_GDP);
title('Evoluzione temporale della risposta del PIL allo shock Tasso');
xlabel('Orizzonte IRF (mesi)');
ylabel('Tempo (Fine finestra)');
zlabel('Risposta (%)');
shading interp; view(-15, 30);