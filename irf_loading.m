function [IRF_loaded, IRF_loaded_boot,Lambda_raw] = irf_loading(X_std, Z, sigma, irf, irfboot,idx_slow, idx_policy)

    % L'equazione è: X = Z * Lambda' + e
    % Coefficients = Regressori \ Dipendente
    % Quindi dobbiamo fare: Z \ X_std
    % Il risultato sarà (m x n), quindi dobbiamo trasporlo per avere (n x m)
    Lambda_std = (Z \ X_std)';

    % Per le variabili "Slow" (HICP, GDP, ecc.), il coefficiente (loading)
    % relativo allo strumento di Policy deve essere ZERO.
    % Questo assicura che al tempo t=0 (impatto), lo shock di policy
    % non entri nell'equazione di queste variabili.
    Lambda_std(idx_slow, idx_policy) = 0;

    % 2. CORREZIONE DIMENSIONI SIGMA
    % Assicurati che sigma sia un vettore colonna (n x 1) per il broadcasting
    if size(sigma, 1) == 1
        sigma = sigma'; 
    end
    
    % 3. RESCALING (De-standardizzazione)
    % Moltiplichiamo ogni riga di Lambda per la deviazione standard della variabile
    % Lambda_std è (n x m), sigma è (n x 1). MATLAB lo fa in automatico.
    Lambda_raw = Lambda_std .* sigma;
    
    % Recupero dimensioni
    [n_z, n_shocks, H] = size(irf); % n_z deve essere uguale a m
    n_x = size(Lambda_raw, 1);      % Numero variabili in X (114)
    
    % --- CALCOLO IRF PUNTUALE ---
    IRF_X = zeros(n_x, n_shocks, H);
    
    for h = 1:H
        % Proiezione: (n x m) * (m x m) = (n x m)
        % Nota: irf(:,:,h) contiene le risposte di Z (righe) agli shock (colonne)
        IRF_X(:, :, h) = Lambda_raw * irf(:, :, h);
    end
    
    % --- CALCOLO IRF BOOTSTRAP ---
    n_bootstrap = size(irfboot, 4);
    IRF_X_Boot = zeros(n_x, n_shocks, H, n_bootstrap);
    
    % Usa il parfor se hai il Parallel Toolbox per velocizzare, altrimenti for
    for b = 1:n_bootstrap
        for h = 1:H
            % Stessa proiezione per ogni replica
            IRF_X_Boot(:, :, h, b) = Lambda_raw * irfboot(:, :, h, b);
        end
    end
    
    % Output
    IRF_loaded = IRF_X;
    IRF_loaded_boot = IRF_X_Boot;

end



