function [FEVD_X, R2X] = calc_favar_fevd(IRFfactor)
% Lambda: (n_vars x n_factors)
% IRF_factors: (n_factors x n_shocks x horizon)

%[n_vars, n_fact] = size(Lambda);
[n_vars, n_shocks, horizon] = size(IRFfactor);

% 2. Eleviamo al quadrato (Varianza parziale istantanea)
Sq_IRF = IRFfactor .^ 2;

% 3. Cumuliamo nel tempo (Varianza parziale accumulata fino ad h)
% Questa è la somma dei quadrati da t=1 a t=h
Cum_Sq_IRF = cumsum(Sq_IRF, 3);

% 4. Calcoliamo la Varianza Totale Spiegata (MSE)
% Sommiamo i contributi di TUTTI gli shock (lungo la dim 2)
Total_Variance = sum(Cum_Sq_IRF, 2);

% 5. Calcoliamo la FEVD (Percentuale)
% Dividiamo il contributo di ogni shock per il totale
% Il risultato sarà (n_vars x n_shocks x horizon)
FEVD_X = zeros(n_vars, n_shocks, horizon);

for k = 1:n_shocks
    FEVD_X(:, k, :) = Cum_Sq_IRF(:, k, :) ./ Total_Variance;
end

% Moltiplichiamo per 100 per avere la percentuale leggibile (0-100%)
%FEVD_X = FEVD_X * 100;
end