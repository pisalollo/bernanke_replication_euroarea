function [HD, baseline] = historical_decomposition(Z, para, res, p)
%WIP
% HISTORICAL_DECOMPOSITION (WIP)
% para: [(1 + n*p) x n] matrix from VARTopicsOLS
%       row 1       = constants
%       rows 2:n+1  = A1 (coefficients for lag 1)
%       rows n+2:2n+1 = A2, etc.
% res:  [T-p x n] reduced-form residuals

[T, n] = size(Z);
T_res = size(res, 1);   % = T - p

%% Step 1: Estrai costante e matrici di lag
const_vec = para(1, :);          % [1 x n]

A_matrices = cell(1, p);
for lag = 1:p
    rows = 1 + (lag-1)*n + 1 : 1 + lag*n;   % blocco di n righe per lag-esimo
    A_matrices{lag} = para(rows, :);          % [n x n]
end

%% Step 2: Cholesky per shock strutturali
Sigma = (res' * res) / (T_res - 1);
C = chol(Sigma, 'lower');                    % Sigma = C*C'
structural_shocks = (C \ res')';             % [T_res x n]

%% Step 3: Coefficienti MA (rappresentazione di Wold)
% Psi(:,:,h) = risposta al passo h agli shock ridotti
Psi = zeros(n, n, T_res);
Psi(:,:,1) = eye(n);

for h = 2:T_res
    for lag = 1:min(h-1, p)
        Psi(:,:,h) = Psi(:,:,h) + A_matrices{lag}' * Psi(:,:,h-lag);
        % NOTA: A_matrices{lag} è [n x n] con convenzione A'
        % VARTopicsOLS stima y = X*para, quindi A_i = para(rows,:)
        % che corrisponde a A_i' nella forma Z_t = A_1 Z_{t-1} + ...
    end
end

% MA strutturale: Theta(h) = Psi(h) * C
Theta = zeros(n, n, T_res);
for h = 1:T_res
    Theta(:,:,h) = Psi(:,:,h) * C;
end

%% Step 4: Historical Decomposition
% HD(i,j,t) = contributo dello shock j alla variabile i al tempo t
% somma sui ritardi s: Theta(i,j,s) * eps(j, t-s)

HD = zeros(n, n, T);   % primi p periodi restano zero

for t_res = 1:T_res
    t_abs = t_res + p;
    for j = 1:n
        contrib = zeros(n, 1);
        for s = 1:t_res   % s=1 corrisponde a h=0 (impatto)
            contrib = contrib + Theta(:,j,s) * structural_shocks(t_res - s + 1, j);
        end
        HD(:, j, t_abs) = contrib;
    end
end

%% Step 5: Baseline deterministica
% baseline(t) = Z(t) - sum_j HD(:,j,t)
baseline = Z' - squeeze(sum(HD, 2));   % [n x T]

end


%per futuro lorenzo: per variabili in X_std (non in Z), ricordati i loadings: HD_obs = Lambda * HD_factors