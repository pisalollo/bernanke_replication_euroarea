function [R2_vector] = calc_favar_r2(X,Z)

[T, N] = size(X);
for i = 1:N

    % Regressione OLS veloce: X_i = a*F + b*Y + e
    b_ols = Z \ X(:,i);
    X_hat = Z * b_ols;

    % Calcolo R2 = 1 - (SSR / SST)
    % Poiché X è standardizzata (media 0), SST è la somma dei quadrati di y
    SS_tot = sum(X(:,i).^2);
    SS_res = sum((X(:,i) - X_hat).^2);

    R2_vector(i) = 1 - (SS_res / SS_tot);
end
