function [Fhat]= bernanke_rotation(r, X_std, X_slow_std, policy_rate)

% A. PCA su tutto X
[~, F_all] = pca(X_std);
F_all = F_all(:,1:r);

% B. PCA su sole variabili lente
[~, F_slow] = pca(X_slow_std);
F_slow = F_slow(:,1:r);

% C. Rotazione Bernanke
beta = [F_slow, policy_rate] \ F_all;
Fhat = F_all - policy_rate * beta(end,:);