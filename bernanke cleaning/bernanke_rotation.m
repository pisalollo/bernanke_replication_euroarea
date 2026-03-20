function [Fhat]= bernanke_rotation(r, X_std, X_slow_std, policy_rate)
% alternative to "bernanke_cleaning.m" this has not the intercept

[~, F_all] = pca(X_std);
F_all = F_all(:,1:r);

[~, F_slow] = pca(X_slow_std);
F_slow = F_slow(:,1:r);

beta = [F_slow, policy_rate] \ F_all;
Fhat = F_all - policy_rate * beta(end,:);