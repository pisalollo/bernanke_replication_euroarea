function [Fhat]= bernanke_cleaning(r, X_std, X_slow_std, policy_rate)
% Extracts r latent factors orthogonal to the policy rate (Bernanke et al. 2005).
%
% INPUTS:
%   r            - number of factors to extract
%   X_std        - full standardized dataset (T x N)
%   X_slow_std   - slow-moving variables only, standardized (T x N_slow)
%   policy_rate  - standardized policy rate (T x 1)
%
% OUTPUT:
%   Fhat         - cleaned factors (T x r), orthogonal to the policy rate


% A. PCA on full X → r principal factors
[~, F_all] = pca(X_std);
F_all = F_all(:, 1:r);
T = size(F_all, 1);

% B. PCA on slow variables only
[~, F_slow] = pca(X_slow_std);
F_slow = F_slow(:, 1:r);


% C. Bernanke rotation: regress F_all on [1, F_slow, rate] and strip out
% the rate's contribution. beta(end,:) isolates that component.
beta = [ones(T,1), F_slow, policy_rate] \ F_all;
Fhat = F_all - policy_rate * beta(end, :);



