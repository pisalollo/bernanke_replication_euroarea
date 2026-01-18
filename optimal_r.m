function [r_opt]= optimal_r(X_std,min_var_explained,titletext)
% Inserire dopo la standardizzazione di X
[~, ~, latent] = pca(X_std);
variance_explained = cumsum(latent) ./ sum(latent);

figure;
yyaxis left
plot(latent(1:20), 'o-', 'LineWidth', 1.5); ylabel('Autovalori');
yyaxis right
bar(variance_explained(1:20), 0.5, 'FaceAlpha', 0.3); ylabel('Varianza Cumulata');
title('Screen Plot per la scelta di r: '+titletext);
xlabel('Numero di Fattori');
grid on;


% Trova r che spiega almeno il 60% della varianza
r_opt = find(variance_explained >= min_var_explained, 1);
fprintf('Numero di fattori consigliato (60%% var): %d\n', r_opt);


