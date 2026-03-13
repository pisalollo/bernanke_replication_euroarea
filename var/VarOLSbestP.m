function [optimal_p, para_optimalp, res_optimalp, aic_values, autoreg]= VarOLSbestP(y,p_max)

autoreg=cell(p_max,2);
%AIC(p)=ln|∑|+2(n^2p+n)/T
optimal_p=1;
n=size(y,2);
T=size(y,1);

for p=1:p_max
    
    
    [para,res]=VARTopicsOLS(y,p);
    sigma = cov(res);
    log_sigma = log(det(sigma));
    aic_values(p) = log_sigma +2*(n^2*p+n) / T;

    %-----------------------------
    %tutti le regressioni fatte
    autoreg{p,1}=para;%% se servisse averli
    autoreg{p,2}=res;%%se servisse averli
    %----------------------------

    if aic_values(p)<=aic_values(optimal_p)
        optimal_p=p;
    end
end

para_optimalp = autoreg{optimal_p, 1};
res_optimalp = autoreg{optimal_p, 2};

