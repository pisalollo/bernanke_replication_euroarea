%VAR OLS estimation of AR(p)
function [para,res]=VARTopicsOLS(data,p)
y=data(p+1:end,:);
T=size(y,1);
X=ones(T,1);
for j=1:p
    X = [X, data(p+1-j:end-j,:)];
end

%para=inv(X'*X)*X'*y;
para=X\y;

res = y - X * para;