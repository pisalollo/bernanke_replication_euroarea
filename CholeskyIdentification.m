function [cirf,cirfvar,CholBoot,cholu,CovU,u] = CholeskyIdentification(y,p,H,c,MaxBoot,ind,labels,opt)

n=size(y,2); %numero di colonne, numero di variabili 

[Y, X] = VarStr(y,c,p);       % yy and XX all the sample: crea lags matrix con costanti , ordine var ecc
T=size(Y,1);%n osservazioni
Bols=inv(X'*X)*X'*Y;
B=[Bols(2:end,:)';eye(n*(p-1)) zeros(n*(p-1),n)];%B companion form
C=Bols(1,:)';%vettore dei costanti
u=Y-X*Bols;%residui del modello
CovU=cov(u);%covariance matrix
u=u';

if max(abs(eig(B)))>=1
    %error('The eigenvalues of B must be less than 1 for stability.');
    cirf=NaN(n,n,H);
    CholBoot=NaN(n,n,H,MaxBoot);
    CovU=[];
    u=[];
    cholu = NaN(T, n);
    cirfvar=[];
    return
end


% impulse response functions
for h=1:H
    irf=B^(h-1);
    cirf(:,:,h)=irf(1:n,1:n)*chol(CovU)'; 
end

cholu=u'/(chol(CovU)')';




%_____________________________________
% Bootstrapping
%_____________________________________
for i=1:MaxBoot
    
    % generate new series of T enght
    for t=1:T+1
        if t==1
            YB(:,1)=X(1,2:end)'; %Contiente tutti i regressori al tempo 1
        else
            a=randi(size(u,2),1); %estraggo numero a caso in modo uniforme
            uu=u(:,a); %prende il residuo in posizione a
            YB(:,t)=[C;zeros(n*p-n,1)]+B*YB(:,t-1)+[uu;zeros(n*p-n,1)];%z= costante + 
        end
    end

    % estimate the new VAR
    yb=[y(1:p,:);YB(1:n,2:end)']; %prendo YB = nxp e gi 
    
    [Yn Xn] = VarStr(yb,c,p);     % yy and XX all the sample, ristima di nuovo tutto
    T=size(Yn,1);
    Bolsn=inv(Xn'*Xn)*Xn'*Yn;
    CovUBoot(:,:,i)=cov(Yn-Xn*Bolsn);
    Bn=[Bolsn(2:end,:)';eye(n*(p-1)) zeros(n*(p-1),n)]; %matrice comp fomr
   
    % impulse response functions
    for h=1:H
        irfBoot=Bn^(h-1);
        % cholesky
        CholBoot(:,:,h,i)=irfBoot(1:n,1:n)*chol(CovUBoot(:,:,i))';%4 dimennioni, variabili,orizzonte e ogni ciclo
    end
end

cirfvar = cirf;
%debugirf='cumindex  serve per indicare le variabili che vogliamo cumulate serve solo per visualizarle comulate ma NON modifica le IRF'


cirfdisplay=cirf;
CholBootdisplay=CholBoot;

cirfdisplay(ind,:,:)=cumsum(cirf(ind,:,:),3);
CholBootdisplay(ind,:,:,:)=cumsum(CholBoot(ind,:,:,:),3);

if(opt)
    k=0;
    figure
    for ii=1:n
        for jj=1:n
            k=k+1;
            subplot(n,n,k),plot(1:H,squeeze(cirfdisplay(ii,jj,:)),'k',...
            1:H,squeeze(prctile(CholBootdisplay(ii,jj,:,:),[16 84],4)),':k'),axis tight

            if ii==1, title(['Shock: ',labels(jj)]); end
            if jj==1, ylabel(labels(ii)); end
        end
    end
end
