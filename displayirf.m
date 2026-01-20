function [disp_irf,disp_irf_boot]=displayirf(irf,irfboot,rows,columns,labels,shocklabels,n_figure,confidence)

disp_irf=irf(rows,columns,:);
disp_irf_boot=irfboot(rows,columns,:,:);

n_rows=size(disp_irf,1);
n_columns=size(disp_irf,2);
H=size(irf,3);

k=0;
    figure(n_figure)
    for i=1:n_rows
        for j=1:n_columns
            k=k+1;
            subplot(n_rows,n_columns,k),
            plot(1:H,squeeze(disp_irf(i,j,:)),'k',...
            1:H,squeeze(prctile(disp_irf_boot(i,j,:,:), ...
            confidence,4)),':k'),axis tight

            if i==1, title(['Shock: ',shocklabels(j)]); end
            if j==1, ylabel(labels(i)); end
        end
    end