function [cumulative_irf, cumulative_irf_boot]=cumulative_irf(irf,irfboot,index)

cumulative_irf=irf;
cumulative_irf_boot=irfboot;

cumulative_irf(index,:,:)=cumsum(irf(index,:,:),3);%somma cumulta variabili specificate in index
cumulative_irf_boot(index,:,:,:)=cumsum(irfboot(index,:,:,:),3);


