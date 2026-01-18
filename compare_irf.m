function compare_irf(irf_s,irfboot_s,labels,shocklabels,legends,area_on)


n_rows=max(size(labels));
n_columns=max(size(irf_s))+1;
n_irf_s=max(size(irf_s));
colors=lines(n_irf_s);
x = 1:size(irf_s{1},3);
x=x(:);

H=size(irf_s{1},3);

k=0;
    figure
    for i=1:n_rows
        for j=1:n_columns
            k=k+1;
            subplot(n_rows,n_columns,k)

            if j==n_columns
                for z=1:n_irf_s
                    band = squeeze(prctile(irfboot_s{z}(i,1,:,:), [16 84], 4)); 
                    %size(band)
                    
                    lower = squeeze(band(:,1));
                    upper = squeeze(band(:,2));
                    lower=lower(:);
                    upper=upper(:);
                    x_area = [x; flipud(x)];          
                    y_area = [upper; flipud(lower)];
        
                    plot(1:H,squeeze(irf_s{z}(i,1,:)),'Color',colors(z,:),'LineWidth', 1.2,'DisplayName',legends(z))
                    hold on
    
                    if area_on
                        Patch = fill(x_area, y_area, colors(z,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');
                        hold on
                    end
    
                    plot(1:H,squeeze(prctile(irfboot_s{z}(i,1,:,:), [16 84],4)),'Color', colors(z,:), 'LineStyle', ':','HandleVisibility', 'off')
                    hold on
                
                end
                
                legend('show')
                axis tight
                %axis manual
                yline(0,'k--','HandleVisibility','off')
                hold off
            else
                band = squeeze(prctile(irfboot_s{j}(i,1,:,:), [16 84], 4)); 
                size(band)
                    
                lower = squeeze(band(:,1));
                upper = squeeze(band(:,2));
                lower=lower(:);
                upper=upper(:);
                x_area = [x; flipud(x)];          
                y_area = [upper; flipud(lower)];
        
                plot(1:H,squeeze(irf_s{j}(i,1,:)),'Color',colors(j,:),'LineWidth', 1.2,'DisplayName',legends(j))
                hold on
    
            
                Patch = fill(x_area, y_area, colors(j,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');
                hold on
            
    
                plot(1:H,squeeze(prctile(irfboot_s{j}(i,1,:,:), [16 84],4)),'Color', colors(j,:), 'LineStyle', ':','HandleVisibility', 'off')
                hold on
                legend('show')
                axis tight
                %axis manual
                yline(0,'k--','HandleVisibility','off')
                hold off
            end

            %if i==1, title(['Shock: ',shocklabels(1)]); end
            if j==1, ylabel(labels(i)); end
        end
    end
    sgtitle(shocklabels(1)) 



  