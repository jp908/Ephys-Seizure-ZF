load('Saved Variables Final/drugnames.mat')
load('Saved Variables Final/norm_spec_all_freq.mat')
load('Saved Variables Final/norm_spec_all_all.mat')
load('Saved Variables Final/fishtocomp_id_mat.mat')

no_comp = size(fishtocomp_id_mat,1);

%% Avg spectra
top_freq = 500;

avg_comp = zeros(no_comp,top_freq);
std_error_comp = zeros(no_comp,top_freq);
x_axis = 1:top_freq;
count = 0;
c2=0;
colour = [0,0,0;0.7098,0.2863,0.0745;0.8549,0.3961,0.2000;0.89,0.745,0.188];
% tiledlayout(3,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
figure('Position',  [100, 100, 500, 400])
ha = tight_subplot(3,2,[.1 .05],[.1 .1],[.1 .1]);
for comp = 1:no_comp
    
    ind_comp = find(fishtocomp_id_mat(comp,:));
    exps = norm_spec_all_all(ind_comp,:);
    
    no_exp = size(exps,1);
    bin_exp = zeros(no_exp,top_freq);
    for idx1 = 1:no_exp
        bin_exp(idx1,:) = bin_spectra(norm_spec_all_freq,exps(idx1,:),1:top_freq);
    end
    
    avg_comp(comp,:) = mean(bin_exp); 
    std_error_comp(comp,:) = std(bin_exp)./sqrt(no_exp); 
end

%% 
for comp = 1:no_comp    
    if any(comp == [1,4,8,11,14,17])
        c2=c2+1;
        count = 1;
        plots = [];
        conc = ["Control"];
        
        %PLOTS CONTROL
        axes(ha(c2))
        errorshadow(x_axis(1:top_freq)',avg_comp(7,1:top_freq)',...
            std_error_comp(7,1:top_freq)',[0.1,0.1,0.1],0.2)
         hold on
        plots(count) = plot(x_axis(1:top_freq),avg_comp(7,1:top_freq),'color',[0.1,0.1,0.1],'LineWidth',1.5);

    end
    count = count + 1;
    if comp~=7
%         count = count + 1;
        axes(ha(c2))
        errorshadow(x_axis(1:top_freq)',avg_comp(comp,1:top_freq)',...
            std_error_comp(comp,1:top_freq)',colour(count,:),0.2)
         hold on
        plots(count) = plot(x_axis(1:top_freq),avg_comp(comp,1:top_freq),'color',...
            colour(count,:),'LineWidth',1.5);
        
        set(gca,'box','off') 
        xlim([0,max(x_axis)])
       
        xticks([1,4,8,15,30,80,150,500])
        
        if comp>13
        xlabel('Frequency (Hz)')
        end
        if any(comp == [1,2,3,8,9,10,14,15,16])
        ylabel('Norm Power')
        end
        head = [char(drugnames{comp}),' Normalised Power Spectra'];
        temp = split(head,' ');   
        title(temp(1),'FontSize',10)
        ax = gca;        
        set(ax, 'XScale', 'log')
        
%         offsetAxes(ax,8,4)
        ylim auto
        y = ylim;
        ylim([0,y(2)])
        ax.XAxis.FontSize = 10;
        ax.YAxis.FontSize = 10;
        
        temp = split(drugnames{comp},' ');
        
        conc = [conc,temp(2)];
        h = legend(plots,conc);
        legend boxoff
        set(h,'FontSize',8);

    end
end

%%
% 
fig_path = strcat('Final Figures\Avg Spectra of Exposure.png');
set(gcf, 'color', 'w');  % makes transparent   
print(gcf,fig_path,'-dpng','-r600');