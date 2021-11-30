%% Avg Activity 
load('Saved Variables Final\AUC.mat')
load('Saved Variables Final\activity timeseries x-axis.mat')
load('Saved Variables Final\drugnames.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')
load('Saved Variables Final\Fish_Info.mat')

%% Combination figure
c = [181,73,19;227,190,48;0.3,0.3,0.3;218,101,51;216, 179, 164;17, 81, 137;161, 202, 204]/255;
% c =c([7,4,3,5,1,2,6],:);
r=size(fishtocomp_id_mat,1);
avg_comp = zeros(r,30);
std_error_comp = zeros(r,30);
select = [1,2,3;4,5,6;7,0,0;8,9,10;11,12,13;14,15,16;17,18,19];
drugnames = {'Aminophylline','Chlorpromazine','Donepezil','Picrotoxin','RST','SB205607','Control'};
p=[];
figure('position',[100,100,400,300])
ha = tight_subplot(1,2,[.15 .03],[.15 .1],[.15 .03]);

for comp = 1:7
    exps = [];
    for idx2 = 1:3
        if select(comp,idx2)~=0
            ind_comp = find(fishtocomp_id_mat(select(comp,idx2),:));
            exps = [exps;AUC(ind_comp,:)];
            
        end
    end

    no_exp = size(exps,1);
    avg_comp(comp,:) = mean(exps); 
    std_error_comp(comp,:) = std(exps)./sqrt(no_exp); 
    
    axes(ha(1))
    p(comp) = plot(x_axis(1:15),avg_comp(comp,1:15),'color',c(comp,:),'LineWidth',4);
    hold on
    errorshadow(x_axis(1:15)',avg_comp(comp,1:15)',std_error_comp(comp,1:15)',c(comp,:),0.1)
    xlabel('~Time (mins)')
    ylabel('AUC of Absolute')
    ylim([0.8,1.5])
    xlim([1,7])
    xticks([1,2,3,4,5,6,7])
    ax1 = gca;
    set(gca,'box','off')
    
    axes(ha(2))
    p(comp) = plot(x_axis(16:30),avg_comp(comp,16:30),'color',c(comp,:),'LineWidth',4);
    hold on
    errorshadow(x_axis(16:30)',avg_comp(comp,16:30)',std_error_comp(comp,16:30)',c(comp,:),0.1)
    xlabel('~Time (mins)')
    ylim([0.8,1.5])
    xlim([17,24])
    xticks([17,18,19,20,21,22,23,24])
    ax1 = gca;
    set(ax1,'box','off')               % gca = get current axis
    ax1.YAxis.Visible = 'off';  
    
end

axes(ha(1))
ax = gca;
ax.FontSize = 10;
axes(ha(2))
ax = gca;
ax.FontSize = 10;
% l = legend(p([1,2,4:6,7,3]),drugnames);
% l.Position = [0.167117703443208 0.539498218316117 0.339999993667006 0.396666655540467];

set(gcf, 'color', 'w');  % makes transparent
fig_path = 'Final Figures/Avg AUC timeseries';
print(gcf,fig_path,'-dpng','-r600');    


