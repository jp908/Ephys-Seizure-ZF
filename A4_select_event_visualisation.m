load('select events\event_expo_select_all_loc.mat')
load('select events\event_expo_select_all_time.mat')
load('select events\event_expo_select_all_spec.mat')
load('select events\select_event_id.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')
load('Saved Variables Final\all_exposure.mat')
load('Saved Variables Final\drugnames.mat')

%% 
representive_events_time = [];
representive_events_time_long = [];
representive_events_spec = [];
representive_events_loc = [];
n = [];
for comp = 1:size(fishtocomp_id_mat,1)
    fish_id = find(fishtocomp_id_mat(comp,:));
    all_events_comp_spec = [];
    all_events_comp_time = [];
    all_events_comp_loc = [];
    for fish = fish_id
        events_spec = event_expo_select_all_spec(select_event_id(fish,:)==1,:);
        events_time = event_expo_select_all_time(select_event_id(fish,:)==1,:);
        events_loc = [repmat(fish,size(events_spec,1),1),event_expo_select_all_loc(select_event_id(fish,:)==1,:)];
%         plot(events')
%         title(drugnames{comp})
        all_events_comp_spec = [all_events_comp_spec;events_spec];
        all_events_comp_time = [all_events_comp_time;events_time];
        if length(events_loc)>1
            all_events_comp_loc = [all_events_comp_loc;events_loc];
        end
    end
    
    x = select_event_id(fish_id,:)==1;
    x(:,sum(x,1)==0)=[];
    x(sum(x,2)==0,:)=[];
    fish_id_mat = x;
%     pause
    avg_event_per_fish = [];
    for idx1 = 1:sum(fish_id_mat,1)
        avg_event_per_fish = [avg_event_per_fish;mean(all_events_comp_spec(fish_id_mat(idx1,:),:),1)];    
    end
    avg_wt_comp = mean(avg_event_per_fish,1);
    n=[n;size(fish_id_mat,1)];
    dist = pdist2(all_events_comp_spec,avg_wt_comp);
 
    id = dist==min(dist);
    id = find(id);
    if length(id)>1
        id = id(randperm(length(id),1));
    end
    
    
    if comp == 7
       [~,id] = sort(dist);
       fish_id_mat = fish_id_mat(:,id);

       idx = id(1);
       f=find(fish_id_mat(:,1));
       for idx1 = 2:size(fish_id_mat,2)
           if any(f == find(fish_id_mat(:,idx1)))
           else
             f = [f,find(fish_id_mat(:,idx1))];
             idx = [idx;id(idx1)];
           end
       end
    end
 
    if comp == 7
        representive_events_time = [representive_events_time;all_events_comp_time(idx(1:3),:)];
        representive_events_spec = [representive_events_spec;all_events_comp_spec(idx(1:3),:)];
        representive_events_loc = [representive_events_loc;all_events_comp_loc(idx(1:3),:)];

        event_long = all_exposure(all_events_comp_loc(idx(1:3),1),...
            (all_events_comp_loc(idx(1:3),2)-5000):(all_events_comp_loc(idx(1:3),2)+5000));
        representive_events_time_long = [representive_events_time_long;event_long];
    else
        representive_events_time = [representive_events_time;all_events_comp_time(id,:)];
        representive_events_spec = [representive_events_spec;all_events_comp_spec(id,:)];
        representive_events_loc = [representive_events_loc;all_events_comp_loc(id,:)];

        event_long = all_exposure(all_events_comp_loc(id,1),...
            (all_events_comp_loc(id,2)-5000):(all_events_comp_loc(id,2)+5000));
        representive_events_time_long = [representive_events_time_long;event_long];
    end
    

end

%% Plot ALL
% 
% x_axis = (0:5000)/500;
% fb(fb<1) = [];
% for comp = 1:19
%     index = 1:2:38;
%     subplot(19,2,index(comp))
%     plot(x_axis,representive_events_time(comp,:))
%     subplot(19,2,index(comp)+1)
%     pcolor(x_axis,fb,reshape(representive_events_wt(comp,:),76,5001))
%     title(drugnames{comp},'FontSize',4)
%     set(gca,'XTick',[])
%     colours = [181,73,19;...
%             227,190,48;...
%             17,81,137;...
%             161,202,204]/255;
%     mycolormap = customcolormap([0,1/4,2/3,1],colours);
% 
%     colormap('jet');
% 
%     shading interp
%     ylim([3 150])
%     caxis([0,0.04])
% end
%    

%%
comp_ind = [1,2,3;4,5,6;7,8,9;10,11,12;13,14,15;16,17,18;19,20,21];
lims = [];
for comp = 1:7
    x = comp_ind(comp,:); 
    x(x==0)=[];
    lim=[];
    for treat = 1:length(x)
        timeseries = representive_events_time_long(x(treat),:);
        maxx = max(timeseries);
        minn = min(timeseries);
        lim = [lim;[maxx,minn]];
    end
    lims = [lims;max(lim(:,1)),min(lim(:,2))];
    
end
lims = [ceil(lims(:,1)*100)/100,floor(lims(:,2)*100)/100];



%% Plot Compounds indvidually

drugnames = {'Aminophylline','Chlorpromazine','Control','Donepezil',...
    'Picrotoxin','RST','SB205607'};
x_axis = (0:10000)/2000;
% x_axis = (0:1000)/1000;
c=0; 
figure('Renderer', 'painters', 'Position', [10 10 675 700])
ha = tight_subplot(14,3,[.005 .005],[.06 .001],[.06 .01]);
for comp = [1,2,4:7,3]
    x = comp_ind(comp,:); 
    x(x==0)=[];
    
    for treat = 1:length(x)
        c=c+1;
        timeseries = representive_events_time_long(x(treat),:);
        
%         plots timeseries
        axes(ha(c))
        plot(x_axis,timeseries)
        set(gca,'XTick',[])
        ylim([lims(comp,2),lims(comp,1)])
        
        if treat>1
            set(gca,'YTick',[])
        else
            ylabel('mV')
        end
        %plots wavelet transform
        [wt,fb] = cwt(timeseries,'amor',2000);
        wt = abs(wt);
        wt(fb<1,:) = [];
        fb(fb<1,:) = [];
        
        axes(ha(c+3))
        pcolor(x_axis,fb,wt) 
        
        xticks([0,2,4])
        if comp ~= 3
            set(gca,'XTick',[])
        else
            xlabel('Time (Secs)')
        end       
        colours = [181,73,19;...
                227,190,48;...
                17,81,137;...
                161,202,204]/255;
        mycolormap = customcolormap([0,1/4,2/3,1],colours);

        colormap('jet');

        shading interp
        ylim([1 500])
        set(gca, 'YScale', 'log')
        yticks([4,15,80,350])
        if treat>1
          set(gca,'YTick',[])
        else
            ylabel('Hz')
        end
        caxis([0,0.04])
    end   
    c=c+3;
end

fig_path = strcat('Final Figures\event visualisation');
print(gcf,fig_path,'-dpng','-r600');
   
