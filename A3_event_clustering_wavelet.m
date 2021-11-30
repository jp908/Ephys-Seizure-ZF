load('Saved Variables Final\event_id_mat_base.mat')
load('Saved Variables Final\event_id_mat_expo.mat')
load('Saved Variables Final\event_time_base_all.mat')
load('Saved Variables Final\event_time_expo_all.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')
load('Saved Variables Final\event_loc_expo_all.mat')
load('Saved Variables Final\event_spec_freq.mat')
load('Saved Variables Final\drugnames.mat')

  %%  wavelet transform events  
all = [event_time_base_all;event_time_expo_all];
all_wts = zeros(size(all,1),14007);
for idx1 = 1:size(all,1)
    event = all(idx1,:);
    [wt,fb] = cwt(event,'amor',2000);
    bwt = abs(wt);
    bwt(fb<1,:) = [];
    fb(fb<1,:) = [];
    bwt(fb>500,:) = [];
    fb(fb>500,:) = [];
    %% binning to different frequency bands

    n = size(bwt,2);
    start = [1,4,8,15,30,80,150];
    stop = [4,7,13,30,80,150,500];
    no_bins = length(start);
    binned = zeros(no_bins,n);
    for idx2 = 1:no_bins
        ind_base = sum([fb>=start(idx2),fb<stop(idx2)],2)==2;
        binned(idx2,:) = mean(bwt(ind_base,:),1);
    end
    wavelet_vec = reshape(binned,1,numel(binned));
    all_wts(idx1,:) = wavelet_vec;
end

event_spec_base_all = all_wts(1:size(event_time_base_all),:);
all_wts(1:size(event_time_base_all),:) = [];
event_spec_expo_all = all_wts;
all_wts = [];

%% cluster baseline vs exposure
fish_no = 0;
t=0;
event_expo_select_all_spec = [];
event_expo_select_all_loc = [];
event_expo_select_all_time = [];
select_event_id = [];
for comp = 1:length(drugnames)
    fish_id = find(fishtocomp_id_mat(comp,:));
    for fish = fish_id
        fish_no = fish_no+1;
        base_event_id = event_id_mat_base(fish,:)==1;
        expo_event_id = event_id_mat_expo(fish,:)==1;
        
        event_base_spec_fish = event_spec_base_all(base_event_id,:);
        event_expo_spec_fish = event_spec_expo_all(expo_event_id,:);
        %normalise to avg baseline
        base_avg = mean(event_base_spec_fish,1);
        event_base_spec_fish = event_base_spec_fish./base_avg;
        event_expo_spec_fish = event_expo_spec_fish./base_avg;
        
        %bin baseline to single hz frequency bands
         n = size(event_base_spec_fish,1);
        start = 1:499;
        stop = 2:500;
        no_bins = length(start);
        event_base_spec_fish_bin = zeros(n,no_bins);
        for idx1 = 1:no_bins
            ind_base = sum([event_spec_freq>=start(idx1);...
                event_spec_freq<stop(idx1)],1)==2;
            event_base_spec_fish_bin(:,idx1) = ...
                mean(event_base_spec_fish(:,ind_base),2);
        end
        %bin baseline to single hz frequency bands
        n = size(event_expo_spec_fish,1);
        start = 1:499;
        stop = 2:500;
        no_bins = length(start);
        event_expo_spec_fish_bin = zeros(n,no_bins);
        for idx1 = 1:no_bins
            ind_base = sum([event_spec_freq>=start(idx1);...
                event_spec_freq<stop(idx1)],1)==2;
            event_expo_spec_fish_bin(:,idx1) = ...
                mean(event_expo_spec_fish(:,ind_base),2);
        end
             
        %%%%%%%%%
        event_expo_loc_fish = event_loc_expo_all(expo_event_id,:);
        event_expo_time_fish = event_time_expo_all(expo_event_id,:);
        
        
        if ~isempty(event_base_spec_fish_bin) && ~isempty(event_expo_spec_fish_bin)
         
            % finds centroid
            [~,centro] = kmeans(event_base_spec_fish_bin,1);   
            % find basleine events that is furthest from the centroid
            base_range = max(pdist2(event_base_spec_fish_bin,centro));
            % finds distance of exposure events from centroid
            clust_ind = pdist2(event_expo_spec_fish_bin,centro)>base_range;
            % selects events 
            event_expo_select_fish_spec = event_expo_spec_fish_bin(clust_ind,:);
            event_expo_select_fish_time = event_expo_time_fish(clust_ind,:);
            event_expo_select_fish_loc = event_expo_loc_fish(clust_ind,:);
            
            % save all to variable
            event_expo_select_all_spec = [event_expo_select_all_spec;event_expo_select_fish_spec];
            event_expo_select_all_loc = [event_expo_select_all_loc;event_expo_select_fish_loc];
            event_expo_select_all_time = [event_expo_select_all_time;event_expo_select_fish_time];
            
           % create identity matrix for selected events by fish 
           if fish_no==1
               id1 = 1;
               id2 = size(event_expo_select_all_spec,1);
               select_event_id(fish_no,id1:id2)=1;
           else
               id1 = sum(sum(select_event_id))+1;
               id2 = size(event_expo_select_fish_spec,1)+id1-1;
               select_event_id(fish_no,id1:id2)=1;
           end
           
           
%            if fish_no > 5
%            pause
%            end
%            [fish_no,sum(clust_ind)]
        elseif isempty(event_base_spec_fish_bin) && ~isempty(event_expo_spec_fish_bin)

          
           event_expo_select_all_spec = [event_expo_select_all_spec;event_expo_spec_fish_bin];
           event_expo_select_all_loc = [event_expo_select_all_loc;event_expo_time_fish];
           event_expo_select_all_time = [event_expo_select_all_time;event_expo_loc_fish];
           % create identity matrix for selected events by fish 
           if fish_no==1
               id1 = 1;
               id2 = size(event_expo_select_all_spec,1);
               select_event_id(fish_no,id1:id2)=1;
           else
               id1 = sum(sum(select_event_id))+1;
               id2 = size(event_expo_spec_fish,1)+id1-1;
               select_event_id(fish_no,id1:id2)=1;
           end
           
        else
          
           
        end
    end
end

save('select events\select_event_id','select_event_id')
save('select events\event_expo_select_all_spec','event_expo_select_all_spec')
save('select events\event_expo_select_all_loc','event_expo_select_all_loc')
save('select events\event_expo_select_all_time','event_expo_select_all_time')





