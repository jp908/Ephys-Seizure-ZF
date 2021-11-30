run('A1_Data_Extraction_downsample')
[r,c] = size(alldat);
ind_analysis = [];

%%
sr = 2000; %Sample rate

%% load 
examples = [4,1;10,8;7,7];
AUC = zeros(157,30);
norm_spec_all_all = zeros(157,16352);
event_time_expo_all=[];
event_time_base_all=[];
event_loc_expo_all=[];
event_loc_base_all=[];
all_exposure = [];
iei = [];
fishtocomp_id_mat = zeros(r,157);
count = 0;
total_num_fish = 0;
for comp = 1:r  
    n = length(alldat(comp,:));
    n_num = 0;
    base_spec_all = zeros(n,6554);
    expo_spec_all = zeros(n,6554);
    event_spec_expo = [];
    event_time_expo = [];
    event_spec_base = [];
    event_time_base = [];
    event_loc_base = [];
    event_loc_expo = [];
    for exp = 1:c  
        if ~isempty(alldat(comp,exp).data)
            n_num = n_num + 1;
            total_num_fish = total_num_fish + 1;
            fishfilecode = alldat(comp,exp).filecode;
            perf = alldat(comp,exp).pt;
            start_drug = perf(2)-100;
            end_drug = perf(2)+100;
            drugname = strcat(alldat(comp,exp).drugname,{' '},string(alldat(comp,exp).conc));
            
            
            fish_info{total_num_fish,1} = [char(string(fishfilecode)),...
                ' ',char(drugname),' [',char(string(comp)),',',...
                char(string(exp)),']'];
            
            raw_signal = alldat(comp,exp).data(:,1:perf(2)+100);
  
   
            [length_sweeps,no_sweeps] = size(raw_signal);
           
            ind_analysis(comp,exp).drugname = drugname;
            ind_analysis(comp,exp).conc = alldat(comp,exp).conc;

            %% Put sweeps side by side to create one long data set - plot raw data
            
            raw1d = flatten(raw_signal);
            
            %% band pass filter - removes irrelevant signal
     
            [bL,aL]=butter(4,500/(sr/2),'low');
            [bH,aH]=butter(4,1/(sr/2),'high');
            filt1d = filtfilt(bL,aL,raw1d);
            filt1d = filtfilt(bH,aH,filt1d);

            filt2d = zeros(length_sweeps,no_sweeps);
            for idx1 = 1:no_sweeps
                id_start = (idx1*length_sweeps - length_sweeps)+1;
                id_end = idx1*length_sweeps;

                filt2d(:,idx1) = filt1d(1,id_start:id_end)';
            end
            
            clear raw_signal raw1d
                  
            %% selects out just drug and baseline period             

            ind = [1:100,perf(2)+1:perf(2)+100];
            filt2d_select = filt2d(:,ind);
            filt1d_select = flatten(filt2d_select);

            clear filt2d filt1d   

            %% Extracting AUC overtime
            x_axis = (1:size(filt1d_select,2))/120000;
            c2 = 0;
            bin = 0.5;
            leng_filt1d = length(filt1d_select)/(2000*60);
            auc = zeros(1,leng_filt1d/bin);
            for idx1 = bin:bin:leng_filt1d
                c2 = c2 + 1;
                ind = sum([x_axis<idx1;x_axis>(idx1-bin)],1)==2;
                signal = filt1d_select(ind);
                auc(c2) = trapz(abs(signal));
            end
    
            norm_auc = auc./mean(auc(1:15));          
            AUC(total_num_fish,:) = norm_auc;

            %% fourier transform set up

            top_freq = 500;            
            params.Fs = sr;
            params.pad = 2; % how points on x-axis of spectra gram
            params.tapers = [2 3]; %smooths spectrum
            params.fpass = [0 top_freq];  %orig [0 6000] %current upper thresh is 95% of the nyqist bound 

            [SpecRaw, ~, norm_spec_all_freq] = mtspecgramc(filt1d_select,[4.5 4.5],params);
                        
            baseline_spec = mean(SpecRaw(1:100,norm_spec_all_freq>1),1);
            exposure_spec = mean(SpecRaw(101:200,norm_spec_all_freq>1),1);  % 7 minute exposure period
            norm_spec_all_freq=norm_spec_all_freq(norm_spec_all_freq>1);
            norm_spec_all_all(total_num_fish,:) = exposure_spec./baseline_spec;  

            %% selecting for baseline period
      
            base_ind = 1:(100*length_sweeps);
            expo_ind = 100*length_sweeps+1:length(filt1d_select);
            baseline = filt1d_select(1,base_ind);
            exposure = filt1d_select(1,expo_ind);
            
            %%  wavelet transform baseline    
        
            [wt,fb] = cwt(baseline,'amor',sr);
            bwt_base = abs(wt);
            bwt_base(fb<1,:) = [];
            fb(fb<1,:) = [];
            bwt_base(fb>500,:) = [];
            fb(fb>500,:) = [];
            %% binning to different frequency bands

            n = size(bwt_base,2);
            start = [1,4,8,15,30,80,150];
            stop = [4,7,13,30,80,150,500];

            no_bins = length(start);
            bwt_base_binned = zeros(no_bins,n);
            for idx1 = 1:no_bins
                ind_base = sum([fb>=start(idx1),fb<stop(idx1)],2)==2;
                bwt_base_binned(idx1,:) = mean(bwt_base(ind_base,:),1);
            end

            %%  wavelet transform exposure period    
%              fb(fb<1) = [];
            clear wt
            
            [wt,fb] = cwt(exposure,'amor',sr);
            bwt_expo = abs(wt);
            bwt_expo(fb<1,:) = [];
            fb(fb<1) = [];
            clear wt
            
            %% binning to different frequency bands

            n = size(bwt_expo,2);
            no_bins = length(start);
            bwt_expo_binned = zeros(no_bins,n);
            for idx1 = 1:no_bins
                ind_expo = sum([fb>=start(idx1),fb<stop(idx1)],2)==2;
                bwt_expo_binned(idx1,:) = mean(bwt_expo(ind_expo,:),1);
            end
            %%
%%%%%%%%%%%%%Plots distance from mean of baseline%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bwt_base_avg_binned = mean(bwt_base_binned,2);
        bwt_base_std_binned = std(bwt_base_binned');
       
        
        distances_from_baseline_expo = pdist2(bwt_base_avg_binned',bwt_expo_binned');
        distances_from_baseline_base = pdist2(bwt_base_avg_binned',bwt_base_binned');
        


        %% Find events of interest

        thresh = 2*std(distances_from_baseline_base);    
        [events_exposure,expo_loc] = select_events_wavelet(exposure,...
            distances_from_baseline_expo,...
            thresh,sr/2);
        [events_baseline,base_loc] = select_events_wavelet(baseline,...
            distances_from_baseline_base,...
            thresh,sr/2);  

        all_exposure = [all_exposure;exposure];

    %% Plotting all event detection figures
%     figure('Renderer', 'painters', 'Position', [10 10 675 350])
%     x_axis = (1:length(distances_from_baseline_base))/2000;
%     plot(x_axis,distances_from_baseline_base)
%     hold on
%     plot(x_axis,distances_from_baseline_expo)
%     yline(thresh,'--k','LineWidth',3)
%     xlabel('Time (secs)')
%     ylabel('Euclidean Distance')
%     legend({'Baseline','Exposure'})
%     box off
%     fig_path = strcat('Final Figures\euclidean distances form baseline');
%     print(gcf,fig_path,'-dpng','-r600');
%     
%     
%     figure('Renderer', 'painters', 'Position', [10 10 675 500])
%     title('Binned Wavelet Transformation')
%     colormap('jet');
%     shading interp
%     
%     subplot(2,1,1)
%     h= imagesc(x_axis,1:8,flipud(bwt_base_binned),'CDataMapping','scaled',[0,0.04]);
%     xlabel('Time (Secs)')
%     title('Baseline')
%     
%     subplot(2,1,2)  
%     h= imagesc(x_axis,1:8,flipud(bwt_expo_binned),'CDataMapping','scaled',[0,0.04]);
%     xlabel('Time (Secs)')
%     title('Exposure')
%     
%     fig_path = strcat('Final Figures\binned wavelet transform');
%     print(gcf,fig_path,'-dpng','-r600');
    
    %%

    event_time_expo = [event_time_expo;events_exposure'];
    event_time_base = [event_time_base;events_baseline'];

    event_loc_expo = [event_loc_expo;expo_loc'];
    event_loc_base = [event_loc_base;base_loc'];
    close all

            %% Makes events id matix
         
            %exposure
            
            if total_num_fish ==1
                id1e = 1;
                id2e = size(events_exposure,2);
            else 
               id1e = id2e+1;
               id2e = id2e+size(events_exposure,2);
            end
            event_id_mat_expo(total_num_fish,id1e:id2e) = 1;
           
            %baseline

            if total_num_fish ==1
                id1b = 1;
                id2b = size(events_baseline,2);
            else 
               id1b = id2b+1;
               id2b = id2b+size(events_baseline,2);
            end
            event_id_mat_base(total_num_fish,id1b:id2b) = 1;
            
        end
    end
    %% Takes events timeseriesper compound and save all into one matrix
    %exposure
    event_time_expo_all = [event_time_expo_all;event_time_expo];
    event_loc_expo_all = [event_loc_expo_all;event_loc_expo];
    %baseline
    event_time_base_all = [event_time_base_all;event_time_base];
    event_loc_base_all = [event_loc_base_all;event_loc_base];

    %% Creates identitiy matrix for all the representatives of each fish
    count = count+1;
    if count ==1
        id11 = 1;
        id22 = n_num;
    else 
       id11 = id22+1;
       id22 = id22+n_num;
    end
    fishtocomp_id_mat(comp,id11:id22) = 1;    
end


%% Make Spectra
top_freq = 500;
            
params.Fs = sr;
params.pad = 2; % how points on x-axis of spectra gram
params.tapers = [2 3]; %smooths spectrum
params.fpass = [0 top_freq];  %orig [0 6000] %current upper thresh is 95% of the nyqist bound 

event_spec_base_all = zeros(size(event_time_base_all,1),2049);
for events = 1:size(event_time_base_all,1)
    [event_spec_base_all(events,:), ~, ~] = mtspecgramc(event_time_base_all(events,:),[1,1],params);
end

event_spec_expo_all = zeros(size(event_time_expo_all,1),2049);
for events = 1:size(event_time_expo_all,1)
    [event_spec_expo_all(events,:), ~, event_spec_freq] = mtspecgramc(event_time_expo_all(events,:),[1,1],params);
end

%%
save('Saved Variables Final/Fish_Info','fish_info')
save('Saved Variables Final/AUC','AUC')
save('Saved Variables Final/all_exposure','all_exposure')
save('Saved Variables Final/iei','iei')
save('Saved Variables Final/norm_spec_all_all','norm_spec_all_all')
save('Saved Variables Final/norm_spec_all_freq','norm_spec_all_freq')
save('Saved Variables Final/fishtocomp_id_mat','fishtocomp_id_mat')
save('Saved Variables Final/event_time_expo_all','event_time_expo_all')
save('Saved Variables Final/event_time_base_all','event_time_base_all')
save('Saved Variables Final/event_loc_expo_all','event_loc_expo_all')
save('Saved Variables Final/event_loc_base_all','event_loc_base_all')
save('Saved Variables Final/event_spec_expo_all','event_spec_expo_all')
save('Saved Variables Final/event_spec_base_all','event_spec_base_all')
save('Saved Variables Final/event_spec_freq','event_spec_freq')
save('Saved Variables Final/event_id_mat_expo','event_id_mat_expo')
save('Saved Variables Final/event_id_mat_base','event_id_mat_base')

