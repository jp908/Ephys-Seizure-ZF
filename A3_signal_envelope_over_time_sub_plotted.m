%% Avg Activity 
run('A1_Data_Extraction_downsample')
load('Saved Variables Final\activity timeseries x-axis.mat')
load('Saved Variables Final\drugnames.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')
load('Saved Variables Final\Fish_Info.mat')
%%
[r,c] = size(alldat);
sr = 2000; %Sample rate
ds_r = 10; %downsample rate

%% Gets all auc information 

count = 0;
n_num = 0;
AUC = [];
for comp = 1:r
    n = length(alldat(comp,:));
    for exp = 1:c
        if ~isempty(alldat(comp,exp).data)
            n_num = n_num + 1;
            fishfilecode = alldat(comp,exp).filecode;
            perf = alldat(comp,exp).pt;
            start_drug = perf(2)-100;
            end_drug = perf(2)+100;
            drugname = strcat(alldat(comp,exp).drugname,{' '},string(alldat(comp,exp).conc));
                      
            fish_info{n_num,1} = [char(string(fishfilecode)),...
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
    

            %% Calculating the AUC of the absolute value of the signal for 
            %the basleine and exposure period 
            signal_base = filt1d_select(1:length(filt1d_select)/2);
            signal_expo = filt1d_select(length(filt1d_select)/2+1:length(filt1d_select));
            auc_hilb_base = zeros(1,15);
            auc_hilb_expo = zeros(1,15);
            x = linspace(1,length(signal_base),16);
            for idx1 = 1:15
                base = signal_base(x(idx1):x(idx1+1));
                auc_hilb_base(idx1) = trapz(abs(hilbert(base)));

                expo = signal_expo(x(idx1):x(idx1+1));
                auc_hilb_expo(idx1) = trapz(abs(hilbert(expo)));

                
            end
            
            auc_hilb_norm = [auc_hilb_base/mean(auc_hilb_base),auc_hilb_expo/mean(auc_hilb_base)];
            AUC = [AUC;auc_hilb_norm];
         
        end
    end
end

%% Combination figure
c = [181,73,19;227,190,48;0.3,0.3,0.3;218,101,51;216, 179, 164;17, 81, 137;161, 202, 204]/255;
% c =c([7,4,3,5,1,2,6],:);
r=size(fishtocomp_id_mat,1);
avg_comp = zeros(r,30);
std_error_comp = zeros(r,30);
select = [1,2,3;4,5,6;7,0,0;8,9,10;11,12,13;14,15,16;17,18,19];
drugnames = {'Aminophylline','Chlorpromazine','Donepezil','Picrotoxin','RST','SB205607','Control'};
p=[];
figure('position',[100,100,500,300])
ha = tight_subplot(7,2,[.01 .03],[.15 .1],[.4 .03]);
c1 = 0;
for comp = [1,2,4:6,7,3]
%     figure('position',[100,100,400,300])
%     ha = tight_subplot(1,2,[.15 .03],[.15 .1],[.15 .03]);
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
    
    c1=c1+1;
    axes(ha(c1))
    p(comp) = plot(x_axis(1:15),avg_comp(comp,1:15),'color',c(comp,:),'LineWidth',4);
    hold on
    yline(1,'--')
    errorshadow(x_axis(1:15)',avg_comp(comp,1:15)',std_error_comp(comp,1:15)',c(comp,:),0.1)
    if c1==13
     xlabel('~Time (mins)')
    end
    if c1==7
     ylabel('AUC of Absolute')
    end
    ylim([0.8,1.45])
    yticks([1,1.4])
    xlim([1,7])
    xticks([1,2,3,4,5,6,7])    
    ax1 = gca;
    ax1.FontSize = 9;
    if c1~=13
       ax1.XAxis.FontSize = 0.00000000000000000001;
    end
    
    set(gca,'box','off')
    
    c1=c1+1;
    axes(ha(c1))
    p(comp) = plot(x_axis(16:30),avg_comp(comp,16:30),'color',c(comp,:),'LineWidth',4);
    hold on
    yline(1,'--')
    errorshadow(x_axis(16:30)',avg_comp(comp,16:30)',std_error_comp(comp,16:30)',c(comp,:),0.1)
    if c1==14
     xlabel('~Time (mins)')
    end
    ylim([0.8,1.45])
    yticks([1,1.4])
    xlim([17,24])
    xticks([17,18,19,20,21,22,23,24])
    ax1 = gca;
    set(ax1,'box','off')               % gca = get current axis
    ax1.YAxis.Visible = 'off';  
    ax1.FontSize = 9;
    if c1~=14
       ax1.XAxis.FontSize = 0.00000000000000000001;
    end
end

l = legend(p([1,2,4:6,7,3]),drugnames);

l.Position = [0.0190721605706343 0.533999998337183 0.271600004043579 0.396666675567627];

set(gcf, 'color', 'w');  % makes transparent
fig_path = 'Final Figures/Avg AUC of hilbert timeseries subplotted';
print(gcf,fig_path,'-dpng','-r600');    


