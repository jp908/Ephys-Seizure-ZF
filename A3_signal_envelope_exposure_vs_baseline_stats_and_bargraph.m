run('A1_Data_Extraction_downsample')

%%
[r,c] = size(alldat);
sr = 2000; %Sample rate
ds_r = 10; %downsample rate

%% Gets all auc information 

count = 0;
n_num = 0;
auc_info = struct('norm', cell(1,r), 'expo', cell(1, r), 'base', cell(1, r));
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
            auc_hilb_base = trapz(abs(hilbert(signal_base)));
        
            
            signal_expo = filt1d_select(length(filt1d_select)/2+1:length(filt1d_select));
            auc_hilb_expo = trapz(abs(hilbert(signal_expo)));
            
            auc_hilb_norm = auc_hilb_expo/auc_hilb_base; 
            
            %% Saves variables into structure
            auc_info(comp).expo(end+1) = auc_hilb_expo;
            auc_info(comp).base(end+1) = auc_hilb_base;
            auc_info(comp).norm(end+1) = auc_hilb_norm;
         
        end
    end
end


%% Does Stats on AUC - paired non parametric test
drugnames = {'Aminophylline','Chlorpromazine','Control','Donepezil',...
    'Picrotoxin','RST','SB205607'};

p_auc = zeros(1,7);
avg_auc_expo = zeros(1,7);
avg_auc_base = zeros(1,7); 
avg_auc_norm = zeros(1,7);  
sem_auc_norm = zeros(1,7);  
select = [1,2,3;4,5,6;7,0,0;8,9,10;11,12,13;14,15,16;17,18,19];
for comp = 1:7
   
    sel = select(comp,:);
    sel(sel==0) = [];
    
    auc_hilb_expo = [];
    auc_hilb_base = [];
    auc_hilb_norm = [];
    for idx1 = 1:length(sel)
        auc_hilb_expo = [auc_hilb_expo,auc_info(sel(idx1)).expo];
        auc_hilb_base = [auc_hilb_base,auc_info(sel(idx1)).base];
        auc_hilb_norm = [auc_hilb_norm,auc_info(sel(idx1)).norm];
    end
  
    p_auc(comp) = signrank(auc_hilb_base,...
        auc_hilb_expo);
    avg_auc_expo(comp) = mean(auc_hilb_expo);
    avg_auc_base(comp) = mean(auc_hilb_base);
    avg_auc_norm(comp) = mean(auc_hilb_norm)*100;    
    sem_auc_norm(comp) = std(auc_hilb_norm)/sqrt(length(auc_hilb_norm))*100;
end

%% correcting for multiple comparisons
FDR = mafdr(p_auc);  %corrected for multiple comparisons using Benjamini and Hochberg (1995)

%% Bar Graph of stats
c = [125,125,125;
    161,202,204;
    17,81,137;
    216,179,164;
    218,101,51;
    227,190,48;
    181,73,19]/255;

names = categorical(cellstr(drugnames));
order = cellstr(names([3,7:-1:4,2,1]));
names = reordercats(names,order);

% figure('position',[100,100,280,300])
figure('position',[100,100,500,400])
b = barh(names,avg_auc_norm);

xlabel('% of Baseline')
title('AUC of Hilbert transform','FontWeight','Normal')
hold on
er = errorbar(avg_auc_norm,names,sem_auc_norm,sem_auc_norm,'horizontal');
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
xlim([85,150])
x = xlim;
temp=0;
for comp = [3,7:-1:4,2,1]
    temp=temp+1;
    b.FaceColor = 'flat';
    b.CData(comp,:) = c(comp,:);
    if p_auc(comp) < 0.05
%            text(avgs_gamma(comp)+sem_gamma(comp)+0.0042*y(2),comp-0.169,'*',...
%                'FontSize',15) %adds stars on significant values
       text(avg_auc_norm(comp)+sem_auc_norm(comp)+0.0025*x(2),temp-0.16,'*',...
           'FontSize',25) %
    end
end
set(gca,'FontSize',20)
set(gca,'TickLength',[0.01 0.01])
set(gca,'box','off')
set(gcf, 'color', 'w');  % makes transparent
fig_path = 'Final Figures/Avg Signal Envelop difference Bar graph stats ';
print(gcf,fig_path,'-dpng','-r600');

    

