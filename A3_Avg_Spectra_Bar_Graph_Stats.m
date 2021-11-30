%% Statistics - normalise exposure to baseline - compare vs control
load('Saved Variables Final\drugnames.mat')
load('Saved Variables Final\norm_spec_all_all.mat')
load('Saved Variables Final\norm_spec_all_freq.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')

%% BIN INTO FREQUENCY BANDS
n = size(norm_spec_all_all,1);
start = [1,4,8,15,30,80,150];
stop = [4,7,13,30,80,150,500];

no_bins = length(start);
norm_spec_all_binned = zeros(n,no_bins);
for idx1 = 1:no_bins
    ind_base = sum([norm_spec_all_freq>=start(idx1);norm_spec_all_freq<stop(idx1)],1)==2;
    norm_spec_all_binned(:,idx1) = mean(norm_spec_all_all(:,ind_base),2);
end

%% Do stats
all_treatment_exposure = [];
no_bins = size(norm_spec_all_binned,2);
r = length(drugnames);
stats = [];
all_p1 = [];
for comp = 1:r
   
   drugname = drugnames{comp};
   control = norm_spec_all_binned(fishtocomp_id_mat(7,:)==1,:);
   exposure = norm_spec_all_binned(fishtocomp_id_mat(comp,:)==1,:);
   
   all_treatment_exposure(comp).dat = exposure;
   
   p = zeros(1,no_bins);
   avg_con_expo = zeros(no_bins,2);
   sem = zeros(no_bins,2);
   percent_increase = zeros(1,no_bins);
   for idx1 = 1:no_bins
     [p(idx1),h,s] = ranksum(control(:,idx1),exposure(:,idx1));
     stats = [stats,s];
     all_p1 = [all_p1,p(idx1)];
     avg_con_expo(idx1,:) = [mean(control(:,idx1)),mean(exposure(:,idx1))];
     sem(idx1,:) = [std(control(:,idx1)),std(exposure(:,idx1))]/sqrt(size(control,1));
     percent_increase(idx1) = (avg_con_expo(idx1,2)-avg_con_expo(idx1,1))/avg_con_expo(idx1,1);
   end
      
   if comp==1
       p_tab = [[{'Drug Name: '},{'Delta'},{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'Fast'},{'Ripples'}];...
       [drugname,string(p)]];
       increase_tab = [[{'Drug Name: '},{'Delta'},{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'Fast'},{'Ripples'}];...
       [drugname,string(percent_increase)]];
       avg_tab = [[{'Drug Name: '},{'Delta'},{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'Fast'},{'Ripples'}];...
       [drugname,string(avg_con_expo(:,2)')]];
       sem_tab = [[{'Drug Name: '},{'Delta'},{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'Fast'},{'Ripples'}];...
       [drugname,string(sem(:,2)')]];
   else
       p_tab = [p_tab;[drugname,string(p)]];
       increase_tab = [increase_tab;[drugname,string(percent_increase)]];
       avg_tab = [avg_tab;[drugname,string(avg_con_expo(:,2)')]];
       sem_tab = [sem_tab;[drugname,string(sem(:,2)')]];
   end
end

save('Saved Variables Final/p_tab','p_tab')
save('Saved Variables Final/increase_tab','increase_tab')
save('Saved Variables Final/avg_tab','avg_tab')


savepath = 'C:\Users\jp422\OneDrive - University of Exeter\MATLAB\Electrophysiology\stattab ephys spectra.xls';

p_tab = cell2table(cellstr(p_tab));
increase_tab = cell2table(cellstr(increase_tab));
avg_tab = cell2table(cellstr(avg_tab));
sem_tab = cell2table(cellstr(sem_tab));

writetable(p_tab,savepath,'Sheet',1);
writetable(increase_tab,savepath,'Sheet',2);
writetable(avg_tab,savepath,'Sheet',3);


%% Correcting for multiple comparisons
all_p = zeros(19,7);
x = p_tab{2:20,2:8};
for idx1 = 1:19
    for idx2 = 1:7
        all_p(idx1,idx2) = str2num(x{idx1,idx2}); 
    end
end
FDR = mafdr(reshape(all_p,1,133));  %corrected for multiple comparisons using Benjamini and Hochberg (1995)
FDR = reshape(FDR,19,7);

x = zeros(19,7);
c=0;
for idx1 = 1:19
    for idx2 = 1:7
        c=c+1;
        x(idx1,idx2) = all_p1(c); 
       
 
        zval(idx1,idx2) = stats(c).zval;
        W(idx1,idx2) = stats(c).ranksum;
    end
end


x(FDR>0.05) = 0;
zval(FDR>0.05) = 0;
W(FDR>0.05) = 0;
FDR(FDR>0.05) = 0;

%% Make big table of results
tab = [];
bandnames = {'δ','θ','α','β','γ','High γ','HFO'};
band_names = [];
drug_names = [];
p_values = [];
corrected_p_values = [];
Test_Statistic = [];
Z_value = [];
concentration = [];
for band = 1:7
    for treat = 1:19
        if x(treat,band)>0
            
           band_names = [band_names;bandnames(band)];
           
           ghj = split(drugnames{treat},' ');
            drug_names = [drug_names;ghj(1)];
            concentration = [concentration;ghj(2)];
            p_values = [p_values;x(treat,band)];
            corrected_p_values = [corrected_p_values;FDR(treat,band)];
            Test_Statistic = [Test_Statistic;W(treat,band)];
            Z_value = [Z_value;zval(treat,band)];
            
        end
    end
end

T = table(band_names,drug_names,concentration,...
    p_values,corrected_p_values,Test_Statistic,Z_value);

%% Bar Graph of stats
c = [0.6,0.6,0.6;...
    161,202,204;161,202,204;161,202,204;
    17,81,137;17,81,137;17,81,137;
    216,179,164;216,179,164;216,179,164;    
    218,101,51;218,101,51;218,101,51;
    227,190,48;227,190,48;227,190,48;
    181,73,19;181,73,19;181,73,19]/255;

for bins = 2:8
    names = categorical(cellstr(avg_tab{2:end,1}));
    order = cellstr(avg_tab{2:end,1});
    order = order([7,19:-1:8,6:-1:1]);
    names = reordercats(names,order);

    avgs_gamma = str2double(avg_tab{2:end,bins});
    sem_gamma = str2double(sem_tab{2:end,bins});
    p_vals_gamma = FDR(:,bins-1)';
    figure('position',[100,100,195,300])
    
    b = barh(names,avgs_gamma);
    hold on
    er = errorbar(avgs_gamma,names,sem_gamma,sem_gamma,'horizontal');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    
    % plotting individual points
    hold on 
    for idx1 = 1:length(names)
        dat_vec = all_treatment_exposure(idx1).dat(:,bins-1);

        for idx2 = 1:length(dat_vec)
            scatter(dat_vec(idx2),names(idx1),10,[0,0,0],'filled',...
                'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
        end
    end
% plotting individual points
    
    x = xlim;
    temp=0;
    for comp = [7,19:-1:8,6:-1:1]
        temp=temp+1;
        b.FaceColor = 'flat';
        b.CData(comp,:) = c(comp,:);
        if p_vals_gamma(comp) < 0.05
%            text(avgs_gamma(comp)+sem_gamma(comp)+0.0042*y(2),comp-0.169,'*',...
%                'FontSize',15) %adds stars on significant values
           text(avgs_gamma(comp)+sem_gamma(comp)+0.0042*x(2),temp-0.1695,'*',...
               'FontSize',15) %
        end
        if comp==7
            xline(avgs_gamma(comp),':','LineWidth',0.9)
        end
    end
%     x = xlim;
%     for comp = [7,19:-1:8,6:-1:1]
%         b.FaceColor = 'flat';
%         b.CData(comp,:) = c(comp,:);
%         if p_vals_gamma(comp) < 0.05
% %            text(avgs_gamma(comp)+sem_gamma(comp)+0.0042*y(2),comp-0.169,'*',...
% %                'FontSize',15) %adds stars on significant values
%            text(avgs_gamma(comp)+sem_gamma(comp)+0.0042*x(2),comp-0.1695,'*',...
%                'FontSize',15) %
%         end
%     end
    
    set(gca,'TickLength',[0 0])
    set(gca,'box','off')
    set(gcf, 'color', 'w');  % makes transparent
    
    
    
    if bins == 2
        xlabel('δ')
        fig_path = 'Final Figures\Bar Graph of Delta Power.png';
    elseif bins == 3
        xlabel('θ')
        fig_path = 'Final Figures\Bar Graph of Theta Power.png';
    elseif bins == 4
        xlabel('α')
        fig_path = 'Final Figures\Bar Graph of Alpha Power.png';
    elseif bins == 5
        xlabel('β')
        fig_path = 'Final Figures\Bar Graph of Beta Power.png';
    elseif bins == 6
        xlabel('γ')
        fig_path = 'Final Figures\Bar Graph of Gamma Power.png';
    elseif bins == 7
        xlabel('High γ')
        fig_path = 'Final Figures\Bar Graph of Fast Power.png';
    elseif bins == 8
        xlabel('HFO')
        xlim([0.5,1.5])
        fig_path = 'Final Figures\Bar Graph of Ripple.png';
    end
    

set(gcf, 'color', 'w');  % makes transparent   
print(gcf,fig_path,'-dpng','-r600');

    
    close all    
end