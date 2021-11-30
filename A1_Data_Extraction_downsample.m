%% Makes Structure
path = 'C:\Users\jp422\OneDrive - University of Exeter\Electrophysiology\Without Perfusion\';

directory_alldat = dir(path);
count = 0;
all_dat = [];
for idx1 = 3:length(directory_alldat)
    path_day = strcat(path,directory_alldat(idx1).name,'\');
    dir_recordings = dir(path_day);
    for idx2 = 3:length(dir_recordings)       
        file_name = dir_recordings(idx2).name;
        if contains(file_name,'abf') && ~contains(file_name,'NOISE')
            count = count + 1;
            file_path = strcat(path_day,file_name);
            all_dat(count).data = downsample(abfload(file_path),5); %DOWNSAMPLING
            all_dat(count).filename = file_name;
            x = split(file_name,'-');
            all_dat(count).filecode = x{1};
        end
    end
end

%% reorganises data

for idx1 = 1:length(all_dat)
    data = all_dat(idx1).data;
    [~,~,sweeps] = size(data);
    reorg_dat = zeros(size(data,1),sweeps);
    for idx2 = 1:sweeps
        reorg_dat(:,idx2) = data(:,1,idx2);
    end
    all_dat(idx1).data = reorg_dat;
end

%% gets information

for idx1 = 1:length(all_dat)
    file_name = all_dat(idx1).filename;
    info = split(file_name,'(');
    
    drug_conc = split(info(2),' ');
    all_dat(idx1).drugname = drug_conc{1};
    all_dat(idx1).conc = sscanf(drug_conc{2},'%f');
    
    perfusion = split(info(3),'-');
    all_dat(idx1).pt = [sscanf(perfusion{1},'%f'),...
        sscanf(perfusion{2},'%f')];
    

end


%% Sort by name
all_nam = {all_dat(:).drugname};
all_nam_ind = unique(all_nam);

all_dat_by_comp = [];
for idx1 = 1:length(all_nam_ind)
    drug_ind = strcmp(all_nam_ind(idx1),all_nam);
    drug_vec = all_dat(drug_ind);
    drug_vec(50).filename = [];
    all_dat_by_comp = [all_dat_by_comp;drug_vec];
end

alldat = all_dat_by_comp;
clearvars -except alldat
%% sort by concentration
alldat_by_conc = alldat(1,1);
[r,c] = size(alldat);
count = 0;
for idx1 = 1:r
    concs_ind = cell2mat({alldat(idx1,:).conc});
    concs = unique(concs_ind);
    for idx2 = 1:length(concs)
        count = count + 1;
        ind = concs(idx2)==concs_ind;
        alldat_by_conc(count,1:sum(ind)) = alldat(idx1,ind);
    end
end

%% Removes ultra high concentrations of rstetra from analysis
for idx1 = 1:size(alldat_by_conc,1)-1
    if contains(alldat_by_conc(idx1,1).drugname,'Rstetrazol') ...
            && alldat_by_conc(idx1,1).conc == 0.625
        alldat_by_conc(idx1,:) = [];
    end
end

%% Combining sb2057 experiments that i paused accidentally

for idx1 = 1:size(alldat_by_conc,1)
    if contains(alldat_by_conc(idx1,1).drugname,'SB205') ...
            && alldat_by_conc(idx1,1).conc == 500
        combineddat = [alldat_by_conc(idx1,4).data,alldat_by_conc(end,3).data];
        alldat_by_conc(idx1,4).data = [];
        alldat_by_conc(idx1,3).data = combineddat;
        alldat_by_conc(idx1,3).pt(2) = 310;
    end
end

clearvars -except alldat_by_conc 
alldat = alldat_by_conc;
clear alldat_by_conc
%% REMOVES NOISEY DATA
alldat(2,5).data = []; %Aminophylline experiment with a couple of hench bits of noise in the baseline

