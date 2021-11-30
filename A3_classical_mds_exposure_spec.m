%% PCA of normalised spectra
load('Saved Variables Final\drugnames.mat')
load('Saved Variables Final\norm_spec_all_all.mat')
load('Saved Variables Final\norm_spec_all_freq.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')
r = size(drugnames,1);

c = [0.709803921568628,0.286274509803922,0.0745098039215686;
    0.709803921568628,0.286274509803922,0.0745098039215686;
    0.709803921568628,0.286274509803922,0.0745098039215686;
    0.890196078431373,0.745098039215686,0.188235294117647;
    0.890196078431373,0.745098039215686,0.188235294117647;
    0.890196078431373,0.745098039215686,0.188235294117647;
    0.00235294117647059,0.00235294117647059,0.00235294117647059;
    0.854901960784314,0.396078431372549,0.200000000000000;
    0.854901960784314,0.396078431372549,0.200000000000000;
    0.854901960784314,0.396078431372549,0.200000000000000;
    0.847058823529412,0.701960784313725,0.643137254901961;
    0.847058823529412,0.701960784313725,0.643137254901961;
    0.847058823529412,0.701960784313725,0.643137254901961;    
    0.0666666666666667,0.317647058823529,0.537254901960784;
    0.0666666666666667,0.317647058823529,0.537254901960784;
    0.0666666666666667,0.317647058823529,0.537254901960784;
    0.631372549019608,0.792156862745098,0.800000000000000;
    0.631372549019608,0.792156862745098,0.800000000000000;
    0.631372549019608,0.792156862745098,0.800000000000000];

figure('position',[100,100,300,400])
dist = pdist(norm_spec_all_all);
[Y,e]  = cmdscale(dist); 
avg_pca_comp = zeros(r,2);

for comp = [1:6,8:19,7]
    ind_comp = find(fishtocomp_id_mat(comp,:));
    avg_pca_comp(comp,:) = mean(Y(ind_comp,1:2)); 
    if comp == 7
        scatter(avg_pca_comp(comp,1),avg_pca_comp(comp,2),100,c(comp,:),'s','filled');
    elseif any(comp == [1,4,8,11,14,17])
        scatter(avg_pca_comp(comp,1),avg_pca_comp(comp,2),80,c(comp,:),'o','filled');
    elseif any(comp == [2,5,9,12,15,18])
        scatter(avg_pca_comp(comp,1),avg_pca_comp(comp,2),80,c(comp,:),'^','filled');
    elseif any(comp == [3,6,10,13,16,19])
        scatter(avg_pca_comp(comp,1),avg_pca_comp(comp,2),80,c(comp,:),'d','filled');
    end
    hold on
end

%Adds shading between points
ind = [1,3;4,6;8,10;11,13;14,16;17,19]; %orig

for idx1 = 1:6   
    shape = fill(avg_pca_comp(ind(idx1,1):ind(idx1,2),1),...
        avg_pca_comp(ind(idx1,1):ind(idx1,2),2),c(ind(idx1,1),:),...
        'linestyle','none');
    alpha(shape,0.2);
end

drugnames = drugnames([1:6,8:19,7]);
legend(drugnames,'Location','northeastoutside');
xlabel('Coord 1')
ylabel('Coord 2')
% title({'Classical Multi-dimensional Scaling of Normalised','Power Spectra Average'})

set(gcf, 'color', 'w');  % makes transparent
fig_path = 'Final Figures\Classical Multi Dimesional Scaling of Normalised Power Spectra.png';
print(gcf,fig_path,'-dpng','-r600');  


%%
%% PCA of normalised spectra
load('Saved Variables Final\drugnames.mat')
load('Saved Variables Final\norm_spec_all_all.mat')
load('Saved Variables Final\norm_spec_all_freq.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')
r = size(drugnames,1);

figure('position',[100,100,300,400])
dist = pdist(norm_spec_all_all);
[Y,e]  = cmdscale(dist); 

for comp = [1:6,8:19,7]
    ind_comp = find(fishtocomp_id_mat(comp,:));
    all_pca_comp = Y(ind_comp,1:2); 
 
    if comp == 7
        scatter(all_pca_comp(:,1),all_pca_comp(:,2),100,c(comp,:),'s','filled','MarkerFaceAlpha',0.5);
    elseif any(comp == [1,4,8,11,14,17])
        scatter(all_pca_comp(:,1),all_pca_comp(:,2),80,c(comp,:),'o','filled','MarkerFaceAlpha',0.5);
    elseif any(comp == [2,5,9,12,15,18])
        scatter(all_pca_comp(:,1),all_pca_comp(:,2),80,c(comp,:),'^','filled','MarkerFaceAlpha',0.5);
    elseif any(comp == [3,6,10,13,16,19])
        scatter(all_pca_comp(:,1),all_pca_comp(:,2),80,c(comp,:),'d','filled','MarkerFaceAlpha',0.5);
    end
    hold on
end

all_pca_all = Y(:,1:2);
select = [1,2,3;4,5,6;7,0,0;8,9,10;11,12,13;14,15,16;17,18,19];
for idx1 = 1:7  
    comp_ind = select(idx1,:);
    comp_ind(comp_ind==0) = [];
    ind = find(sum(fishtocomp_id_mat(comp_ind,:),1));
    
    selected_points = all_pca_all(ind,:);
    selected_points = selected_points(boundary(selected_points),:);
       
    shape = fill(selected_points(:,1),...
        selected_points(:,2),c(select(idx1,1),:),...
        'linestyle','none');
    alpha(shape,0.2);
end
xlabel('Coord 1')
ylabel('Coord 2')
% title({'Classical Multi-dimensional Scaling of Normalised Power','Spectra Individual Larva'})

set(gcf, 'color', 'w');  % makes transparent
fig_path = 'Final Figures\Classical Multi Dimesional Scaling of Normalised Power Spectra all fish.png';
print(gcf,fig_path,'-dpng','-r600');  