load('select events\event_expo_select_all_loc.mat')
load('select events\event_expo_select_all_time.mat')
load('select events\event_expo_select_all_wt.mat')
load('select events\select_event_id.mat')
load('Saved Variables Final\fishtocomp_id_mat.mat')
load('Saved Variables Final\all_exposure.mat')
load('Saved Variables Final\drugnames.mat')

%% 
n_per_fish_struct = [];
n = zeros(19,1);
sem = zeros(19,1);
p = zeros(19,1);
h = zeros(19,1);
for comp = [7,1:6,8:19]
    fish_id = find(fishtocomp_id_mat(comp,:));

    x = select_event_id(fish_id,:)==1;
    x(:,sum(x,1)==0)=[];
%     x(sum(x,2)==0,:)=[];
    
    if comp==7
        con = sum(x,2);
    end
    n_per_fish = sum(x,2);
    n_per_fish_struct(comp).n = n_per_fish;
    [p(comp),h(comp),stats(comp)] = ranksum(con,n_per_fish);   
    sem(comp) = std(n_per_fish)/sqrt(length(n_per_fish));

    n(comp)=mean(sum(x,2));
    
    
end
%% correct for multiple comparisons

p = mafdr(p);

%%

c = [0.6,0.6,0.6;...
    161,202,204;161,202,204;161,202,204;
    17,81,137;17,81,137;17,81,137;
    216,179,164;216,179,164;216,179,164;    
    218,101,51;218,101,51;218,101,51;
    227,190,48;227,190,48;227,190,48;
    181,73,19;181,73,19;181,73,19]/255;

names_orig = [];
for idx1 = 1:19
    names_orig = [names_orig;drugnames{idx1}];
end

names = categorical(names_orig);
order = cellstr(names_orig);
order = order([7,19:-1:8,6:-1:1]);
names = reordercats(names,order);

figure('position',[100,100,400,300])

b = barh(names,n);
xlabel('Mean Event Number')
hold on
er = errorbar(n,names,sem,sem,'horizontal');
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca,'TickLength',[0 0])


x = xlim;
temp=0;
for comp = [7,19:-1:8,6:-1:1]
    temp=temp+1;
    b.FaceColor = 'flat';
    b.CData(comp,:) = c(comp,:);
    if p(comp) < 0.05
       text(n(comp)+sem(comp)+0.0042*x(2),temp-0.18,'*',...
           'FontSize',15)
    end
end
box off
fig_path = strcat('Final Figures/no events');
print(gcf,fig_path,'-dpng','-r600');
