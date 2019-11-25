clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

% input_filename = 'neutral_stochastic_static_GUD_X01_surface_transport';
% input_filename = 'neutral_stochastic_static_GUD_X01_weighted_transport';
% input_filename = 'neutral_stochastic_static_GUD_X17_weighted_transport';
input_filename = 'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_CRASHED';

pathname   = '~/GitHub/EPMD/Output/';
matObj  = matfile([pathname input_filename '.mat']);

if exist([pathname input_filename])~=7
    mkdir([pathname input_filename]);
end

ocean       = matObj.ocean;
t_occupied  = matObj.t_occupied;

i_lastyr    = matObj.yrs_saved;
disp([num2str(i_lastyr) ' years evaluated.'])
%% 

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

K = ocean.forcing_PCapacity(:,end);
T = ocean.forcing_temp;
t_occ=full(t_occupied(:,1:numel(ocean.sample_points)));
   
disp(['ocean is ' num2str(100.*nnz(t_occ)/numel(t_occ),'%2.0f') '% connected.'])

t_occ(~t_occ(:))=NaN;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=figure(1);
clf

mean_abundance = mean(ocean.forcing_PCapacity,2);
[ax] = plot_vector(mean_abundance,'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
caxis([0 25]);
% th=title(['Prochlorococcus']);
colormap(parula)
ch=colorbar('Location','SouthOutside');
ch.TickLabels={'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'};
drawnow

sname=[pathname input_filename '/abundance.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2=figure(2);
clf

prc=95;
t_immigration = prctile(t_occ,prc,2);
t_emmigration = prctile(t_occ,prc,1);

[ax] = plot_vector(t_immigration,'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
% th=title(['Prochlorococcus']);
ch=colorbar;
colormap(flipud(turbo))
caxis(log10([2 50]));
ch.Ticks=log10([1 2 5 10 20 50 100 200]);
ch.TickLabels={'1','2','5','10','20','50','100','200'};

scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),25,log10(t_emmigration),'filled')
scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),25,'k','LineWidth',0.1)

sname=[pathname input_filename '/connection_times_map.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3=figure(3);
f3.Position = [75 130 560 715];
clf

D           = t_occ(ocean.sample_points,:);
% D(isnan(D)) = inf;
indx        = find(speye(size(D)));
D(indx)     = 0; 

basins=ocean.basins;
basins(find(ocean.lat<-60))=5; % southern ocean
basins(find(ocean.basins==1 & ocean.lon<180 & ocean.lat<50  & ocean.lat>30))=6; % mediterranean

basins(basins==1)=8;
basins(basins==2)=10;
basins(basins==3)=11;
basins(basins==4)=7;
basins(basins==5)=12;
basins(basins==6)=9;
basin_names{7}='Arctic';
basin_names{8}='Atlantic';
basin_names{9}='Med. Sea';
basin_names{10}='Indian';
basin_names{11}='Pacific';
basin_names{12}='S. Ocean';


[ibasins,ind] = sort(basins(ocean.sample_points));
pcolor(log10(D(ind,ind)))
shading flat
hold on
isdv=find(diff(ibasins))';
axlim=xlim';
ones(2,numel(isdv)).*axlim;
plot([axlim],repmat(isdv,2,1),'k-','LineW',1.5)
plot(repmat(isdv,2,1),[axlim],'k-','LineW',1.5)
set(gca,'LineW',1.5,'layer','top')

for i=unique(ibasins)'
    ii=find(ibasins==i);
    text(mean(ii),axlim(2).*1.01,basin_names(i),...
         'HorizontalAlignment','l',...
         'Rotation',90)
    text(axlim(2).*1.01,mean(ii),basin_names(i),...
         'HorizontalAlignment','l')
     set(gca,'XTick',[],'YTick',[])
end
xlabel('source')
ylabel('destination')
axis xy

ch=colorbar('Location','SouthOutside');
ch.Label.String={' ';'Connectance Time'};
colormap(flipud(turbo))
caxis(log10([1/365 100]));
ch.Ticks=log10([1/365 7/365 1/12 1 10 100]);
ch.TickLabels={'day','week','month','year','decade','century'};
ch.FontSize=11;
axis square


sname=[pathname input_filename '/connection_matrix.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5=figure(5);
clf

scatter(t_emmigration,t_immigration(ocean.sample_points),15,ibasins,'filled')

box on
xlabel('Emigration time')
ylabel('Immigration time')
axis square
axis([0 100 0 100])
hold on
plot([0 100],[0 100],'k-')

sname=[pathname input_filename '/imm_vs_em.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4 = figure(4);
hold on

[N,edges] = histcounts(t_occ(:),0:(1/12):i_lastyr);

clr=rand(1,3);
plot([0 edges],[0 cumsum(N)./numel(t_occ) sum(N)./numel(t_occ)],'Color',clr,'LineW',2)
ylim([0 1])
box on
xlabel('Time (years)')
ylabel('Connectance')


sname=[pathname input_filename '/cummulative_connections.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(6)

seed_ID=1;

% get abundance data

for iyr = i_lastyr
    clf
    x  = cell2mat(matObj.x(iyr,1)) .* ocean.ann_abundance;
    
    x_i=x(:,seed_ID);
    
    [ax] = plot_vector(x_i,'log',mygrid,ocean);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    caxis([0 25])
    hold on
    title(['Year ' num2str(iyr)])
    scatterm(ocean.lat(ocean.sample_points(seed_ID)),ocean.lon(ocean.sample_points(seed_ID)),25,'m')
    colorbar
    drawnow
    
    sname=[pathname input_filename '/Seed_' num2str(seed_ID,'%03i') '_Year_'  num2str(iyr,'%03i') '.png'];
    set(gcf,'Color','w')
    export_fig(sname,'-r300')
end


return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
clf
D           = t_occ(ocean.sample_points,:);
D(isnan(D)) = inf;
indx        = find(speye(size(D)));
D(indx)     = 0; 

triuD=triu(D);

tree = linkage(triuD);

Dv=squareform(triuD);

leafOrder = optimalleaforder(tree,Dv,'Criteria','group','Transformation','inverse');

[H, nodes, outperm] = dendrogram(tree,341,'ReOrder',leafOrder);
























