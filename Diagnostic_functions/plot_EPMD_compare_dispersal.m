clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

input_filename = {'neutral_stochastic_static_GUD_X01_surface_transport',...
                  'neutral_stochastic_static_GUD_X01_weighted_transport'};

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

for i=1:2
    pathname   = '~/GitHub/EPMD/Output/';
    matObj  = matfile([pathname input_filename{i} '.mat']);

    ocean       = matObj.ocean;
    t_occupied  = matObj.t_occupied;

    i_lastyr    = matObj.yrs_saved;
    disp([num2str(i_lastyr) ' years evaluated.'])
    %% 

    K = ocean.forcing_PCapacity(:,end);
    T = ocean.forcing_temp;
    t_occ{i}=full(t_occupied(:,1:numel(ocean.sample_points)));

    disp(['ocean is ' num2str(100.*nnz(t_occ{i})/numel(t_occ{i}),'%2.0f') '% connected.'])

    t_occ{i}(~t_occ{i}(:))=NaN;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=figure(1);
clf

prc=95;
t_immigration = prctile(t_occ{1},prc,2)./prctile(t_occ{2},prc,2);
t_emmigration = prctile(t_occ{1},prc,1)./prctile(t_occ{2},prc,1);

[ax] = plot_vector(t_immigration,'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
% th=title(['Prochlorococcus']);
ch=colorbar;
colormap(flipud(redblue))
caxis(log10([1/10 10]));
ch.Ticks=log10([1/10 1/5 1/3 1/2 1 2 3 5 10]);
ch.TickLabels={'1/10','1/5','1/3','1/2','1','2','3','5','10'};

scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),25,log10(t_emmigration),'filled')
scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),25,'k','LineWidth',0.1)

sname=['../Figures/surface_vs_depth_map.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3=figure(3);
f3.Position = [75 130 560 715];
clf

D           = t_occ{1}(ocean.sample_points,:)./ t_occ{2}(ocean.sample_points,:);
D(isnan(D)) = inf;

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
imagesc(log10(D(ind,ind)))
hold on
isdv=find(diff(ibasins))';
axlim=xlim';
% ones(2,numel(isdv)).*axlim;
plot([axlim],repmat(isdv,2,1),'k-','LineW',1.5)
plot(repmat(isdv,2,1),[axlim],'k-','LineW',1.5)

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
colormap(flipud(redblue))
caxis(log10([1/10 10]));
ch.Ticks=log10([1/10 1/5 1/3 1/2 1 2 3 5 10]);
ch.TickLabels={'1/10','1/5','1/3','1/2','1','2','3','5','10'};
ch.FontSize=11;
axis square

% sname=[pathname input_filename '/connection_matrix.png'];
% set(gcf,'Color','w')
% export_fig(sname,'-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%f5=figure(5);
f5=figure(5);

clrs={'r','b'};
for i=1:2
    
    t_immigration = prctile(t_occ{i},prc,2);
    t_emmigration = prctile(t_occ{i},prc,1);
    
    scatter(t_emmigration,t_immigration(ocean.sample_points),15,clrs{i},'filled')
    hold on
end
box on
xlabel('Emigration time')
ylabel('Immigration time')
axis square
axis([0 100 0 100])
hold on
plot([0 100],[0 100],'k-')
set(gcf,'Color','w')

sname=['../Figures/imm_vs_em.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
