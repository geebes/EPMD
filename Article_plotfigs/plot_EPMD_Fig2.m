clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

input_filename = {   'neutral_stochastic_static_GUD_X01_surface_transport',...
                     'neutral_stochastic_static_GUD_X01_weighted_transport'}; 

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);
%%
for i=1:2
    pathname   = '~/GitHub/EPMD/Output/';
    matObj  = matfile([pathname input_filename{i} '.mat']);

    clear ocean t_occupied
    ocean       = matObj.ocean;
    t_occupied  = matObj.t_occupied;

    i_lastyr{i} = matObj.yrs_saved;
    disp([num2str(i_lastyr{i}) ' years evaluated.'])
 

    T = ocean.forcing_temp;
    t_occ{i}=full(t_occupied(:,1:numel(ocean.sample_points)));
    
    ctime=t_occ{i}(:);
    ctime(isnan(ctime))=inf;
    
    disp(input_filename{i})
    disp(['ocean is ' num2str(100.*nnz(t_occ{i})/numel(t_occ{i}),'%2.0f') '% connected.'])
    disp(':::::::::::::::::::::::::::')

    t_occ{i}(~t_occ{i}(:)) = NaN;
    
    K{i}  = mean(ocean.forcing_PCapacity,2);
    xx{i} = cell2mat(matObj.x(i_lastyr{i},1));
    
    run_options{i} = matObj.run_options;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2=figure(2);
f2.Position = [209 1 930 1344];
clf


for i=1:2
    tmp=t_occ{i};
    tmp(isnan(tmp))=100;
    prc=95;
    t_immigration{i} = prctile(tmp,prc,2);
    t_emmigration{i} = prctile(tmp,prc,1);
    
    subplot(3,1,i)
    [ax] = plot_vector(t_immigration{i},'log',mygrid,ocean);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    ax.Position(2)=ax.Position(2)+(i-1).*0.04;
    
    % transport vectors
    [x,y,u,v] = get_circ_vectors(ocean); % second input is coarse graining resolution
    iplot=randsample(60646,1e4);
    h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');
    
    % th=title(['Prochlorococcus']);
    ch=colorbar;
    colormap(flipud(turbo))
    caxis(log10([2 100]));
    ch.Ticks=log10([1 2 5 10 20 50 100]);
    ch.TickLabels={'1','2','5','10','20','50','100 years'};
    ch.FontSize=18;
    ch.Position(2)=ch.Position(2)+0.025;
    ch.Position(4)=ch.Position(4)-0.050;
    t=title([char(96+i) ')'],'FontSize',18,'FontWeight','normal');
    t.Position(1)=t.Position(1)-2;
    t.Position(2)=t.Position(2)-0.3;
    
    scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),100,log10(t_emmigration{i}),'filled')
    scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),100,'k','LineWidth',0.1)
    
end
sp=subplot(313);
[ax] = plot_vector(t_immigration{2}./t_immigration{1},'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
ax.Position(2)=ax.Position(2)+0.08;

% % transport vectors
% [x,y,u,v] = get_circ_vectors(ocean); % second input is coarse graining resolution
% iplot=randsample(60646,1e4);
% h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');

ch=colorbar;
colormap(sp,flipud(redblue))
caxis(log10([0.1 10]));
ch.Ticks=log10([1/10 1/5 1/2 1 2 5 10]);
ch.TickLabels={'\div10','\div5','\div2',' ','\times2','\times5','\times10'};
ch.FontSize=18;
ch.Position(2)=ch.Position(2)+0.025;
ch.Position(4)=ch.Position(4)-0.050;
t=title([char(96+3) ')'],'FontSize',18,'FontWeight','normal');
t.Position(1)=t.Position(1)-2;
t.Position(2)=t.Position(2)-0.3;

scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),100,log10(t_emmigration{2}./t_emmigration{1}),'filled')
scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),100,'k','LineWidth',0.1)

sname=[pathname '../Figures/Figure_2.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
