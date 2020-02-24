clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

input_filename = {   'neutral_stochastic_static_GUD_X01_surface_transport',...
                     'neutral_stochastic_static_GUD_X01_weighted_transport'};
  
lgnd = {'Neutral, surface transport',...
        'Neutral, depth-integrated transport'};

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);
%%
for i=1:numel(input_filename)
    pathname   = '~/GitHub/EPMD/Output/';
    matObj  = matfile([pathname input_filename{i} '.mat']);

    clear t_occupied
    ocean{i}    = matObj.ocean;
    t_occupied  = matObj.t_occupied;

    i_lastyr{i} = matObj.yrs_saved;
    disp([num2str(i_lastyr{i}) ' years evaluated.'])
 

    T = ocean{i}.forcing_temp;
    t_occ{i}=full(t_occupied(:,1:numel(ocean{i}.sample_points)));
    
    ctime=t_occ{i}(:);
    ctime(isnan(ctime))=inf;
    
    disp(input_filename{i})
    disp(['ocean is ' num2str(100.*nnz(t_occ{i})/numel(t_occ{i}),'%2.0f') '% connected.'])
    disp(':::::::::::::::::::::::::::')

    t_occ{i}(~t_occ{i}(:)) = NaN;
    
    K{i}  = mean(ocean{i}.forcing_PCapacity,2);
    xx{i} = cell2mat(matObj.x(i_lastyr{i},1));
    
    run_options{i} = matObj.run_options;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf,'Position',[123 1 1879 1344])
clf
dotsz=200;
set(0,'defaultAxesFontSize',15)

for i=1:2
    t_occupied=t_occ{i};
    t_occupied(isnan(t_occupied))=100;
    prc=95;
    t_immigration{i} = prctile(t_occupied,prc,2);
    t_emmigration{i} = prctile(t_occupied,prc,1);
    
    subplot(3,1,i)
    [ax] = plot_vector(t_immigration{i},'log',mygrid,ocean{i});
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    caxis([0 25]);
    
    scatterm(ocean{i}.lat(ocean{i}.sample_points),ocean{i}.lon(ocean{i}.sample_points),dotsz,log10(t_emmigration{i}),'filled')
    scatterm(ocean{i}.lat(ocean{i}.sample_points),ocean{i}.lon(ocean{i}.sample_points),dotsz,[0 0 0]+0.25,'LineWidth',0.1)

    % transport vectors
    [x,y,u,v] = get_circ_vectors(ocean{i}); % second input is coarse graining resolution
    iplot=randsample(60646,1e4);
    h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');

    colormap(flipud(turbo));
    caxis(log10([2 100]));
    ch=colorbar('Location','EastOutside');
    ch.Ticks=log10([1 2 5 10 20 50 100]);
    ch.TickLabels={'1','2','5','10','20','50','100 years'};
    
    text(-2.6,1.3,['(' char(96+i) ')'],'FontSize',15,'FontWeight','bold')
    drawnow
    
end
ax.Position(2)=ax.Position(2)+0.025;

hsp=subplot(3,1,3);
[ax] = plot_vector(t_immigration{2}./t_immigration{1},'log',mygrid,ocean{i});
ax.Position(2)=ax.Position(2)+0.05;
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
caxis([0 25]);

scatterm(ocean{i}.lat(ocean{i}.sample_points),ocean{i}.lon(ocean{i}.sample_points),dotsz,log10(t_emmigration{2}./t_emmigration{1}),'filled')
scatterm(ocean{i}.lat(ocean{i}.sample_points),ocean{i}.lon(ocean{i}.sample_points),dotsz,[0 0 0]+0.25,'LineWidth',0.1)

colormap(hsp,flipud(redblue));
caxis(log10([0.1 10]));
ch=colorbar('Location','EastOutside');
ch.Ticks=log10([1/10 1/5 1/2 1 2 5 10]);
ch.TickLabels={'\div10','\div5','\div2','\times1','\times2','\times5','\times10'};

text(-2.6,1.3,['(c)'],'FontSize',15,'FontWeight','bold')
drawnow
    
sname=['~/GitHub/EPMD_description/EPMD_Figures/Figure_2.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
