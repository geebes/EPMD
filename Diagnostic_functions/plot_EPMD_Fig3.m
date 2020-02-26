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
% f2.Position = [209 1 930 1344];
clf
axlim=50;
fntsz=16;


for i=1:2
    tmp=t_occ{i};
    tmp(isnan(tmp))=100;
    prc=95;
    t_immigration = prctile(tmp,prc,2);
    t_emmigration = prctile(tmp,prc,1);
    
    if i==1
        h=scatter(t_emmigration,t_immigration(ocean.sample_points),55,'k','LineWidth',1);
        hold on
        h=scatter(t_emmigration,t_immigration(ocean.sample_points),50,...
                    abs(ocean.lat(ocean.sample_points)),'filled');
        colormap(redblue)
        ch=colorbar;
        set(get(ch,'Title'),'String','Absolute Latitude','Rotation',90,...
                            'Position',[60 347.5400/2 0])
        caxis([0 90])
    else
        h=scatter(t_emmigration,t_immigration(ocean.sample_points),50,'k','Marker','.');
        box on
        xlabel('Emigration time (years)','FontSize',fntsz)
        ylabel('Immigration time (years)','FontSize',fntsz)
        set(gca,'FontSize',fntsz)
        axis square
        axis([0 axlim 0 axlim])
        plot([0 axlim],[0 axlim],'k-','LineWidth',1)
    end


    hold on
end

sname=[pathname '../Figures/Figure_3.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
