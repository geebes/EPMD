clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

input_filename = 'neutral_stochastic_static_GUD_X01_surface_transport';
% input_filename = 'neutral_stochastic_static_GUD_X01_weighted_transport';
% input_filename = 'neutral_stochastic_static_GUD_X17_weighted_transport';
% input_filename = 'selective_dispersal_stochastic_static_GUD_X01_weighted_transport';
% input_filename = 'selective_dispersal_stochastic_static_GUD_X17_weighted_transport';

pathname   = '~/GitHub/EPMD/Output/';
matObj  = matfile([pathname input_filename '.mat']);

if exist([pathname input_filename])~=7
    mkdir([pathname input_filename]);
end

ocean       = matObj.ocean;
run_options = matObj.run_options;
t_occupied  = matObj.t_occupied;

i_lastyr    = matObj.yrs_saved;
disp([num2str(i_lastyr) ' years evaluated.']) 

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

mean_temp = mean(ocean.forcing_temp,2);
[ax] = plot_vector(mean_temp,'lin',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
caxis([-2 36]);
th=title(['Temperature']);
colormap(parula)
ch=colorbar('Location','SouthOutside');
ch.Ticks=[-2:2:36];
drawnow

[u,v] = get_circ_vectors(ocean);
h=quiverm(lat,lon,v,u,'k');

% sname=[pathname input_filename '/abundance.png'];
% set(gcf,'Color','w')
% export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%























