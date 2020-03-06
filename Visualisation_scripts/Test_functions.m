% Setup path to all subdirectories of EPMD
addpath(genpath('~/GitHub/EPMD'))

clear

%% load simulation data and metadata
filename = 'neutral_stochastic_static_GUD_X01_surface_transport';
pathname = '../Output/';

[ocean,run_options,t_occ,x,i_lastyr] = load_EPMD_output(filename,pathname);

%% load geo grid data
[coastlat coastlon mygrid land] = load_geo_grid('~/GitHub/EPMD/nctiles_grid/');

%% extract abundance data for year = plot_yr
plot_yr   = 100;
abundance = x{plot_yr,1};

%% plot global abundance map from output vector
cscale='log';
plot_vector(abundance(:,94),cscale,mygrid,ocean,'mollweid');

%% plot global 95th percentile dispersal times
figure(1)
clf

prc=95;
t_occ(isnan(t_occ))=100; % set unconnected points to 100 years
t_immigration = prctile(t_occ,prc,2);
t_emigration  = prctile(t_occ,prc,1);
cscale='log';

% Background map of immigration times
[ax] = plot_vector(t_immigration,cscale,mygrid,ocean,'mollweid');
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!

% Overlaid scatterplot of emigration times
scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),100,log10(t_emigration),'filled')
scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),100,'k','LineWidth',0.1)

caxis(log10([2 100]))
ch=colorbar;
ch.Ticks=log10([2 5 10 20 50 100]);
ch.TickLabels={'2','5','10','20','50','100'};
colormap(flipud(turbo))