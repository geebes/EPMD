% Setup path to all subdirectories of EPMD
addpath(genpath('~/GitHub/EPMD'))
clear

fnames = {'neutral_stochastic_static_GUD_X01_surface_transport',...
          'neutral_stochastic_static_GUD_X01_weighted_transport',...
          'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport',...
          'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01',...
          'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1'};

%% load simulation data and metadata
filename=fnames{2};
pathname = '../Output/';

[ocean,run_options,t_occ,x,i_lastyr] = load_EPMD_output(filename,pathname);


%% load geo grid data and Tara sites
[coastlat coastlon mygrid land] = load_geo_grid('~/GitHub/EPMD/nctiles_grid/');

%% Load Tara site indices

% Load Tara station coordiantes
load Tara_sites.mat

% Surface samples only...
insurf          = find(Tara.depth=='SRF');
Tara.lat        = Tara.lat(insurf);
Tara.lon        = Tara.lon(insurf);
Tara.depth      = Tara.depth(insurf);
Tara.StationID	= Tara.StationID(insurf);

% locations of sites sampled in size fraction #2
loaddata = load('StatID_sz2.mat');
StatID_sz2 = loaddata.StatID_sz2;


[~,ia,~] = intersect(Tara.StationID,StatID_sz2);

% get coordinates
Tlat = Tara.lat(ia);
Tlon = Tara.lon(ia);
Tlon(Tlon<0) = Tlon(Tlon<0)+360;

% Find nearest points in model grid
EucDist = sqrt( (Tlat - ocean.lat').^2 + (Tlon - ocean.lon').^2 );
[~,Tara_ind] = min(EucDist,[],2);

%% extract abundance data for year = plot_yr
plot_yr   = 100;

x_freq = x{plot_yr,1};
switch run_options.seed_dist
    case 'nonadaptive_dispersal'
        % sum across phenotypes for each seed population
        tmp = reshape(full(x_freq),60646,run_options.nphen,[]);
        x_freq  = squeeze(sum(tmp,2));
end

x_freq=x_freq(:,1:numel(ocean.sample_points)); % extract seed populations only
                                               % i.e. no global resident

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

%% plot dispersal of seed population n
figure(2)
clf

n=8;
cscale='lin';
dotsz=75;

X = x_freq(:,n).*ocean.ann_abundance;

[ax] = plot_vector(X,cscale,mygrid,ocean,'mollweid');
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!

smpl_pnts = ocean.sample_points;
scatterm(ocean.lat(smpl_pnts(n)),ocean.lon(smpl_pnts(n)),dotsz,'m','filled')
scatterm(ocean.lat(smpl_pnts(n)),ocean.lon(smpl_pnts(n)),dotsz,'k')

ch=colorbar;
caxis([0 max(caxis)]);
colormap(turbo)

%% t-SNE analysis of model data at Tara sites

x=x_freq./sum(x_freq,2); % normalise so sum x is 1 at each site
S=full(x(Tara_ind,:)); % extract Tara sites
S=S./sum(S,2); % normalise so sum x is 1 at each site

opts = statset('MaxIter',1e5);
perplexity=20;
disp(['Analysing ' filename])
disp('Performing 2D t-SNE analysis')
C2D  = tsne(S,'Algorithm','exact',...
              'NumDimensions',2,...
              'Distance',@bray_curtis,...
              'Perplexity',perplexity,...
              'Options',opts);
disp('Performing 3D t-SNE analysis')
C3D  = tsne(S,'Algorithm','exact',...
              'NumDimensions',3,...
              'Distance',@bray_curtis,...
              'Perplexity',perplexity,...
              'Options',opts);
disp('----------------------------------')

% Normalise in range 0-1 (scale to max value in all dimensions)
C2D = (C2D - min(C2D(:))) ./ (max(C2D(:)) - min(C2D(:)));
C3D = (C3D - min(C3D(:))) ./ (max(C3D(:)) - min(C3D(:)));

%%
figure(3)
clf

subplot(211)
ax=axesm('MapProjection','mollweid','frame','on','MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
scatterm(ocean.lat(Tara_ind),ocean.lon(Tara_ind),dotsz,C3D,'filled')
scatterm(ocean.lat(Tara_ind),ocean.lon(Tara_ind),dotsz,'k')



% Tara ocean basin index and colour scale
load Tara_Basin_sz2.mat
load Basin8.mat

sh=subplot(223);
shpos=sh.Position;
scatter(C2D(:,1),C2D(:,2),dotsz,grp2idx(Tara_Basin),'filled')
hold on
scatter(C2D(:,1),C2D(:,2),dotsz,'k')
axis square
set(gca,'XTick',[],'YTick',[])
box on
axis([-0.05 1.05 -0.05 1.05])
ch=colorbar;
ch.Ticks=linspace(1.4375,7.5625,8);
ch.TickLabels=cellstr(unique(Tara_Basin));
ch.Location='westoutside';
sh.Colormap=bsnclr;
caxis([1 8])
sh.Position=shpos;

sh=subplot(224);
scatter(C2D(:,1),C2D(:,2),dotsz,C3D,'filled')
hold on
scatter(C2D(:,1),C2D(:,2),dotsz,'k')
axis square
set(gca,'XTick',[],'YTick',[])
box on
axis([-0.05 1.05 -0.05 1.05])









