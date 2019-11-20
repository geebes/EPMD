clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

input_filename = 'selective_dispersal_1-17_stochastic_static_GUD_X01_weighted_transport';

pathname   = '~/GitHub/EPMD/Output/';
matObj  = matfile([pathname input_filename '.mat']);

if exist([pathname input_filename])~=7
    mkdir([pathname input_filename]);
end

ocean       = matObj.ocean;
t_occupied  = matObj.t_occupied;
run_options = matObj.run_options;

i_lastyr    = matObj.yrs_saved;
disp([num2str(i_lastyr) ' years evaluated.'])
%% 

grid_load('~/GitHub/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

K = ocean.forcing_PCapacity(:,end);
T = ocean.forcing_temp;

nlineages = run_options.nlineages;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_unq = unique(run_options.T_opt);
T_anc = reshape(repmat(T_unq',1,run_options.nphen)',1,[]);

annual_temp = ocean.ann_theta;
[tsort,i]=sort(annual_temp);

iyr=i_lastyr;
    
% get abundance data
x  = cell2mat(matObj.x(iyr,1)) .* ocean.ann_abundance;

% reshape for locations, phenotypes, lineages
x1 = reshape(full(x),60646,run_options.nphen,run_options.nlineages);

x2 = squeeze(sum(x1,2)); % integrate within each ancestral lineage
x3 = squeeze(sum(x1,3)); % integrate within each phenotype

figure(99)
clf
for i = 1:nlineages
    clear x
    x=x2(:,i);
    
    subplot(3,6,i)
    [ax] = plot_vector(x,'log',mygrid,ocean);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    caxis([0 25])
    hold on
    if i<nlineages
        scatterm(ocean.lat(ocean.sample_points(i)),ocean.lon(ocean.sample_points(i)),25,'m')
    end
    drawnow
    
end




