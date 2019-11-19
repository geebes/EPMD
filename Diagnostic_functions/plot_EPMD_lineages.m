clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

input_filename = 'lineages_stochastic_static_GUD_X01_weighted_transport';

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



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_unq = unique(run_options.T_opt);


for iyr=1:i_lastyr
    
    % get abundance data
    x  = cell2mat(matObj.x(iyr,1));

    % reshape for locations, phenotypes, lineages
    x1 = reshape(full(x),60646,run_options.npopn,run_options.nlineages);

    % Integrate lineages
    lineage_int = squeeze(sum(x1.*ocean.ann_abundance,1));
    lineage_int(lineage_int<=0) = NaN;
    
    % plot globally intgrated fraction by lineage and phenotype
    subplot(2,2,iyr)
    contourf(T_unq,T_unq,log10(lineage_int),[10:25])
    shading flat
    ylabel('Contemporary phenotype')
    xlabel('Ancestral phenotype')
    caxis([10 25])
    colorbar
    
    drawnow
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





















