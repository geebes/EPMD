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
T_anc = reshape(repmat(T_unq',1,run_options.npopn)',1,[]);

annual_temp = ocean.ann_theta;
[tsort,i]=sort(annual_temp);

for iyr=1:10
    
    % get abundance data
    x  = cell2mat(matObj.x(iyr,1)) .* ocean.ann_abundance;

    % reshape for locations, phenotypes, lineages
    x1 = reshape(full(x),60646,run_options.npopn,run_options.nlineages);
    
    
    x2 = squeeze(sum(x1,2)); % integrate within each ancestral lineage
    x3 = squeeze(sum(x1,3)); % integrate within each phenotype
   
    binedges=[-inf T_unq(1:end-1)+diff(T_unq)./2 inf];
    
    N2=zeros(run_options.npopn);
    N3=zeros(run_options.npopn);
    for ipop = 1:run_options.npopn % for each water temperature bin
        ii=find(annual_temp>binedges(ipop) & annual_temp<=binedges(ipop+1));
        
        N2(:,ipop)=sum(x2(ii,:),1)';
        N3(:,ipop)=sum(x3(ii,:),1)';
    end
    
    
    % plot environmental temperature against ancestral phenotype
    figure(1)
    subplot(10,2,2*iyr)
    N2(N2<=0)=NaN;
    pcolor(T_unq,T_unq,log10(N2))
    if iyr==10; xlabel('Environmental temperature'); end
    ylabel('Ancestral Phenotype')
    title(['Year ' num2str(iyr)])
    shading flat
    colorbar
    
    % plot environmental temperature against phenotype
    subplot(10,2,2*iyr-1)
    N3(N3<=0)=NaN;
    pcolor(T_unq,T_unq,log10(N3))
    if iyr==10; xlabel('Environmental temperature'); end
    ylabel('Phenotype')
    title(['Year ' num2str(iyr)])
    shading flat
    colorbar
    
    drawnow
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





















