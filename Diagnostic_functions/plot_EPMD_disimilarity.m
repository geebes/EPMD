clear
% clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_filename = {'neutral_deterministic_static_GUD_X01_surface_transport',...
                  'neutral_stochastic_static_GUD_X01_surface_transport',...
                  'neutral_stochastic_static_GUD_X01_weighted_transport',...
                  'neutral_stochastic_seasonal_GUD_X01_weighted_transport',...
      'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1',... 
      'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01',... % running on belafonte
       'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport',...
                  ...
                  'neutral_stochastic_static_GUD_X17_surface_transport',... % Needs more than 100 year run
                  'neutral_stochastic_static_GUD_X17_weighted_transport',... % Needs more than 100 year run
                  'neutral_stochastic_seasonal_GUD_X17_weighted_transport',...  % Needs more han 100 year run
      'selective_dispersal_stochastic_static_GUD_X17_weighted_transport'}; % Needs more han 100 year run
  
input_filename = input_filename{7};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathname   = '~/GitHub/EPMD/Output/';
matObj  = matfile([pathname input_filename '.mat']);

ocean       = matObj.ocean;
run_options = matObj.run_options;
t_occupied  = matObj.t_occupied;

i_lastyr    = matObj.yrs_saved;
disp([num2str(i_lastyr) ' years evaluated.']) 

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

ocean.ann_abundance = mean(ocean.forcing_PCapacity,2);

K = ocean.forcing_PCapacity(:,end);
T = ocean.forcing_temp;
t_occ=full(t_occupied(:,1:numel(ocean.sample_points)));
   
disp(['ocean is ' num2str(100.*nnz(t_occ)/numel(t_occ),'%2.0f') '% connected.'])

t_occ(~t_occ(:))=NaN;


if exist([pathname input_filename])~=7
    mkdir([pathname input_filename]);
end


%%

basins=ocean.basins;
basins(ocean.lat<-40)=5;

lon=ocean.lon(ocean.sample_points);
lat=ocean.lat(ocean.sample_points);


nclust = 20;
cmap=lhsdesign(nclust,3);

tdist=t_occ(ocean.sample_points,:);

tdist=(tdist+tdist')./2;

for iyr=i_lastyr

x1  = cell2mat(matObj.x(iyr,1));
% get taxonomic composition at all sites
x = x1(ocean.sample_points,1:end-1);
% row = sites, column = taxa

BC=zeros(numel(ocean.sample_points));

for i=1:numel(ocean.sample_points)
    for j=1:numel(ocean.sample_points)
        xi=x(i,:)';
        xj=x(j,:)';
        
        C =sum(min([xi xj],[],2));
        Si=sum(xi);
        Sj=sum(xj);
        
        BC(i,j) = 1 - 2.*C ./ (Si+Sj);
        if i==j
            BC(i,j) = 0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
clf
set(gcf,'color','w')
plot(tdist(:),BC(:),'k.')
axis([0 100 0 1])
title(iyr)
xlabel('Travel time')
ylabel('BC dissimilarity')
hold on
plot([iyr iyr],[0 1],'k-')
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
clear y Z
y = squareform(BC);
Z = linkage(y,'complete');
T = cluster(Z,'maxclust',nclust);

clf
ax=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
scatterm(lat,lon,50,cmap(T,:),'filled')
title(iyr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%

for iyr=i_lastyr

    x1  = cell2mat(matObj.x(iyr,1));
    lnx=log(x1);
    lnx(isinf(lnx))=0;
    shnn=-sum(x1.*lnx,2);

    figure(3)
    clf
    
    [ax] = plot_vector(shnn,'lin',mygrid,ocean);
%     geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    title(iyr)
    colorbar
    caxis([0 0.25])
    drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%









