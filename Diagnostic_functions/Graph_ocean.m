clear
% clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_filename = {   'neutral_stochastic_static_GUD_X01_surface_transport',...
                     'neutral_stochastic_static_GUD_X01_surface_transport_slow',...
                     'neutral_stochastic_static_GUD_X01_weighted_transport',...
                     'neutral_stochastic_static_GUD_X01_weighted_transport_slow',...
       'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport',...
       'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport_slow',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01_slow',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1_slow',...
    ...
                     'neutral_stochastic_static_GUD_X17_surface_transport',...
                     'neutral_stochastic_static_GUD_X17_surface_transport_slow',...
                     'neutral_stochastic_static_GUD_X17_weighted_transport',...
                     'neutral_stochastic_static_GUD_X17_weighted_transport_slow',...
       'nonadaptive_dispersal_stochastic_static_GUD_X17_weighted_transport',...
       'nonadaptive_dispersal_stochastic_static_GUD_X17_weighted_transport_slow',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.01',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.01_slow',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.1',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.1_slow'};
  
input_filename = input_filename{3};
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

%%

G = digraph(ocean.B,'OmitSelfLoops');

%%
figure(1)
clf

x=ocean.ann_theta;
ax=axesm ('eqaazim','frame','on','FlineWidth',0.5,'Origin',[-46.283333,-86.666665,0]);

ax.Position=ax.Position+[-0.1 -0.2 2 4].*0.01;

fld=vector2gcmfaces(x,ocean.iface,ocean.ix,ocean.iy,ocean.iz);

[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);

Xmap(Xmap==0)=NaN;

sh=surfacem(lat, lon, Xmap);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!

axis off

colorbar

set(gcf,'Color','w')

return
%%

clear x y phi* lam* c k
% Eurasian Pole of Inaccesability
phi1=deg2rad(-46.283333); % latitude
lam0=deg2rad(-86.666667); % longitude

phi = deg2rad(ocean.lat);
lam = deg2rad(ocean.lon);

c = acos(sin(phi1).*sin(phi) + cos(phi1).*cos(phi).*cos(lam-lam0));

k = c./sin(c);

x = k .* cos(phi) .* sin(lam-lam0);
y = k .* (cos(phi1).*sin(phi) - sin(phi1).*cos(phi).*cos(lam-lam0));


figure(2)
clf
scatter(x,y,10,ocean.ann_theta,'filled')
hold on
viscircles([0,0],pi.*0.9,'Color','k')
axis tight
axis off

phi = zeros(1,1000);
lam = linspace(-pi,pi,1000);
c = acos(sin(phi1).*sin(phi) + cos(phi1).*cos(phi).*cos(lam-lam0));
k = c./sin(c);
x1 = k .* cos(phi) .* sin(lam-lam0);
y1 = k .* (cos(phi1).*sin(phi) - sin(phi1).*cos(phi).*cos(lam-lam0));

hold on
plot(x1,y1,':');



%%
plot(G,'Layout','force','XStart',x,'YStart',y)



















