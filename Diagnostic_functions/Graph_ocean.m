clear
% clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/baw103/GitHub/EPMD/TM_data/pre-rolled/surface_transport.mat')
grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);
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

% transport vectors
ocean.forcing_PCapacity = 1;
[x,y,u,v] = get_circ_vectors(ocean); % second input is coarse graining resolution
iplot=randsample(60646,1e4);
h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');

geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!

axis off

set(gcf,'Color','w')

return
%%

G = digraph(ocean.B,'OmitSelfLoops');

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



















