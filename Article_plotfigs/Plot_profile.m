clear
% clc
cd ~/GitHub/EPMD/TM_data
addpath(genpath('~/GitHub/EPMD'))

load('../TM_data/full_matrix.mat');
load('../GUD_forcing/GUD_X01_abundance.mat');

%%
ii=find(ocean.lon==220.5);

A   = ocean.A(ii,ii);
lat = ocean.lat(ii);
z   = ocean.z(ii);
vol = ocean.volume(ii);
con = mean(abundance(ii,:),2)./vol;

Ainv=A+speye(size(A));

T   = ocean.ann_theta(ii);

lat2 = Ainv*(lat.*con);
z2   = Ainv*(z.*con);

% calculate velocity vectors
u = (lat2-lat);
v = (z2  -z  );

% saturate at 95th percentile
% u = sign(u).*min(abs(u),prctile(abs(u),99));
% v = sign(v).*min(abs(v),prctile(abs(v),99));



clf
quiver(lat,z,u,v)
axis ij
% set(gca,'YScale','Log')
% ylim([0 200])



















