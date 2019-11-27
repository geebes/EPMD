function [x,y,u,v] = get_circ_vectors(ocean)



% get grid coords and inverse of transport matrix
lon = ocean.lon;
lat = ocean.lat;
TMi = ocean.B';

% convert lat and lon to cartesian coordinates to avoid ege effects
az = deg2rad(lon); % convert to radians
el = deg2rad(lat); % convert to radians
[X,Y,Z] = sph2cart(az,el,1);

% Apply inverse of transport matrix
X2 = TMi*X;
Y2 = TMi*Y;
Z2 = TMi*Z;

% convert back to lat and lon for plotting
[az2,el2,r2] = cart2sph(X2,Y2,Z2);

lon2 = rad2deg(az2); % convert to degrees
lat2 = rad2deg(el2); % convert to degrees

% recentre on 0 meridian
lon2(lon2<=0)=lon2(lon2<=0)+360;

% calculate velocity vectors
v = lat2-lat;
u = lon2-lon;

% saturate at 99th percentile
v = sign(v).*min(abs(v),prctile(abs(v),99));
u = sign(u).*min(abs(u),prctile(abs(u),99));

y = lat;
x = lon;

end