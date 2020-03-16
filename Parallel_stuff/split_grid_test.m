clear

clc

% Setup path to all subdirectories of EPMD
addpath(genpath('~/GitHub/EPMD'))
clear

% Candidate simulation IDs
fnames = {'neutral_stochastic_static_GUD_X01_weighted_transport_NULL'};
% Choose simulation ID
filename=fnames{1};

%% load simulation data and metadata

pathname = 'Output/';

[ocean,run_options,t_occ,x,i_lastyr] = load_EPMD_output(filename,pathname);


%% load geo grid data and Tara sites
[coastlat coastlon mygrid land] = load_geo_grid('~/GitHub/EPMD/nctiles_grid/');
%%
B=ocean.B;

nvect=60646;
n=2;
indx=round(linspace(1,nvect,n+1));

i1=indx(1):indx(2);
[j1 ~] = find(B(:,i1));
j1=unique(j1);

i2=indx(2)+1:indx(3);
[j2 ~] = find(B(:,i2));
j2=unique(j2);

x=ocean.lon;

figure(1)
clf
plot_vector(x,'lin',mygrid,ocean,'mollweid');
colorbar
% caxis([-90 90])
caxis([0 360])
title(1)
drawnow

nt=1000;
x1=repmat(x,1,nt);
for t=2:nt
    x1(:,t) = B*x1(:,t-1);
    
    clf
    plot_vector(x1(:,t),'lin',mygrid,ocean,'mollweid');
    colorbar
%    caxis([-90 90])
    caxis([0 360])
    title(t)
    drawnow
    
end
%%  

nt=1000;
xx=repmat(x,nt,1);
for t=2:nt
    xx(t,i1) = xx(t-1,j1)*B(j1,i1);
    xx(t,i2) = xx(t-1,j2)*B(j2,i2);
    
    clf
    plot_vector(xx(t,:),'lin',mygrid,ocean,'mollweid');
    colorbar
    caxis([-90 90])
    title(t)
    drawnow
    
end



