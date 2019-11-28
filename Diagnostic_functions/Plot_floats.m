clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;


run_options.TM_scheme       = 'surface_transport'; % 'surface_transport' or 'GUD_X01_weighted_transport', or similar
run_options.seedseed        = 2;            % seed for global seeding sites

% Initialise ocean structural array
disp('Allocate ocean metadata and generate seeding points')
[ocean] = allocate_ocean(run_options,run_options.seedseed,'quad');

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

A=ocean.B;
B=speye(size(A));
B=B(:,ocean.sample_points);
B=repmat(B,n,1);
B=reshape(B,60646,94*n);



%% calculate float dispersal
n=10;
tmax=4*365*10;

disp('Calculating float dispersal')

indx=zeros(tmax,size(B,2));
[~,I]=max(B,[],1);
indx(1,:)=I;
for t=2:tmax
    B2=(A*B).*rand(size(B));
    [~,I]=max(B2,[],1);
    
    indx(t,:)=I;
    B=sparse(I,1:numel(I),1,size(B,1),size(B,2));

end
    

%%

disp('Plotting')

clf
ax=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
colormap(flipud(gray))
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
h=plotm(lat,lon,'k','LineW',1);

tplot=4*365*10;

lat=ocean.lat(indx(1:tplot,:));
lon=ocean.lon(indx(1:tplot,:));

lat=lat+randn([1 n*94]).*0.3;
lon=lon+randn([1 n*94]).*0.3;


drawnow

clrs=lhsdesign(94,3);
clrsn=repmat(clrs(:),1,n)';
clrsn=reshape(clrsn,n*94,3);
clear rgba cd
for i=1:numel(h)
    rgba(1,:)=linspace(1,0.25,tplot).*clrsn(i,1);
    rgba(2,:)=linspace(1,0.25,tplot).*clrsn(i,2);
    rgba(3,:)=linspace(1,0.25,tplot).*clrsn(i,3);
    rgba(4,:)=linspace(0,1,tplot);
    cd = uint8(rgba*255);
    set(h(i).Edge, 'ColorBinding','interpolated', 'ColorData',cd)
end


lat=ocean.lat(ocean.sample_points);
lon=ocean.lon(ocean.sample_points);
scatterm(lat,lon,25,clrs,'filled');
scatterm(lat,lon,25,'k');
title(t)



sname=[pathname input_filename '/floats.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')