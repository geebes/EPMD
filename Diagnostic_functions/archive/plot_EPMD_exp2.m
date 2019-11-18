clear
clc

%% Prepare TM and grid metadata
disp('Loading MITgcm ECCO v4 Transport Matrix data')
load('~/GitHub/EPMD/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/Matrix13/Data/boxes.mat');
load('~/GitHub/EPMD/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/Matrix13/TMs/matrix_nocorrection_annualmean.mat')
load('~/GitHub/EPMD/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/grid.mat');
load('~/GitHub/EPMD/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/GCM/basin_mask.mat')
grid_load('~/GitHub/EPMD/Transport_Matrices/MITgcm/MITgcm_ECCO_v4_data/Data_v1/nctiles_grid/',5,'nctiles')
gcmfaces_global
dt_sec  = 60*60*6; % 6 hours
nsubtime=4;
is=1;

load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);
load ~/GitHub/EPMD/Transport_Matrices/MITgcm/MITgcm_ECCO_v4_data/GUD_spec_abundance.mat
load ~/GitHub/EPMD/Transport_Matrices/MITgcm/MITgcm_ECCO_v4_data/Theta.mat

        
Ib = find(izBox==1); % surface boundary
Ii = find(izBox~=1); % interior
i=ixBox;
j=iyBox;
k=izBox;
volume=volb(Ib);
iface=boxfacenum(Ib);
ix=ixBox(Ib);
iy=iyBox(Ib);
iz=izBox(Ib);

xlon=Xboxnom(Ib); % full surface grid
ylat=Yboxnom(Ib); % full surface grid

bmask=zeros(1,nb);
bathy=zeros(1,numel(Ib));
for i=1:5
    ii=sub2ind(size(bathy_basin_mask{i}),ixBoxFace{i},iyBoxFace{i},izBoxFace{i});
    jj=boxnumglob{i};
    bmask(jj)=bathy_basin_mask{i}(ii);

end
bmask=bmask(Ib);
bmask(ylat<-40)=5;
basin_names{5}='Southern';
bmask(find(bmask'==1 & ylat>30 & ylat<47 & xlon<40))=6;
basin_names{6}='Mediterranian';

%  get water column depth 
clear depth tmp
depth=[];
for i=1:5
    tmp=sum(dz{i},3);
    tmp=tmp(find(tmp));
    depth=[depth;tmp];
end


% z(find(ideep{1}))

diag_fcns=diagnostic_functions;
%% Transport matrix
I = speye(nb,nb);
A = I+dt_sec*Aexpms;
dt_rat = dt_sec/deltaT;

% Pre-process TMs
disp('Correcting TM for mass conservation')
Af = no_negatives(Aexpms,volb);
Bf = extract_surf(Af,Ib,Ii,volb);

B=Bf.*dt_sec;%.*1e-3;
I=speye(size(B));
B=B+I;


xn=B'*xlon;
% xn(xn>180)=xn(xn>180)-360;
yn=B'*ylat;
u=xn-xlon;
v=yn-ylat;



%% 

[ matObj ] = diag_fcns.open_matobj('Experiment_2/EPMD_output_old.mat');

T_opt=diag_fcns.load_variable(matObj,'T_optima',1);
Temperature=diag_fcns.load_variable(matObj,'Temperature',1);

tseries_lon=diag_fcns.load_variable(matObj,'tseries_lon',1);
tseries_lat=diag_fcns.load_variable(matObj,'tseries_lat',1);

iyr = diag_fcns.last_year(matObj)

x = diag_fcns.load_variable(matObj,'x',iyr);

nphen=numel(unique(T_opt));
nxr=numel(Ib);
nxc=size(x,2);
nspc=nxc./nphen;

nday=365;

% reshape for each seed population 
x3=reshape(full(x),nxr,nphen,nspc);
% takes original [nxc*nxc] matrix
% reads columnwise (1 column = global grid for 1 phenotype/species)
% reads nxr locations into dimension 1
% reads nphen phenotypes into dimension 2
% reads nspc ancestral species into dimension 3


% sum across all phenotypes in each 'species'
x2=squeeze(sum(x3,3));

K = diag_fcns.load_variable(matObj,'CarryingCapacity',1);
X = x.*K(:,end);
abundance=X;
frequency=x;

% load timescales
if ~isempty(cell2mat(matObj.T_occupied))
    t1 = diag_fcns.load_variable(matObj,'t_occupied',1);
    tN=t1;
    tN(tN==0)=NaN;
end
% load index of seed points
load siteindex.mat

% if ~isempty(cell2mat(matObj.RGB(1,1)))
%     tmp=diag_fcns.load_variable(matObj,'RGB',iyr);
%     rgb=reshape(full(tmp),nxr,nxc,3);
% end

clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot global distribution of each ancestral species

if true
    f100=figure(100);
    clf
    f100.Position=[52 1 2509 1344];
    
        
    c=0;
    for N=1:nspc
        c=c+1;
        subplot(9,9,c)
        
        
        ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
        ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
        
        xx=squeeze(x3(:,:,N)); % extract ancestral species N (dimension 3)
        typeN=sum(xx,2);       % sum across phenotypes
        
        fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
        [lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
        Xmap(Xmap==0)=NaN;
        Xmap = (Xmap);
        sh=surfacem(lat, lon, Xmap);
        caxis([0 1]);
%         
%         contr=sum(xx>0,2);
%         fld=vector2gcmfaces(contr,iface,ix,iy,iz);
%         [lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
%         Xmap(Xmap==nphen+1)=NaN; % set uninhabited locations to NaN
%         Xmap(Xmap==0)=NaN;
%         surfacem(lat, lon, Xmap); % map
%         caxis([0 77]);

        geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
        
        pos=[0.1    0.7673    0.75    0.1577];
        
        th=title(['T_{ancestral} = ' num2str(T_opt(N)) '^\circC']);
        th.Position(1:2)=th.Position(1:2)-[2 0.65];
        axis off;
        drawnow
    end
    
%     set(gcf,'Color','none');
%     
%     fname=['figures/Topt_abundances.png'];
%      export_fig(fname,'-r200')
    
end


%% Plot global distribution of each phenotype
if true
    f101=figure(101);
    clf
    f101.Position=[52 1 2509 1344];
            
    c=0;
    for N=1:nspc
        c=c+1;
        subplot(9,9,c)
        
        ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
        ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
        
        xx=squeeze(x3(:,N,:)); % extract phenotype N (dimension 2)
        typeN=sum(xx,2);       % sum across ancestral species
        
        fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
        [lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
        Xmap(Xmap==0)=NaN;
        Xmap = log10(Xmap);
        sh=surfacem(lat, lon, Xmap);
        geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
        caxis([-10 0]);
%         caxis([0 1]);
        
        pos=[0.1    0.7673    0.75    0.1577];
        
        th=title(['T_{opt} = ' num2str(T_opt(N)) '^\circC']);
        th.Position(1:2)=th.Position(1:2)-[2 0.65];
        axis off;
        drawnow
    end
    
    
%     set(gcf,'Color','none');
%     
%     fname=['figures/Topt_abundances.png'];
%      export_fig(fname,'-r200')
    
end

%% Plot different phenotypes of one example ancestral species
if true
    f102=figure(102);
    clf
    f102.Position=[52 1 2509 1344];
    
    iseed=39;
    
    for N=1:nxc/nphen
        subplot(9,9,N)
        
        ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
        ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
        
        typeN=squeeze(x3(:,N,iseed));
        fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
        [lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
        Xmap(Xmap==0)=NaN;
        Xmap = log10(Xmap);
        sh=surfacem(lat, lon, Xmap);
        geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
        caxis([-6 0]);
%         caxis([0 0.1]);
        
        pos=[0.1    0.7673    0.75    0.1577];
        
        th=title(['T_{opt} = ' num2str(T_opt(N)) '^\circC']);
        th.Position(1:2)=th.Position(1:2)-[2 0.65];
        if N==iseed
            th.Color='r';
        end
        axis off;
        drawnow
    end
    
    
    
%     set(gcf,'Color','none');
%     
%     fname=['figures/Topt_abundances.png'];
%     export_fig(fname,'-r200')
    
end

%% Diversity
figure(111)
threshold=1e-12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functional diversity (diversity of different phenotypes in each location)
typeN   =  squeeze(sum(x3,2)); % sum across ancestral species (dimension 2))
shannon = -nansum(typeN.*log(typeN),2);
typeN   =  sum(typeN>threshold,2); % count number of nonzero phenotypes in each location

subplot(321)
fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
Xmap(Xmap==0)=NaN;
sh=surfacem(lat, lon, Xmap);
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
title('Number of coexisting phenotypes')
colorbar

subplot(322)
fld=vector2gcmfaces(shannon,iface,ix,iy,iz);
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
Xmap(Xmap==0)=NaN;
sh=surfacem(lat, lon, Xmap);
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
title('Shannon diversity of phenotypes')
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taxonomic diversity (diversity of different ancestral species in each location)
typeN=squeeze(sum(x3,3)); % sum across phenotypes (dimension 3))
shannon = -nansum(typeN.*log(typeN),2);
typeN=sum(typeN>threshold,2); % count number of nonzero ancestral species in each location

subplot(323)
fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
Xmap(Xmap==0)=NaN;
Xmap = (Xmap);
sh=surfacem(lat, lon, Xmap);
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
colorbar
title('Number of coexisting ancestral species')

subplot(324)
fld=vector2gcmfaces(shannon,iface,ix,iy,iz);
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
Xmap(Xmap==0)=NaN;
sh=surfacem(lat, lon, Xmap);
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
title('Shannon diversity of ancestral species')
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taxonomic diversity (diversity of different ancestral species in each location)
shannon = -nansum(x.*log(x),2);
typeN=squeeze(sum(x>threshold,2)); 

subplot(325)
fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
Xmap(Xmap==0)=NaN;
Xmap = (Xmap);
sh=surfacem(lat, lon, Xmap);
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
colorbar
title('Total diversity')

subplot(326)
fld=vector2gcmfaces(shannon,iface,ix,iy,iz);
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
Xmap(Xmap==0)=NaN;
sh=surfacem(lat, lon, Xmap);
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
title('Total shannon diversity')
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Time series
figure(4)
tseries_lon = [-064 -158 -145 +062 -140 -180 -019 +068 -170 ];
tseries_lat = [+032 +023 +050 +016 +000 -076 +047 -051 -61.5];
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
scatterm(tseries_lat,tseries_lon,10,'r','filled')
for i=1:9
    textm(tseries_lat(i),tseries_lon(i),num2str(i))
end
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
drawnow


xyears=zeros(9,nxc,iyr*nday);
for yr=1:iyr
    days=(yr-1).*365+(1:nday);
      xyears(:,:,  days)=diag_fcns.load_variable(matObj,'tseries_x',yr);
end
%%
figure(2)

xsite1=squeeze(xyears(1,:,:));
xsite=reshape(xsite1,nphen,nspc,[]); % reshape to separate phenotypes and ancestral species

minlogx=-15;
for t=1:10:max(days)
    cla
    contourf(T_opt,T_opt,log10(xsite(:,:,t))',minlogx:0)
    hold on
    plot([minmax(T_opt)],[minmax(T_opt)],'m-')
    caxis([minlogx 0])
    title((i./365))
    xlabel('phenotype')
    ylabel('ancestor')
    title(floor([t./365 rem(t,365)]))
    drawnow
end

%%


f3=figure(3);
f3.Position=[1000 682 1136 656];
nclr=16;

clf
clear xx yy zz

% set(gcf,'DefaultAxesColorOrder',lhsdesign(nphen,3))
cmap = repmat(flipud(redblue(nphen)),nphen,1);
set(gcf,'DefaultAxesColorOrder',cmap)

for i=1:9
    xsite1=squeeze(xyears(i,:,:));        % get local timeseries data
    xsite=reshape(xsite1,nphen,nspc,[]); % reshape to separate phenotypes and ancestral species
    
    xsum=squeeze(sum(xsite,1)); % sum across phenotypes (biomass in each 'species' through time)
    cumxsum=cumsum(xsum,1);
    
    ax=subplot(3,3,i);
    ax.Position(3:4)=ax.Position(3:4).*1.25;
    ax.Position(1)=ax.Position(1)-0.1;
    
    % plot succession of ancestral species
    bsline=zeros(1,size(xsite,3));
    for j=1:nphen
        tmp = squeeze(xsite(:,j,:));
        h=area([bsline;tmp]','EdgeColor','none');
        h(1).FaceColor='none';
        bsline=csxsum(j,:);
        hold on
        ylim([0 1])
        drawnow
    end
    
    % plot phenotypes
    hold on
    csxsite1=cumsum(xsite1)';
    h=plot(csxsite1,'LineWidth',0.5);
    set(h, {'color'},num2cell(cmap,2));
    % 
%     plot(cumsum(xsite1)',':','color',[.5 .5 .5],'LineWidth',0.25);
    plot(cumsum(xsum)','color',[.5 .5 .5],'LineWidth',1)
    
    set(gca,'XTick',[],'YTick',[])
    axis([[0 iyr].*nday 0 1])
    
    drawnow
end

% fname=['figures/RGB_timeseries_sites_' num2str(iyr) '.png'];
% export_fig(fname,'-r200')









