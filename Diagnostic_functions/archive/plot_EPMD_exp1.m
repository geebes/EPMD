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

[ matObj ] = diag_fcns.open_matobj('Experiment_0/circulation_only/EPMD_output.mat');

T_opt=diag_fcns.load_variable(matObj,'T_optima',1);
Temperature=diag_fcns.load_variable(matObj,'Temperature',1);

tseries_lon=diag_fcns.load_variable(matObj,'tseries_lon',1);
tseries_lat=diag_fcns.load_variable(matObj,'tseries_lat',1);

iyr = diag_fcns.last_year(matObj);

tmp= diag_fcns.load_variable(matObj,'x',iyr);

nphen=numel(unique(T_opt));
nxr=numel(Ib);
nxc=size(tmp,2);

nday=365;
%
x  = reshape(tmp,nxr,nxc);

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


%%
if true
    f102=figure(102);
    clf
    f102.Position=[52 1 2509 1344];
        
    c=0;
    for N=1:nxc
        c=c+1;
        subplot(9,9,c)
        
        ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
        ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.01;
        
        typeN=x(:,N);
        fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
        [lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
        Xmap(Xmap==0)=NaN;
%         Xmap = log10(Xmap);
        Xmap = Xmap;
        sh=surfacem(lat, lon, Xmap);
        geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
%         caxis([-6 0]);
        caxis([0 1]);
        
        pos=[0.1    0.7673    0.75    0.1577];
        
        th=title(['T_{opt} = ' num2str(T_opt(N)) '^\circC']);
        th.Position(1:2)=th.Position(1:2)-[2 0.65];
        axis off;
        drawnow
    end
%     sp1=subplot(9,9,[nxc+1 81]);
%     ch1=colorbar('SouthOutside');
%     set(ch1,...
%             'XTick',log10([1 23 5 10 20 30]),...
%             'XTickLabel',{'1','2','3','5','10','20','30'})
%     ch1.Position=[0.5723 0.135 0.3326 0.02];
%     ch1.FontSize=15;
%     caxis([0 0.33]);
%     colormap(sp1,winter);
%     axis off
    
    
    
%     set(gcf,'Color','none');
%     
%     fname=['figures/Topt_abundances.png'];
%     export_fig(fname,'-r200')
    
end

%% Time series

f3=figure(3);
f3.Position=[1000 682 1136 656];
nclr=16;

xyears=zeros(9,nxc,iyr*nday);
for yr=1:iyr
    days=(yr-1).*365+(1:nday);
      xyears(:,:,days)=diag_fcns.load_variable(matObj,'tseries_x',yr);
end

[~,ilatsrt]=sort(abs(ylat))

clf
clear xx yy zz

cmap=hot(nxc);
set(gcf,'DefaultAxesColorOrder',cmap)

for i=1:9
    xsite=squeeze(xyears(i,:,:));
    
%     xsum=squeeze(sum(xsite,1));

    ax=subplot(3,3,i);
    ax.Position(3:4)=ax.Position(3:4).*1.25;
    ax.Position(1)=ax.Position(1)-0.1;
    
    % plot phenotypes
    h=area(xsite','EdgeColor',[1 1 1].*0.75,'LineWidth',1)
%     set(gca,'XScale','log')
    
%     set(gca,'XTick',[],'YTick',[])
%     xlim([95 99].*nday)
%     ylim([0 1])

    axis tight
    
    
    drawnow
end

% fname=['figures/RGB_timeseries_sites_' num2str(iyr) '.png'];
% export_fig(fname,'-r200')








