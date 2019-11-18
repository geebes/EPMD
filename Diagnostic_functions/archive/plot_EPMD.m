clear
clc

%% Prepare TM and grid metadata
disp('Loading MITgcm ECCO v4 Transport Matrix data')
load('~/GitHub/EPMD/TM_data/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/Matrix13/Data/boxes.mat');
load('~/GitHub/EPMD/TM_data/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/Matrix13/TMs/matrix_nocorrection_annualmean.mat')
load('~/GitHub/EPMD/TM_data/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/grid.mat');
load('~/GitHub/EPMD/TM_data/Transport_Matrices/MITgcm/MITgcm_ECCO_v4/GCM/basin_mask.mat')
grid_load('~/GitHub/EPMD/TM_data/Transport_Matrices/MITgcm/MITgcm_ECCO_v4_data/Data_v1/nctiles_grid/',5,'nctiles')
gcmfaces_global
dt_sec  = 60*60*6; % 6 hours
nsubtime=4;
is=1;

load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);
load ~/GitHub/EPMD/TM_data/Transport_Matrices/MITgcm/MITgcm_ECCO_v4_data/GUD_spec_abundance.mat
load ~/GitHub/EPMD/TM_data/Transport_Matrices/MITgcm/MITgcm_ECCO_v4_data/Theta.mat

        
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

[ matObj ] = diag_fcns.open_matobj('EPMD_output.mat');

T_opt=diag_fcns.load_variable(matObj,'T_optima',1);
Temperature=diag_fcns.load_variable(matObj,'Temperature',1);

tseries_lon=diag_fcns.load_variable(matObj,'tseries_lon',1);
tseries_lat=diag_fcns.load_variable(matObj,'tseries_lat',1);

iyr = diag_fcns.last_year(matObj);
% iyr=9

tmp= diag_fcns.load_variable(matObj,'x',iyr);
nxr=numel(Ib);
nxc=size(tmp,2);
nday=365;

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



if ~isempty(cell2mat(matObj.RGB(1,1)))
    tmp=diag_fcns.load_variable(matObj,'RGB',iyr);
    rgb=reshape(full(tmp),nxr,nxc,3);
end

clc

%%

f1=figure(1);
clf
f1.Position=[250 630 1502 715];

ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);

N=60;
typeN=frequency(:,N);
fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
Xmap(Xmap==0)=NaN;
sh=surfacem(lat, lon, log10(Xmap));
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
caxis([-10 0]);
colormap(ax1,winter);
pos=[0.1    0.7673    0.75    0.1577];

axis off;
ch1=colorbar('SouthOutside','FontSize',48);


set(gcf,'Color','none');

fname=['figures/X_' num2str(T_opt(N)) '_' num2str(iyr) '.png'];
export_fig(fname,'-r200')

%%
if false
    f101=figure(101);
    clf
%     f101.Position=[52 1 2509 1344];
        
    ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
    ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.005;
    
    t_emmigrant=prctile(t1,95,1);
    t_immigrant=prctile(t1,95,2);
    
    fld=vector2gcmfaces(t_immigrant,iface,ix,iy,iz);
    [lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
    Xmap(Xmap==0)=NaN;
    Xmap = log10(Xmap);
    sh=surfacem(lat, lon, Xmap);
    geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    
    scatterm(ylat(iK),xlon(iK),50,'k','filled')
    scatterm(ylat(iK),xlon(iK),40,log10(t_emmigrant),'filled')
    
    caxis([log10([3 100])]);
    colormap(ax1,flipud(parula));
%     pos=[0.1    0.7673    0.75    0.1577];
%     
%     
%     %         th=title(['T_{opt} = ' num2str(T_opt(N)) '^\circC']);
%     %         th.Position(1:2)=th.Position(1:2)-[2 0.65];
%     axis off;
%     drawnow
%     
    ch1=colorbar('SouthOutside');
    set(ch1,...
            'XTick',log10([3 5 10 20 30 50 100]),...
            'XTickLabel',{'3','5','10','20','30','50','100'})
%     ch1.Position=[0.5723 0.135 0.3326 0.02];
    ch1.FontSize=15;
%     caxis([0 0.33]);
% %     colormap(sp1,winter);
%     axis off
%     
%     set(gcf,'Color','none');
%     
%     fname=['figures/Dispersal_Times.png'];
%     export_fig(fname,'-r200')
    
end

%%
% figure(66)
% clf
% [N X]=hist(reshape(tN(:),[],1),10000);
% loglog(X,numel(t1)-cumsum(N),'LineWidth',2)
% hold on
% 
% axis([1/365 100 1 1e8])
% 
% set(gca,'XTick',[[1 7 30]./365 1 10 100],...
%         'XTickLabel',{'day','week','month','year','decade','century'},...
%         'XMinorTick','off','XGrid','on','XMinorGrid','off',...
%         'YGrid','off')

%%
if true
    f102=figure(102);
    clf
    f102.Position=[52 1 2509 1344];
    
    [~,ii]=sortrows([ylat(iK) xlon(iK)],[1 2]);
    ii=ii(round(linspace(1,nxc,64)))';
    c=0;
    for N=ii
        c=c+1;
        subplot(8,8,c)
        
        ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
        ax1.Position=ax1.Position+[-0.1 -0.2 2 4].*0.0125;
        
        typeN=t1(:,N);
        fld=vector2gcmfaces(typeN,iface,ix,iy,iz);
        [lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);
        Xmap(Xmap==0)=NaN;
        Xmap = log10(Xmap);
        sh=surfacem(lat, lon, Xmap);
        geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
        caxis(log10([1/12 30]));
        colormap(ax1,flipud(parula));
        pos=[0.1    0.7673    0.75    0.1577];
        
        scatterm(ylat(iK(N)),xlon(iK(N)),10,'r','filled')
        
%         th=title(['T_{opt} = ' num2str(T_opt(N)) '^\circC']);
%         th.Position(1:2)=th.Position(1:2)-[2 0.65];
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
%%
f2=figure(2);
f2.Position=[361 630 1502 715];
clf

ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);

cscale=max(abs(rgb(:)));
rgbN=(1+rgb./cscale)./2;                 % normalise

[~,dominant]=max(X,[],2);
clear rgbdom
for i=1:size(x,2)
    ii=find(dominant==i);
    rgbdom(ii,:)=squeeze(rgbN(ii,i,:));
end

clear rgbclr
for i=1:3
    fld=vector2gcmfaces(rgbdom(:,i),iface,ix,iy,iz);
    [lon lat rgbclr(:,:,i)]=convert2pcol(mygrid.XC,mygrid.YC,fld);
end

np=numel(lat);
geoshow(lat,lon,rgbclr);
% s=scatterm(lat(:)+randn(np,1)/10,lon(:)+randn(np,1)./10,50,reshape(rgbclr,[],3),'filled','o');
% s.Children.MarkerFaceAlpha = .5;

cntr=vector2gcmfaces(dominant,iface,ix,iy,iz);
[lon lat cntr]=convert2pcol(mygrid.XC,mygrid.YC,cntr);
lon(isnan(lon))=-999;lat(isnan(lat))=-999;cntr(isnan(cntr))=-1e9;
contourm(lat,lon,cntr,0:5:75,'color',ones(1,3).*1.0,'LineWidth',0.25);

scatterm(tseries_lat,tseries_lon,150,'r')
th=textm(tseries_lat,tseries_lon,{'1','2','3','4','5','6','7','8','9',},'HorizontalAlignment','center');

%% Time series

%%
f3=figure(3);
f3.Position=[1000 682 1136 656];
nclr=16;

  xyears=zeros(9,nxc  ,iyr*nday);
RGByears=zeros(9,nxc,3,iyr*nday);
for yr=max(1,iyr-4):iyr
    days=(yr-1).*365+(1:nday);
      xyears(:,:,  days)=diag_fcns.load_variable(matObj,'tseries_x',yr);
    RGByears(:,:,:,days)=diag_fcns.load_variable(matObj,'tseries_RGB',yr);
end
RGByears=permute(RGByears,[1 2 4 3]);
RGByears=(1+RGByears./cscale)./2;                 % normalise

clf
clear xx yy zz
for i=1:9
    ax=subplot(3,3,i);
    ax.Position(3:4)=ax.Position(3:4).*1.25;
    ax.Position(1)=ax.Position(1)-0.1;
    
    xx=repmat(padarray(1:iyr.*nday,[0 1],NaN,'post'),1,nxc);
    yy=reshape(padarray(squeeze(xyears(i,:,:))',[1 0],NaN,'post'),1,[]);
    zz(1,:,1)=reshape(padarray(squeeze(RGByears(i,:,:,1))',[1 0],NaN,'post'),1,[]);
    zz(1,:,2)=reshape(padarray(squeeze(RGByears(i,:,:,2))',[1 0],NaN,'post'),1,[]);
    zz(1,:,3)=reshape(padarray(squeeze(RGByears(i,:,:,3))',[1 0],NaN,'post'),1,[]);
%     scatter(xx,yy,15,zz,'filled')
    surface([xx;xx],[yy;yy],zeros(2,size(zz,2)),[zz;zz],...
        'facecol','no',...
        'edgecol','flat',...
        'linew',3);
    
    hold on
%     plot(1:nday,Temperature(tseries_ind(i),:)./100,'k','LineWidth',2)
%     set(gca,'YScale','log')

    box on
    set(gca,'XTick',[],'YTick',[])
    xlim([max(1,iyr-4) iyr].*nday)
end

fname=['figures/RGB_timeseries_sites_' num2str(iyr) '.png'];
export_fig(fname,'-r200')

%%

f4=figure(4);
f4.Position=[950 682 1136 656];

cmp=colormap(jet(size(x,2)));
set(gcf,'DefaultAxesColorOrder',cmp)

colormap(cmp);
clf
for i=1:9
    ax=subplot(3,3,i);
    ax.Position(3:4)=ax.Position(3:4).*1.25;
    ax.Position(1)=ax.Position(1)-0.1;
    
    xx=repmat(padarray(1:iyr.*nday,[0 1],NaN,'post'),1,nxc);
    yy=reshape(padarray(squeeze(  xyears(i,:,:))',[1 0],NaN,'post'),1,[]);
    zz=reshape(padarray(repmat(T_opt,iyr.*nday,1),[1 0],NaN,'post'),1,[]);
%     scatter(xx,yy,15,zz,'filled')
    surface([xx;xx],[yy;yy],zeros(2,size(zz,2)),[zz;zz],...
        'facecol','no',...
        'edgecol','flat',...
        'linew',3);
    caxis([-2 36]);
    
    box on
    set(gca,'XTick',[],'YTick',[])
    xlim([max(1,iyr-4) iyr].*nday)
end

fname=['figures/Toptima_timeseries_sites_' num2str(iyr) '.png'];
export_fig(fname,'-r200')
    
%%
f5=figure(5);
f5.Position=[850 682 1500 656];

clf
set(gcf,'DefaultAxesColorOrder',cmp)
for i=1:9
    subplot(3,3,i)
    h=bar(squeeze(xyears(i,:,:))',1,'stacked');
    hold on
    plot(cumsum(squeeze(xyears(i,:,:)),1)','k','LineWidth',1)

    axis([[max(1,iyr-4) iyr].*nday 0 1])

end

fname=['figures/Toptima_fractionation_sites_' num2str(iyr) '.png'];
export_fig(fname,'-r200')
%%
f6=figure(6);
f6.Position=[800 682 1500 656];

clf
set(gcf,'DefaultAxesColorOrder',cmp)

for i=1:9
    ax=subplot(3,3,i);
    ax.Position(3:4)=ax.Position(3:4).*1.25;
    ax.Position(1)=ax.Position(1)-0.1;
    
    stacks  =padarray(cumsum(squeeze(xyears(i,:,:)),1)',[0 1],'pre')';
    RGBsite =squeeze(RGByears(i,:,:,:));
    RGBsite =permute(repmat(RGBsite,1,1,1,2),[4 1 2 3]);
    
    for j=1:nxc
        xx=[1:iyr.*nday;1:iyr.*nday];
        yy=[stacks(j,:);stacks(j+1,:)];
        zz=zeros(2,iyr.*nday);
        cc=squeeze(RGBsite(:,j,:,:));
        surface(xx,yy,zz,cc,'edgecol','none');
        hold on
    end
    
    plot(cumsum(squeeze(xyears(i,:,:)),1)','LineWidth',0.5);
    axis tight
    box on
    set(gca,'XTick',[],'YTick',[])
    xlim([max(1,iyr-4) iyr].*nday)
    title(i)
    drawnow
end


fname=['figures/RGB_fractionation_sites_' num2str(iyr) '.png'];
export_fig(fname,'-r200')

%% subsample grid points

if true % if true, even spacing
         % if false, Tara
    % generate equally spaced points on surface of a sphere
    mesh='quad';
    switch mesh
        case 'tri'
            TR=SubdivideSphericalMesh(IcosahedronMesh,3);
            xyz=TR.X;
        case 'quad'
            fv=SubdivideSphericalMesh(QuadCubeMesh,2);
            xyz=fv.vertices;
    end
    % convert to spherical coordinates
    [theta,phi,~]=cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
    % convert to degrees
    xs=rad2deg(theta);
    ys=rad2deg(phi);
    % adjust longitude range
    xs(xs<0)=xs(xs<0)+360;
else
    load('Tara_sites.mat');
    xs = tara_latlon(:,2);
    ys = tara_latlon(:,1);
end
% find sample points nearest to each model gridpoint
Kall = dsearchn([xs ys],[xlon ylat]);
% remove duplicates
[K,~,Kall] = unique(Kall);
% xs(K) & ys(K) are grid coordinates nearest to sample grid
% find index of those coordinates in model grid
K = dsearchn([xlon ylat],[xs(K),ys(K)]);

% remove all sites where water column depth less than maxd metres
maxd=200;
K=K(depth(K)>maxd);

xk=xlon(K);
yk=ylat(K);

nK=size(K,1);

%% Calculate and draw phylogeny

f7=figure(7);
f7.Position=[50 773 975 572];
clf

[tree D D2 t_index c_index] = diag_fcns.make_tree(matObj,iyr,1e6,K);
T_optima=repmat(T_opt,numel(Ib),1);

nclustrs=nxc;
[ taxID, cmap, outperm, Tclust, RGBclust, H, Dxy] = ...
    diag_fcns.draw_tree(tree,D,t_index,nclustrs,T_optima,rgbN,X);
xl=get(gca,'XLim');





% Draw genetic distance matrix
ax4=subplot(133);
ax4.Position=[0.4200 0.0680 0.4000 0.8990];

Dsq=squareform(D);
% Dxy=1:size(Dsq,1);
h=pcolor(Dxy,Dxy,Dsq(outperm,outperm));
colormap(parula)
caxis(xl);
h.EdgeAlpha=0.1;
axis square
hold on
axs=axis;
shading flat

clear cc Cbasin
cc(:,1,:)=cmap;              % add extra dimension 
bmap=dec2bin(1:6,3)-'0';    % index for ocean basin colours
ii=repmat(bmask',1,nxc);    % pad array to use t_index
tmp=bmap(ii(t_index(outperm)),:); % extract basin IDs
Cbasin(:,1,:)=tmp;              % add extra dimension 


% Matrix colorbars
axis tight
ax=axis;
axmax=ax(2);
axmin=ax(1);
axrng=diff(ax(1:2));
barwdth=axrng/20;
[xg yg]=meshgrid([0 barwdth],Dxy);
zg=yg.*0;
surf(-xg,yg         ,zg,cc      ,'EdgeAlpha',0.01)
surf(xg+axmax   ,yg         ,zg,Cbasin  ,'EdgeAlpha',0.01)
surf(yg         ,-fliplr(xg),zg,RGBclust,'EdgeAlpha',0.01)
surf(yg         ,xg+axmax   ,zg,Tclust  ,'EdgeAlpha',0.01)
axis square
axis tight
shading flat
axis ij
ax4.XTick=[];
ax4.YTick=[];
axis off
rectangle('Position',[axmin axmin axrng axrng],'LineWidth',0.75)
rectangle('Position',[axmin-barwdth axmin axrng+2*barwdth axrng],'LineWidth',0.75)
rectangle('Position',[axmin axmin-barwdth axrng axrng+2*barwdth],'LineWidth',0.75)
% ax4.Position=[0.5    0.0464  0.3003    0.9075];
text(axmax/2,axmax*22/20,'Optimal Temperature','HorizontalAlign','c','FontSize',12)
text(axmax/2,-axmax*2/20,'RGB gene','HorizontalAlign','c','FontSize',12)
text(axmax*22/20,axmax/2,'Ocean Basin','HorizontalAlign','c','FontSize',12,'Rotation',90)

% mark out top nKclusters
nKclusters = 5;
Kclust = cluster(tree,'maxclust',nKclusters,'Criterion','distance');
for ik=1:nKclusters
    linelimits=[axmin-barwdth axmax+barwdth];
    plot(linelimits,[1 1].*min(min(yg(Kclust(outperm)==ik,:))),'k-',...
         'LineWidth',1)
    plot([1 1].*min(min(yg(Kclust(outperm)==ik,:))),linelimits,'k-',...
         'LineWidth',1)
end
     
% cbh=colorbar('SouthOutside');
% ax4.Colormap=flipud(parula(64));
% cbh.Position=[ 0.25 0.2 0.9 0.55];


% rectangle('Position',[axmin+axrng.*[0.0075 0.6] axrng.*[0.3 0.3925]],...
%           'EdgeColor','none',...
%           'FaceColor',[1 1 1 0.66],...
%           'Curvature',0.1)
% for i=1:5
%     text(axmax/25,0.6.*axmax+i*axmax./15,basin_names{i},'color',bmap(i,:),'FontSize',18,'FontWeight','bold')
% end

set(gcf,'Color','none');
fname=['figures/Phylogeny_' num2str(iyr) '.png'];
export_fig(fname,'-r200')
set(gcf,'Color','w');




%% Shortest path distance
% bi-directional graph of non-zero elements in transport matrix 
% Edges are weighted by the inverse values, such that stronger fluxes have higher weights
B2=spones(B); % binary transport matrix
% (i.e. connected or not - useful for shortest)

G=digraph(1./B,'OmitSelfLoops');
sp=distances(G,K,K);
sp0=sp./4/365;
spAll=triu(sp0)+tril(sp0);
spUp =triu(sp0)+triu(sp0)';
spLo =tril(sp0)'+tril(sp0);
sp2=(spAll+spAll')./2;


%% % Plot global dispersal timescales and clusters 
f22=figure(22);
f22.Position=[72 36 1515 1309];
clf

nDclusters=10;
Clrmap=lhsdesign(3,nDclusters,'iterations',1)';
Clrmap=[033 119 176
        240 142 059
        048 159 057
        212 042 045
        174 122 159
        139 086 076
        225 120 191
        127 127 179
        187 188 059
        031 189 204]./255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,3,1)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-245 115]);
vect=median(spAll,1);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
scatterm(yk,xk,50,vect,'filled')
ch=colorbar('SouthOutside');caxis([0 10])
ch.TickLabels{end} = [ch.TickLabels{end} ' years'];
axis tight
axis off
title('(a) Median global dispersal time')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,3,2)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-245 115]);
vect=median(spAll,2);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
scatterm(yk,xk,50,vect,'filled')
ch=colorbar('SouthOutside');caxis([0 10])
ch.TickLabels{end} = [ch.TickLabels{end} ' years'];
axis tight
axis off
title('(b) Median global arrival time')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,3,3)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-245 115]);
vect=(median(spAll,1)+median(spAll,2)')./2;
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
scatterm(yk,xk,50,vect,'filled')
ch=colorbar('SouthOutside');caxis([0 10])
ch.TickLabels{end} = [ch.TickLabels{end} ' years'];
axis tight
axis off
title('(c) Median Bidirectional global connection time')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

treeUp = linkage(squareform(spUp),'average');
leafOrderUp = optimalleaforder(treeUp,spUp);
sp4=subplot(4,3,4);
cutoff = median([treeUp(end-(nDclusters-1),3) treeUp(end-(nDclusters-2),3)]);
clustID = cluster(treeUp,'MaxClust',nDclusters);
hd=dendrogram(treeUp,0,'ColorThreshold',cutoff,'ReOrder',leafOrderUp);
for i=1:numel(hd)
    if ~ismember(hd(i).Color,[0 0 0],'rows')
        ii=round(median(hd(i).XData));
        hd(i).Color=Clrmap(clustID(leafOrderUp(ii)),:);
    end
end
TdistUp = cluster(treeUp,'MaxClust',nDclusters);
ylim([0 10])
% axis off
sp4.Position=[0.1300 0.59 0.2134 0.1577];
title('(d) dispersal clusters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

treeLo = linkage(squareform(spLo),'average');
leafOrderLo = optimalleaforder(treeLo,spLo);
sp5=subplot(4,3,5);
cutoff = median([treeLo(end-(nDclusters-1),3) treeLo(end-(nDclusters-2),3)]);
clustID = cluster(treeLo,'MaxClust',nDclusters);
hd=dendrogram(treeLo,0,'ColorThreshold',cutoff,'ReOrder',leafOrderLo);
for i=1:numel(hd)
    if ~ismember(hd(i).Color,[0 0 0],'rows')
        ii=round(median(hd(i).XData));
        hd(i).Color=Clrmap(clustID(leafOrderLo(ii)),:);
    end
end
TdistLo = cluster(treeLo,'MaxClust',nDclusters);
ylim([0 10])
% axis off
sp5.Position=[0.4108 0.59 0.2134 0.1577];
title('(e) Arrival clusters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tree2 = linkage(squareform(sp2),'average');
leafOrder2 = optimalleaforder(tree2,sp2);
sp6=subplot(4,3,6);
cutoff = median([tree2(end-(nDclusters-1),3) tree2(end-(nDclusters-2),3)]);
clustID = cluster(tree2,'MaxClust',nDclusters);
hd=dendrogram(tree2,0,'ColorThreshold',cutoff,'ReOrder',leafOrder2);
for i=1:numel(hd)
    if ~ismember(hd(i).Color,[0 0 0],'rows')
        ii=round(median(hd(i).XData));
        hd(i).Color=Clrmap(clustID(leafOrder2(ii)),:);
    end
end
Tdist2 = cluster(tree2,'MaxClust',nDclusters);
ylim([0 10])
% axis off
sp6.Position=[0.691 0.59 0.2134 0.1577];
title('(f) Bidirectional connectivity clusters')
set(gca,'XDir','Reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp7=subplot(4,3,7);
imagesc(spUp(leafOrderUp,leafOrderUp));
axis square
hold on
distanceBoundary=find(diff(TdistUp(leafOrderUp)));
axis tight
xl=xlim;
yl=ylim;
clear clrs
clrs(:,1,:)=bmap(bmask(K(leafOrderUp)),:);              % add extra dimension
surf([xl(2)+[0 numel(K)./20]].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs);
shading flat
clear clrs
tmp=flipud(redblue(nclustrs));
clrs(:,1,:)=tmp(discretize(Temperature(K(leafOrderUp),end),nclustrs),:); % add extra dimension
surf([-numel(K)./20 xl(1)].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs);
clear clrs
clrs(:,1,:)=Clrmap(TdistUp(leafOrderUp),:);              % add extra dimension
surf((1:numel(K))'*[1 1]-0.5,[xl(1)+[0 -numel(K)./20]].*ones(size(K,1),2),K*[0 0],clrs);
shading flat
ch=colorbar('SouthOutside');caxis([0 15]);
for i=1:nDclusters-1
    plot(xlim,[distanceBoundary(i) distanceBoundary(i)]+0.5,'w-')
    plot([distanceBoundary(i) distanceBoundary(i)]+0.5,ylim,'w-')
end
rectangle('Position',[0 0 [1 1]*numel(K)],'LineWidth',1)
rectangle('Position',[-numel(K)./20 0 [1.1 1]*numel(K)],'LineWidth',1)
rectangle('Position',[0 -numel(K)./20 [1 1.05]*numel(K)],'LineWidth',1)
axis off
sp7.Position=sp4.Position([1 2 3 3])+[-0.029 -0.275 +0.058.*[1 1]];
ch.Position=ch.Position+[0 0.015 0 0];
ch.TickLabels{end} = [ch.TickLabels{end} ' years'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp8=subplot(4,3,8);
imagesc(spLo(leafOrderLo,leafOrderLo))
axis square
hold on
distanceBoundary=find(diff(TdistLo(leafOrderLo)));
axis tight
xl=xlim;
yl=ylim;
clear clrs
clrs(:,1,:)=bmap(bmask(K(leafOrderLo)),:);              % add extra dimension
surf([xl(2)+[0 numel(K)./20]].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs)
shading flat
clear clrs
tmp=flipud(redblue(nclustrs));
clrs(:,1,:)=tmp(discretize(Temperature(K(leafOrderLo),end),nclustrs),:);             % add extra dimension
surf([-numel(K)./20 xl(1)].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs)
clear clrs
clrs(:,1,:)=Clrmap(TdistLo(leafOrderLo),:);              % add extra dimension
surf((1:numel(K))'*[1 1]-0.5,[xl(1)+[0 -numel(K)./20]].*ones(size(K,1),2),K*[0 0],clrs)
shading flat
ch=colorbar('SouthOutside');caxis([0 15]);
for i=1:nDclusters-1
    plot(xlim,[distanceBoundary(i) distanceBoundary(i)]+0.5,'w-')
    plot([distanceBoundary(i) distanceBoundary(i)]+0.5,ylim,'w-')
end
rectangle('Position',[0 0 [1 1]*numel(K)],'LineWidth',1)
rectangle('Position',[-numel(K)./20 0 [1.1 1]*numel(K)],'LineWidth',1)
rectangle('Position',[0 -numel(K)./20 [1 1.05]*numel(K)],'LineWidth',1)
axis off
sp8.Position=sp5.Position([1 2 3 3])+[-0.029 -0.275 +0.058.*[1 1]];
ch.Position=ch.Position+[0 0.015 0 0];
ch.TickLabels{end} = [ch.TickLabels{end} ' years'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp9=subplot(4,3,9);
leafOrder2r=fliplr(leafOrder2);
imagesc(sp2(leafOrder2r,leafOrder2r))
axis square
hold on
distanceBoundary=find(diff(Tdist2(leafOrder2r)));
axis tight
xl=xlim;
yl=ylim;
clear clrs
clrs(:,1,:)=bmap(bmask(K(leafOrder2r)),:);              % add extra dimension
surf([xl(2)+[0 numel(K)./20]].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs);
shading flat
clear clrs
tmp=flipud(redblue(nclustrs));
clrs(:,1,:)=tmp(discretize(Temperature(K(leafOrder2r),end),nclustrs),:);             % add extra dimension
surf([-numel(K)./20 xl(1)].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs);
clear clrs
clrs(:,1,:)=Clrmap(Tdist2(leafOrder2r),:);              % add extra dimension
surf((1:numel(K))'*[1 1]-0.5,[xl(1)+[0 -numel(K)./20]].*ones(size(K,1),2),K*[0 0],clrs);
shading flat
ch=colorbar('SouthOutside');caxis([0 15]);
for i=1:nDclusters-1
    plot(xlim,[distanceBoundary(i) distanceBoundary(i)]+0.5,'w-')
    plot([distanceBoundary(i) distanceBoundary(i)]+0.5,ylim,'w-')
end
rectangle('Position',[0 0 [1 1]*numel(K)],'LineWidth',1)
rectangle('Position',[-numel(K)./20 0 [1.1 1]*numel(K)],'LineWidth',1)
rectangle('Position',[0 -numel(K)./20 [1 1.05]*numel(K)],'LineWidth',1)
axis off
sp9.Position=sp6.Position([1 2 3 3])+[-0.029 -0.275 +0.058.*[1 1]];
ch.Position=ch.Position+[0 0.015 0 0];
ch.TickLabels{end} = [ch.TickLabels{end} ' years'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,3,10)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-245 115]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
scatterm(yk,xk,50,Clrmap(TdistUp,:),'filled')
axis tight
axis off
title('(g) Dispersal clusters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,3,11)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-245 115]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
scatterm(yk,xk,50,Clrmap(TdistLo,:),'filled')
axis tight
axis off
title('(h) Arrival clusters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,3,12)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-245 115]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
scatterm(yk,xk,50,Clrmap(Tdist2,:),'filled')
axis tight
axis off
title('(i) Bidirectional clusters')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname=['figures/Transport_clusters_' num2str(iyr) '.png'];
export_fig(fname,'-r300')

%%
PhenoFreq = zeros(numel(K),nxc);
TaxonFreq = zeros(numel(K),nclustrs);
TaxonID   = zeros(numel(K),nclustrs);

tmp = frequency(t_index); % frequency at sample sites

iphen=sub2ind([numel(K) nxc],c_index(:,1),c_index(:,2)); % index for sample sites only
PhenoFreq(iphen)=tmp;

TaxonID(iphen)=taxID;

for i=1:numel(K) % for each sample site
    for j=find(TaxonID(i,:)) % for all extant phenotypes at that site
        TaxonFreq(i,TaxonID(i,j)) = TaxonFreq(i,TaxonID(i,j)) + PhenoFreq(i,j);  
    end
end
      




clear minTaxonFreq minPhenoFreq
for i=unique(taxID)' % for all genetic species
    % for each species, calculate minimum between all pairs of locations
    minTaxonFreq(:,:,i)=min(TaxonFreq(:,i),TaxonFreq(:,i)'); 
    minPhenoFreq(:,:,i)=min(PhenoFreq(:,i),PhenoFreq(:,i)'); 
end
Ctaxo=sum(minTaxonFreq,3);
Cphen=sum(minPhenoFreq,3);
% BC = sum of only the lesser counts for each species found in both sites
% We do not need to normalise by S_i+S_j because we are using normalised abundances
% S_i+S_j = 2, so BC=1-2.*C./(S_i+S_j) = 1-C;
BCtaxo=1-Ctaxo;
BCtaxo=BCtaxo.*(1-eye(size(BCtaxo))); % ensure diagonals are exactly zero
BCtaxo=(BCtaxo+BCtaxo')./2;
BCtaxo(BCtaxo<0)=0;
BCphen=1-Cphen;
BCphen=BCphen.*(1-eye(size(BCphen))); % ensure diagonals are exactly zero
BCphen=(BCphen+BCphen')./2;
BCphen(BCphen<0)=0;
% ALL THESE CORRECTIONS ARE A BIT DODGY
% SHOULD NOT NEED TO DO THEM!!!


%%

figure(33)
clf

binedges=0:0.5:10;
nbin=numel(binedges)-1;
[Nbin,edges,ibin] = histcounts(spAll(:),binedges);
BCdistTaxo=zeros(max(Nbin),nbin).*NaN;
BCdistPhen=zeros(max(Nbin),nbin).*NaN;
for i=1:nbin
    BCdistTaxo(1:Nbin(i),i)=BCtaxo(ibin==i);
    BCdistPhen(1:Nbin(i),i)=BCphen(ibin==i);
end

subplot(221)
boxplot(BCdistPhen,'outliersize',1,'symbol','.','color',ones(1,3).*0.25,'plotstyle','compact')
set(gca,'XTick',0.5:2:nbin+1,'XTickLabel',binedges(1:2:nbin+1))
xlabel('Distance (Years)')
title('Phenotypic BC dissimilarity')
subplot(222)
boxplot(BCdistTaxo,'outliersize',1,'symbol','.','color',ones(1,3).*0.25,'plotstyle','compact');
set(gca,'XTick',1:2:nbin+1,'XTickLabel',binedges(1:2:nbin+1))
xlabel('Distance (Years)')
title('Taxonomic BC dissimilarity')

subplot(223)
tree = linkage(squareform(BCphen),'average');
leafOrder = optimalleaforder(tree,BCphen);
pcolor(BCphen(leafOrder,leafOrder))
title('Phenotypic dissimilarity')
colorbar
axis square
hold on
clear clrs
clrs(:,1,:)=bmap(bmask(K(leafOrder)),:);              % add extra dimension
surf([-numel(K)./20 xl(1)].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs)
clear clrs
clrs(:,1,:)=Clrmap(Tdist2(leafOrder),:);              % add extra dimension
surf([xl(2)+[0 numel(K)./20]].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs)
clear clrs
tmp=flipud(redblue(nclustrs));
clrs(:,1,:)=tmp(discretize(Temperature(K(leafOrder),end),nclustrs),:)             % add extra dimension
surf((1:numel(K))'*[1 1]-0.5,[xl(2)+[0 numel(K)./20]].*ones(size(K,1),2),K*[0 0],clrs)
shading flat
axis tight
axis ij
axis off

subplot(224)
tree = linkage(squareform(BCtaxo),'average');
leafOrder = optimalleaforder(tree,BCtaxo);
pcolor(BCtaxo(leafOrder,leafOrder))
title('Taxonomic dissimilarity')
colorbar
axis square
hold on
clear clrs
clrs(:,1,:)=bmap(bmask(K(leafOrder)),:);              % add extra dimension
surf([-numel(K)./20 xl(1)].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs)
clear clrs
clrs(:,1,:)=Clrmap(Tdist2(leafOrder),:);              % add extra dimension
surf([xl(2)+[0 numel(K)./20]].*ones(size(K,1),2),(1:numel(K))'*[1 1]-0.5,K*[0 0],clrs)
clear clrs
tmp=flipud(redblue(nclustrs));
clrs(:,1,:)=tmp(discretize(Temperature(K(leafOrder),end),nclustrs),:)             % add extra dimension
surf((1:numel(K))'*[1 1]-0.5,[xl(2)+[0 numel(K)./20]].*ones(size(K,1),2),K*[0 0],clrs)
shading flat
axis tight
axis ij
axis off



%%
f8=figure(8);
f8.Position=[1301 1 1019 1344];
clf

normfreq=X./max(sum(X,2));
alp=normfreq(t_index(outperm)).*2;
rndx=randn(size(t_index)).*2;
rndy=randn(size(t_index)).*2;

ax2=subplot(311);
axm=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(axm,lat,lon,zeros(size(rgbclr.*0)))
geoshow(axm,lat,lon,rgbclr)
% hs=scatterm(ylat(K),xlon(K),40,'k');
% hs=scatterm(ylat(K),xlon(K),35,bmap(bmask(K),:),'filled');
geoshow(axm, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
axis off
ax2.Position=[0.01    0.66    1    0.33];

ax3=subplot(312);
axm=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(axm,lat,lon,zeros(size(lon))+0.9)
xx=repmat(xlon,[1 nxc]);
yy=repmat(ylat,[1 nxc]);
clear clr
for i=1:3
    tmp=rgbN(:,:,i);
    clr(:,i)=tmp(t_index(outperm))';
end
hs=scatterm(yy(t_index(outperm))+rndx,xx(t_index(outperm))+rndy,alp.*200,clr,'filled');
geoshow(axm, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
axis off
ax3.Position=[0.01    0.33    1    0.33];

ax4=subplot(313);
axm=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(axm,lat,lon,zeros(size(lon))+0.9)
scatterm(yy(t_index(outperm))+rndx,xx(t_index(outperm))+rndy,alp.*200,cmap,'filled');
geoshow(axm, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
axis off
ax4.Position=[0.01    0.00    1    0.33];

set(gcf,'Color','none');

fname=['figures/Phylogenetic_maps_' num2str(iyr) '.png'];
export_fig(fname,'-r200')

%% Bar chart (decoupling of function and taxonomy)
f10=figure(10);
f10.Position=[78 536 1100 809];
clf
set(gcf,'DefaultAxesColorOrder',flipud(redblue(nxc)));

% sort sample site by temperature
[~,I]=sort(Temperature(K,end));
nsite=numel(K);

[~, mreptuo] = sort(outperm); % inverse index for outperm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:nsite
    xcoord=repmat([i-1 i],nxc+1,1)';
    ycoord=cumsum(full(frequency(K(I([i i])),:)),2);
    ycoord=[[0;0] ycoord];
    
    sh1=subplot(311);
    % reset colour array
    clr   =rgbN(K(I(i)),:,:);
    clr   =repmat(clr,2,1,1);
    surf(xcoord,ycoord,zeros(size(xcoord)),clr);
    hold on
    
    sh2=subplot(312);
    % reset colour array
    clr=zeros(1,nxc,3).*NaN;
    ii=find(c_index(:,1)==I(i)); % index of terminal nodes at site i
    clr(1,c_index(ii,2),:) = cmap(mreptuo(ii),:); % get colours for site i
    clr   =repmat(clr,2,1,1);
    surf(xcoord,ycoord,zeros(size(xcoord)),clr);
    hold on
    
end


sh1=subplot(311)
view(2)
axis tight
title('Taxonomy (RGB)')
set(gca,'YTick',[0 1],'YTickLabel',{'0%','100%'},...
        'XTick',[])
sh1.Position(2)=sh1.Position(2)-0.05;
    
sh2=subplot(312);
view(2)
axis([0 nsite 0 1])
title('Taxonomy (OTU Clusters)')
set(gca,'YTick',[0 1],'YTickLabel',{'0%','100%'},...
        'XTick',[])    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sh3=subplot(313);
sh3.Position(2)=sh3.Position(2)+0.05;


bar(frequency(K(I),:),1,'stacked','EdgeColor','k')
axis tight
title('Phenotype (Optimal Temperature)')
set(gca,'YTick',[0 1],'YTickLabel',{'0%','100%'},...
        'XTick',[])
set(gca,'YTick',[0 1],...
        'YTickLabel',{'0%','100%'},...
        'YColor','k')

pos=sh3.Position;

colormap(flipud(redblue(nxc-1)))
caxis(T_opt([1 end]))
hc=colorbar('SouthOutside');
hc.Ticks=T_opt([1:2:end]);
sh3.Position=pos;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'Color','none');
fname=['figures/Decoupling_' num2str(iyr) '.png'];
export_fig(fname,'-r200')
set(gcf,'Color','w');


%%

richness_phen=sum(frequency~=0,2);

[I,J]=ind2sub([nxr nxc],t_index);
matr = sparse(I,taxID,frequency(t_index),nxr,nclustrs); % matrix showing frequency for each genetic cluster and sampled location
richness_geno=sum(matr(K,:)~=0,2);

fld=vector2gcmfaces(richness_phen,iface,ix,iy,iz);
[lon lat R]=convert2pcol(mygrid.XC,mygrid.YC,fld);

fld=vector2gcmfaces(Temperature(:,end),iface,ix,iy,iz);
[lon lat T]=convert2pcol(mygrid.XC,mygrid.YC,fld);
lon(isnan(lon))=-999;lat(isnan(lat))=-999;T(isnan(T))=-1e9;

f9=figure(9)
f9.Position=[1063 1 568 1344];
clf

subplot(411)
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
surfacem(lat, lon, R);
caxis([0 nclustrs])
colormap(ax1,winter);
caxis([0 60])
pos=[0.1    0.7673    0.75    0.1577]
axis off
ch1=colorbar('WestOutside');
title(ch1,'Functional Richness',...
          'HorizontalAlignment','l',...
          'FontSize',12)
ax1.Position=pos;

subplot(412)
ax2=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
contourm(lat, lon, T, [0:5:35],'-','LineWidth',0.5);
geoshow(ax2, land, 'FaceColor', [0.7 0.7 0.7]);
colormap(ax2,jet(7));
caxis([-2.5 32.5])
axis off
ch2=colorbar;
title(ch2,'Temperature (^\circC)',...
                      'HorizontalAlignment','r',...
                      'FontSize',12)
ax2.Position=pos;

ax3=subplot(412);
sh=scatter(Temperature(K,end),richness_phen(K),50,abs(ylat(K)),'filled');
sh.MarkerFaceAlpha=1;
sh.MarkerFaceAlpha=0.25;
sh.MarkerEdgeColor='k';
colormap(ax3,redblue)
% xlabel('Temperature ^\circC')
ylabel('Functional Richness')
colorbar
caxis([0 80])
ax3.Position=ax3.Position+[-0.05 0.03 -0.02 0];

% ax4=subplot(413);
% sh=scatter(Temperature(K,end),richness_geno,50,abs(ylat(K)),'filled');
% sh.MarkerFaceAlpha=1;
% sh.MarkerFaceAlpha=0.25;
% sh.MarkerEdgeColor='k';
% colormap(ax4,redblue)
% ylabel('Taxonomic Richness')
% colorbar
% caxis([0 80])
% ax4.Position=ax4.Position+[-0.05 0.06 -0.02 0];

set(gcf,'Color','none');

fname=['figures/Latitudinal_diversity_' num2str(iyr) '.png'];
export_fig(fname,'-r200')


%%

figure(11)
% Make symmetric (get rid of tiny numerical error)
D2=mean(cat(3,D2,D2'),3);
tree = linkage(squareform(D2),'average');
leafOrder = optimalleaforder(tree,D2);

dendrogram(tree,0,'ColorThreshold',0.5,'Orientation','left','ReOrder',leafOrder)




set(gcf,'Color','w');

%%

% Calculate Bray-Curtis at K sites
siteX=full(abundance(K,:));
clear minX
for i=1:nxc % nxc = n phenotypic species
    % for each species, calculate minimum between all pairs of locations
    minX(:,:,i)=min(siteX(:,i),siteX(:,i)'); 
end
C=sum(minX,3);
% BC = sum of only the lesser counts for each species found in both sites
% We do not need to normalise by S_i+S_j because we are using normalised abundances
% S_i+S_j = 2, so BC=1-2.*C./(S_i+S_j) = 1-C;
BCtaxo=1-C;
BCtaxo=BCtaxo.*(1-eye(size(BCtaxo))); % ensure diagonals are exactly zero

tree = linkage(squareform(BCtaxo),'average');
leafOrder = optimalleaforder(tree,BCtaxo);
dendrogram(tree,0,'ColorThreshold',0.5,'Orientation','left','ReOrder',leafOrder)
nclust=10;
T = cluster(tree,'MaxClust',nclust);

figure(12)
clf
ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
scatterm(ylat(K), xlon(K), 50, T,'filled');
geoshow(ax1, land, 'FaceColor', [0.7 0.7 0.7]);
colormap(lhsdesign(nclust,3))

axis off


%%
figure(13)
clf
sh=scatter(Temperature(K,end),richness_phen(K),50,abs(ylat(K)),'filled');
sh.MarkerFaceAlpha=1;
sh.MarkerFaceAlpha=0.25;
sh.MarkerEdgeColor='k';
colormap(ax3,redblue)
xlabel('Temperature ^\circC')
ylabel('Richness')



set(gcf,'Color','none');



























return
%%


T   = squeeze(Temperature(isite,:));
x   = squeeze(x(isite,:,:));
clr = squeeze(rgb(isite,:,:,:));
clr=permute(clr,[1 3 2]);

clr=(1+clr./cscale)./2;

subplot(211)
imagesc(1:nday,T_opt,x)
hold on
plot(1:nday,T,'k','LineWidth',3)
axis xy

subplot(212)
imagesc(1:nday,T_opt,clr)
hold on
contour(1:nday,T_opt,x,10.^[-3:0],'k','LineWidth',0.75)
plot(1:nday,T,'k','LineWidth',3)
axis xy



%%
while true
for i=1:nday
    data=x(:,i);
    
    clf
    h=bar(T_opt,data);
    
    rgb=squeeze(clr(:,i,:));
    h.FaceColor = 'flat';
    h.CData=rgb;
    
    hold on
    plot(T(i),0,'ko')
    
    axis([-2 36 1e-15 1])
    title(i)
    set(gca,'YScale','Log')
    drawnow
end
end
%%






























return
%%
if plotmaps && rem(dy,nday)==0
    toc
    c=c+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot metrics
    figure(fh3)
    clf
    
    subplot(221)
    ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',lonrng);
    ax1.Colormap=wntr;
    
    clr=log10(sum(X,2));
    
    cx=[0 23];
    cx=linspace(cx(1),cx(2),size(wntr,1));
    clr(find(clr>0 & clr<cx(2)))=cx(2); % Move smallest non zero abundances to second color bin
    
    fld=vector2gcmfaces(clr,iface,ix,iy,iz);
    [lon lat clr]=convert2pcol(mygrid.XC,mygrid.YC,fld);
    surfm(lat,lon,clr);
    colorbar
    caxis([cx(1) cx(end)])
    %                 shading flat;
    %                 geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    axis off;
    
    subplot(222)
    ax1=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',lonrng);
    ax1.Colormap=wntr;
    
    clr=log10(sum(x,2));
    
    cx=[-25 0];
    cx=linspace(cx(1),cx(2),size(wntr,1));
    clr(find(clr>0 & clr<cx(2)))=cx(2); % Move smallest non zero abundances to second color bin
    
    fld=vector2gcmfaces(clr,iface,ix,iy,iz);
    [lon lat clr]=convert2pcol(mygrid.XC,mygrid.YC,fld);
    surfm(lat,lon,clr);
    colorbar
    caxis([cx(1) cx(end)])
    %                 shading flat;
    %                 geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    axis off;
    
    subplot(223)
    ax2=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',lonrng);
    ax2.Colormap=parula;
    tmp=x.*log(x);
    tmp(~isfinite(tmp))=0;
    shannon = -sum(tmp,2);
    
    fld=vector2gcmfaces(shannon,iface,ix,iy,iz);
    [lon lat clr]=convert2pcol(mygrid.XC,mygrid.YC,fld);
    surfm(lat,lon,clr);
    colorbar
    %             caxis([cx(1) cx(end)])
    %                 shading flat;
    %                 geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    axis off;
    
    
    title(['Year ' sprintf('%3d', yr) '; Day ' sprintf('%3d', dy)],'FontSize',20);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot allele frequency
    figure(fh1)
    clf
    colormap(wntr)
    for i=1:nxc
        subplot(9,9,i);
        ax=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',lonrng);
        ax.Position(3:4)=ax.Position(3:4).*1.25;
        
        clr=(x(:,i));
        
        cx=[0 1];
        cx=linspace(cx(1),cx(2),size(wntr,1));
        clr(find(clr>0 & clr<cx(2)))=cx(2); % Move smallest non zero abundances to second color bin
        
        fld=vector2gcmfaces(clr,iface,ix,iy,iz);
        [lon lat clr]=convert2pcol(mygrid.XC,mygrid.YC,fld);
        surfm(lat,lon,clr);
        
        
        %                 shading flat;
        %                 geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
        axis off;
        caxis([0 1]);
        title(['T_{opt} = ' num2str(T_opt(i)) '^{\circ}C'],'FontSize',14);
    end
    text(-24,15,['Year ' sprintf('%3d', yr) '; Day ' sprintf('%3d', dy)],'FontSize',20)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot RGB
    figure(fh2)
    colormap(cmap);
    
    tmp=reshape(full(rgb),nxr,nxc,[]);      % reshape to [nxr x nxc x 3]
    cscale=max(cscale,max(abs(tmp(:))));    % normalise relative to largest GLOBAL rgb value (throughout run)
    tmp=(1+tmp./cscale)./2;                 % normalise
    tmp=max(1,ceil(tmp.*16));               % snap to color space grid
    for i=1:nxc
        subplot(9,9,i);
        ax=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',lonrng);
        ax.Position(3:4)=ax.Position(3:4).*1.25;
        
        clr=squeeze(tmp(:,i,:));
        clr(X(:,i)==0,:)=8;
        clr=sub2ind([nclr nclr nclr],clr(:,1),clr(:,2),clr(:,3));
        
        fld=vector2gcmfaces(clr,iface,ix,iy,iz);
        [lon lat clr]=convert2pcol(mygrid.XC,mygrid.YC,fld);
        surfm(lat,lon,clr);
        %                 shading flat;
        %                 geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
        axis off;
        caxis([1 16^3]);
        title(['T_{opt} = ' num2str(T_opt(i)) '^{\circ}C'],'FontSize',14);
    end
    toc
    text(-24,15,['Year ' sprintf('%3d', yr) '; Day ' sprintf('%3d', dy)],'FontSize',20)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fh1.InvertHardcopy = 'off';
    fh2.InvertHardcopy = 'off';
    fh3.InvertHardcopy = 'off';
    fname1=['figures/Frequency_' num2str(c,'%05i') '.png'];
    fname2=['figures/RGB_genes_' num2str(c,'%05i') '.png'];
    fname3=['figures/Invasion_' num2str(c,'%05i') '.png'];
    %             tic;print(fname, '-dpng','-r300');toc
    tic;saveas(fh1,fname1);toc
    tic;saveas(fh2,fname2);toc
    tic;saveas(fh3,fname3);toc
    disp('--------------------------------')
    tic
end % end plotmaps

if plotseries && rem(dy,28)==0
    figure(10)
    clf
    
    subplot(211)
    imagesc(1:di,T_opt,xday(:,1:di))
    hold on
    plot(1:di,theta_tseries(1:di),'k-','LineWidth',2)
    axis xy
    %             caxis([0 1])
    colorbar('SouthOutside')
    
    subplot(212)
    clr=rgbday(:,1:di,:);
    clrscl=max(abs(clr(:)));
    clr=(1+clr./clrscl)./2;
    imagesc(1:di,T_opt,clr)
    hold on
    plot(1:di,theta_tseries(1:di),'k-','LineWidth',2)
    axis xy
    drawnow
    
    
    figure(99)
    di=(yr-1).*nday+dy;
    clrscl=max(abs(rgbday(:)));
    
    for i=1:di
        data=xday(:,i);
        h=bar(T_opt,data);
        
        clr=squeeze(rgbday(:,i,:));
        clr=(1+clr./clrscl)./2;
        h.FaceColor = 'flat';
        h.CData=clr;
        
        hold on
        plot(T_opt,s(isite,:).*data','k','LineWidth',2)
        
        ylim([0 1])
        title(i)
        drawnow
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%