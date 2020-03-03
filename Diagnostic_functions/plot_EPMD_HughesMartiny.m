clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;
load('Basin8.mat');

input_filename = {...
    'neutral_stochastic_static_GUD_X01_weighted_transport',...
    'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport',...
    'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01',...
    'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1'};

%%


grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

%%
disp('::::::::::::::::::::::::::::::::::')
for i=1:4
    pathname   = '~/GitHub/EPMD/Output/';
    matObj  = matfile([pathname input_filename{i} '.mat']);
    
    clear  t_occupied
    ocean{i}    = matObj.ocean;
    t_occupied  = matObj.t_occupied;
    
    i_lastyr{i} = matObj.yrs_saved;
    disp([num2str(i_lastyr{i}) ' years evaluated.'])
        
    T = ocean{i}.forcing_temp;
    t_occ{i}=full(t_occupied(:,:));
    
    ctime=t_occ{i}(:);
    ctime(isnan(ctime))=inf;
    
    disp(input_filename{i})
    disp(['ocean is ' num2str(100.*nnz(t_occ{i})/numel(t_occ{i}),'%2.0f') '% connected.'])
    disp('::::::::::::::::::::::::::::::::::')
    
    t_occ{i}(~t_occ{i}(:)) = NaN;
    
    K{i}  = mean(ocean{i}.forcing_PCapacity,2);
    xx{i} = cell2mat(matObj.x(i_lastyr{i},1));
    
    run_options{i} = matObj.run_options;
end

%%
pathname   = '~/GitHub/EPMD_description/';
load Tara_sites.mat

insurf = find(Tara.depth=='SRF');

loaddata = load('StatID_sz2.mat');
StatID_sz2 = loaddata.StatID_sz2;

loaddata = load('~/GitHub/EPMD_description/EPMD_Figures/Tara_Basin_sz2.mat');
Tara_Basin = loaddata.Tara_Basin;


[~,ia,~] = intersect(Tara.StationID,StatID_sz2);

Tlat = Tara.lat(ia);
Tlon = Tara.lon(ia);
Tlon(Tlon<0) = Tlon(Tlon<0)+360;

EucDist = sqrt( (Tlat - ocean{1}.lat').^2 + (Tlon - ocean{1}.lon').^2 );
[~,Tara_ind] = min(EucDist,[],2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cluster diagrams
fG=figure(1);
% f2.Position = [209 1 930 1344];
clf
axlim=50;
fntsz=16;
ttls={'Neutral','Selective','Adaptive 1%','Adaptive 10%'};

%% Richter 'Geographic distance'
load('~/GitHub/EPMD/Richter_2020/supplementary_table_16.geographic_distance.mat');
Geographic_data      = Data;
Geographic_distance  = Geographic_data{:,2:end};
Geographic_distance(isnan(Geographic_distance))=0;
Geographic_stations  = char(Geographic_data.VarName1);
Geographic_stationID = str2num(Geographic_stations(:,1:3));

[~,ia,~] = intersect(StatID_sz2,Geographic_stationID);
[~,ib,~] = intersect(StatID_sz2,Tara.StationID(insurf));

Gdist=Geographic_distance(ia,ia);
Gdist=(Gdist+Gdist')./2;

G1D=cmdscale(Gdist,3);
G1D(G1D<prctile(G1D(:),05))=prctile(G1D(:),05);
G1D(G1D>prctile(G1D(:),95))=prctile(G1D(:),95);
G1D=normalize(G1D,'range');
% Environmental distance = temperature difference
Edist=sqrt((ocean{1}.ann_theta(Tara_ind)-ocean{1}.ann_theta(Tara_ind)').^2);
E1D=cmdscale(Edist,3);
E1D=normalize(E1D,'range');
%%

%% t-SNE analysis

disp('----------------------------------')
for i=1:4
    x=xx{i};    
    if i==2
        % selective experiment has lots of empty columns...
        x=reshape(full(x),60646,77,[]);
        x=sparse(squeeze(sum(x,2)));
    end
    x=x(:,1:94);   % remove global resident
    xxx{i}=x;      % save full frequency array
    x=x./sum(x,2); % normalise so sum x is 1 at each site
    S{i}=full(x(Tara_ind,:)); % extract Tara sites
    S{i}=S{i}./sum(S{i},2); % normalise so sum x is 1 at each site
    
    for isites=1:size(S{i},1)
        Mdist{i}(:,isites) = bray_curtis(S{i}(isites,:),S{i});
    end

    opts = statset('MaxIter',1e5);
    perplexity=20;
    disp(['Analysing ' input_filename{i}])
    disp('Performing 2D t-SNE analysis')
    C2D{i}  = tsne(S{i},'Algorithm','exact',...
                        'NumDimensions',2,...
                        'Distance',@bray_curtis,...
                        'Perplexity',perplexity,...
                        'Options',opts);
    disp('Performing 3D t-SNE analysis')
    C3D{i}  = tsne(S{i},'Algorithm','exact',...
                        'NumDimensions',3,...
                        'Distance',@bray_curtis,...
                        'Perplexity',perplexity,...
                        'Options',opts);
    disp('----------------------------------')
    
%     C2D{i}=normalize(C2D{i},'range');
%     C3D{i}=normalize(C3D{i},'range');
    C2D{i} = (C2D{i} - min(C2D{i}(:))) ./ (max(C2D{i}(:)) - min(C2D{i}(:)));
    C3D{i} = (C3D{i} - min(C3D{i}(:))) ./ (max(C3D{i}(:)) - min(C3D{i}(:)));
end

%% 
f1=figure(1);
f1.Position=[375 108 1128 1237];
clf

for i=1:4
    
    dotsz=75;
    fsize=15;
    
    sh=subplot(4,4,i*4-3);
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,grp2idx(Tara_Basin),'filled')
    set(gca,'XTick',[],'YTick',[]),axis square
    shpos=sh.Position;
    if i==1
        ch=colorbar;
        ch.Location='westoutside'
        ch.Ticks=linspace(1.4375,7.5625,8);
        ch.TickLabels=cellstr(unique(Tara_Basin));
        ch.FontSize=fsize;
    end
    caxis([1 8])
    hold on
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,'k')
    title({[ttls{i}],['(' char(96+4*(i-1)+1) ') Ocean Basins']},'FontSize',fsize)
    box on
    axis([-0.05 1.05 -0.05 1.05])
    sh.Colormap=bsnclr;
    sh.Position=shpos;
    sh.Position(2)=sh.Position(2)+0.03.*(i-1);
    sh.Position(1)=sh.Position(1)+0.03;
    
    sh=subplot(4,4,i*4-2);
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,E1D,'filled')
    set(gca,'XTick',[],'YTick',[]),axis square
    hold on
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,'k')
    sh.Colormap=redblue;
    rho   = corr(Mdist{i}(:),Edist(:),'Type','Spearman');
    rho_p = partialcorr(Mdist{i}(:),Edist(:),Gdist(:),'Type','Spearman');
    title({['(' char(96+4*(i-1)+2) ') Environmental distance']},'FontSize',fsize)
    box on
    axis([-0.05 1.05 -0.05 1.05])
    sh.Position(2)=sh.Position(2)+0.03.*(i-1);
    sh.Position(1)=sh.Position(1);
    
    sh=subplot(4,4,i*4-1);
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,G1D,'filled')
    set(gca,'XTick',[],'YTick',[]),axis square
    hold on
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,'k')
%     text(C2D{i}(:,1),C2D{i}(:,2),{basinname{Tara_ind}});
    rho   = corr(Mdist{i}(:),Gdist(:),'Type','Spearman');
    rho_p = partialcorr(Mdist{i}(:),Gdist(:),Edist(:),'Type','Spearman');
    title({['(' char(96+4*(i-1)+3) ') Geographic distance']},'FontSize',fsize)
    box on
    axis([-0.05 1.05 -0.05 1.05])
    sh.Position(2)=sh.Position(2)+0.03.*(i-1);
    sh.Position(1)=sh.Position(1)-0.03;
    
    sh=subplot(4,4,i*4);
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,C3D{i},'filled')
    set(gca,'XTick',[],'YTick',[]),axis square
    hold on
    scatter(C2D{i}(:,1),C2D{i}(:,2),dotsz,'k')
%     text(C2D{i}(:,1),C2D{i}(:,2),{basinname{Tara_ind}});
    rho   = corr(Mdist{i}(:),Gdist(:),'Type','Spearman');
    rho_p = partialcorr(Mdist{i}(:),Gdist(:),Edist(:),'Type','Spearman');
    title({['(' char(96+4*(i-1)+4) ') Community dissimilarity']},'FontSize',fsize)
    box on
    axis([-0.05 1.05 -0.05 1.05])
    sh.Position(2)=sh.Position(2)+0.03.*(i-1);
    sh.Position(1)=sh.Position(1)-0.06;
    
    set(gcf,'color','w')
    drawnow
end





sname=[pathname '/EPMD_Figures/EPMD_Full_cluster_maps.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')





%% Maps
f2=figure(2);
f2.Position=[1253 8 1307 1337];
clf
dotsz=75;
fontsz=18;

smpl_pnts = ocean{1}.sample_points;
n_eg=8;
clear shpos
for i=1:4
    
    sh=subplot(4,3,i*3-2);
    X=xxx{i}(:,n_eg).*ocean{i}.ann_abundance;
    [ax] = plot_vector(X,'lin',mygrid,ocean{1});
    scatterm(ocean{i}.lat(smpl_pnts(n_eg)),ocean{i}.lon(smpl_pnts(n_eg)),dotsz,'m','filled')
    scatterm(ocean{i}.lat(smpl_pnts(n_eg)),ocean{i}.lon(smpl_pnts(n_eg)),dotsz,'k')
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
    th=title(['(' char(96+i*3-2) ')'],'FontSize',fontsz);
    th.Position(2)=th.Position(2)+0.09;
    axis off
    drawnow
    sh.Position(1)=sh.Position(1)-0.05;
    sh.Position(2)=sh.Position(2)+(i-1)*0.065;
    sh.Position(3)=sh.Position(3).*1.2;
    shpos{i}=sh.Position;
    sh.Colormap=turbo;
    ch=colorbar('Location','westoutside');
    ch.Position(1) = ch.Position(1) - 0.03;
    ch.Position(2) = ch.Position(2) + ch.Position(4) .* 0.1;
    ch.Position(4) = ch.Position(4) .* 0.8;
    ch.FontSize=fontsz;
    caxis([0 max(caxis)]);
    sh.Position=shpos{i};
    
    sh=subplot(4,3,i*3-1);
    ax=axesm('MapProjection','mollweid','frame','on','MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
    scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,C3D{i},'filled')
    scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,'k')
    th=title(['(' char(96+i*3-1) ')'],'FontSize',fontsz);
    th.Position(2)=th.Position(2)+0.09;
    axis off
    sh.Position=shpos{i};
    sh.Position(1)=shpos{i}(1)+0.275;
    
    sh=subplot(4,3,i*3);
    scatter(C2D{i}(:,1),C2D{i}(:,2),50,grp2idx(Tara_Basin),'filled')
    hold on
    scatter(C2D{i}(:,1),C2D{i}(:,2),50,'k')
    set(gca,'XTick',[],'YTick',[]),axis square
    caxis([1 8])
    th=title(['(' char(96+i*3) ')'],'FontSize',fontsz);
    box on
    axis([-0.05 1.05 -0.05 1.05])
    if i==1
        ch=colorbar;
        ch.Location='eastoutside';
        ch.Ticks=linspace(1.4375,7.5625,8);
        ch.TickLabels=cellstr(unique(Tara_Basin));
        ch.FontSize=fontsz;
    end
    sh.Colormap=bsnclr;
%     sh.Position=shpos;
    sh.Position(4)=shpos{i}(4).*0.7;
    sh.Position(1)=shpos{i}(1)+0.54;
    sh.Position(2)=shpos{i}(2)+0.0325;
end

sname=[pathname '/EPMD_Figures/EPMD_cluster_maps.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

return
%%

figure(22)

subplot(311)
clrs = bsnclr(grp2idx(Tara_Basin),:);
ax=axesm('MapProjection','mollweid','frame','on','MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,clrs,'filled')
scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,'k')
title({'Ocean Basins'},'FontSize',20)
axis off


subplot(312)
ax=axesm('MapProjection','mollweid','frame','on','MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,G1D,'filled')
scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,'k')
title({'Geographic clustering'},'FontSize',20)
axis off

sh=subplot(313);
sh.Colormap=redblue;
ax=axesm('MapProjection','mollweid','frame','on','MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,E1D,'filled')
scatterm(ocean{i}.lat(Tara_ind),ocean{i}.lon(Tara_ind),dotsz,'k')
title({'Environmental clustering'},'FontSize',20)
axis off

set(gcf,'color','w')


%%

%% plot dendrogram and metadata
figure(3)
clf
for i=1:4
    % plot dendrogram
    subplot(4,3,(i-1)*3+[1 2])
    tree = linkage(squareform(Mdist{i}),'Ward');
    leafOrder = optimalleaforder(tree,Mdist{i});
    [H,T,outperm] = dendrogram(tree,0,...
        'Reorder',leafOrder,...
        'Orientation','left');
    
    set(gca,'YTickLabel',cellstr(Tara_Basin(outperm)),...
            'FontSize',5)
    % get the current tick labeks
    ticklabels = get(gca,'YTickLabel');
    % prepend a color for each tick label
    ticklabels_new = cell(size(ticklabels));
    for ii = 1:length(ticklabels)
        jj=outperm(ii);
        r=bsnclr(grp2idx(Tara_Basin(jj)),1);
        g=bsnclr(grp2idx(Tara_Basin(jj)),2);
        b=bsnclr(grp2idx(Tara_Basin(jj)),3);
        ticklabels_new{ii} = ['\color[rgb]{' num2str(r) ',' num2str(g) ',' num2str(b) '} ' ticklabels{ii}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_new);
        
    for ii=1:size(H)
        H(ii).Color='k';
        H(ii).LineWidth=1;
    end
    title({'Community Distance Clusters',ttls{i}},'FontSize',10);
    box on
    
    % plot environmental variables
    ncol=9;
    
    sh=subplot(4,ncol,(i-1)*ncol+ncol-2);
    clear clr
    clr(1,:,:)=C3D{i}(outperm,:);
    clr=permute(clr,[2 1 3]);
    imagesc(clr)
    axis xy
    ylim([0 size(Tara_ind,1)+1])
    if i==4
        set(gca,'XTick',1,...
            'XTickLabel',{'Community'})
        xtickangle(90)
    else
        set(gca,'XTick',[])
    end
    set(gca,'YTick',[])
    
    sh=subplot(4,ncol,(i-1)*ncol+ncol-1);
    clear clr
    clr(1,:,:)=G1D(outperm,:);
    clr=permute(clr,[2 1 3]);
    imagesc(clr)
    axis xy
    ylim([0 size(Tara_ind,1)+1])
    if i==4
        set(gca,'XTick',1,...
            'XTickLabel',{'Geographic'})
        xtickangle(90)
    else
        set(gca,'XTick',[])
    end
    set(gca,'YTick',[])
    
    sh=subplot(4,ncol,(i-1)*ncol+ncol);
    imagesc(ocean{1}.ann_theta(Tara_ind(outperm)))
    sh.Colormap=flipud(redblue(64));
    axis xy
    ylim([0 size(Tara_ind,1)+1])
    if i==4
        set(gca,'XTick',1,...
            'XTickLabel',{'Temperature'})
        xtickangle(90)
    else
        set(gca,'XTick',[])
    end
    set(gca,'YTick',[])
end
sname=[pathname 'Figures/Dendrogram_' Biodiversity_type '_' num2str(Size_class) '.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

%%
figure(4)
clf
for i=1:4
    % plot correlations
    clear fG xiG fE xiE
    
    [fG,xiG,bwG{i}]=ksdensity([Gdist(:),Mdist{i}(:)],[Gdist(:),Mdist{i}(:)],'PlotFcn','contour','Bandwidth',[0.5 0.05]);   
    
    subplot(4,2,i*2-1)
    scatter(xiG(:,1),xiG(:,2),10,fG,'filled')
    rho   = corr(Mdist{i}(:),Gdist(:),'Type','Spearman');
    rho_p = partialcorr(Mdist{i}(:),Gdist(:),Edist(:),'Type','Spearman');
    title({ttls{i},['Rank Correlation = ' num2str(rho)],['Partial Rank Correlation = ' num2str(rho_p)]})
    xlabel({'Geographic distance'},'FontSize',20)
    ylabel({'Community distance'},'FontSize',20)
    
    
    [fE,xiE,bwE{i}]=ksdensity([Edist(:),Mdist{i}(:)],[Edist(:),Mdist{i}(:)],'PlotFcn','contour','Bandwidth',[0.5 0.05]);
    
    subplot(4,2,i*2)
    scatter(xiE(:,1),xiE(:,2),10,fE,'filled')
    rho   = corr(Mdist{i}(:),Edist(:),'Type','Spearman');
    rho_p = partialcorr(Mdist{i}(:),Edist(:),Gdist(:),'Type','Spearman');
    title({ttls{i},['Rank Correlation = ' num2str(rho)],['Partial Rank Correlation = ' num2str(rho_p)]})
    xlabel({'Environmental distance'},'FontSize',20)
    ylabel({'Community distance'},'FontSize',20)
    
end


%%
clc
loaddata = load('Mdist_sz2.mat');
Mdist{5} = loaddata.Mdist;

figure(5)
for i=1:4
    rhoP = corr(Mdist{5}(:),Mdist{i}(:),'Type','Pearson');
    rhoS = corr(Mdist{5}(:),Mdist{i}(:),'Type','Spearman');
    disp([i rhoP rhoS])
    
    subplot(4,1,i)
    [f,xi,bw{i}]=ksdensity([Mdist{5}(:),Mdist{i}(:)],[Mdist{5}(:),Mdist{i}(:)],'PlotFcn','contour','Bandwidth',[0.5 0.05]);   
    scatter(xi(:,1),xi(:,2),10,f,'filled')
%     title({ttls{i},['Rank Correlation = ' num2str(rho)],['Partial Rank Correlation = ' num2str(rho_p)]})
    xlabel({'Tara distance'})
    ylabel({'EPMD distance'})
    
end



return
%% % Full map clustering - very slow
clc
opts = statset('MaxIter',10000);

figure(55)
clf
clear xp yp in_smpl sind1 xlon ylat rgb xdata

% initial analysis with fewer sample points 
in_smpl=randperm(60646,1000);

xlon = ocean{1}.lon(in_smpl);
ylat = ocean{1}.lat(in_smpl);

i=2

xdata = full(xxx{i}(in_smpl,1:94));  
    
rgb = tsne(xdata,...
    'Algorithm','barneshut',... % 'exact' or 'barneshut'
    'Distance', @bray_curtis,...
    'NumDimensions',3,...
    'Options',opts,...
    'Verbose',1,'NumPrint',1000);

rgb = normalize(rgb,'range');

%%

[X,Y,Z] = sph2cart(deg2rad(ocean{1}.lon),deg2rad(ocean{1}.lat),1);

clear F
for ic=1:3
    disp(ic)
    F{ic} = scatteredInterpolant(X(in_smpl),Y(in_smpl),Z(in_smpl),rgb(:,ic),'natural');
end
clear Vq
for ic=1:3
    disp(ic)
    Vq(:,ic) = F{ic}(X,Y,Z);
end

% subplot(312)
[ax] = plot_image(Vq,mygrid,ocean{1})
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
axis off
drawnow

sname=['Global_rgb1.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

return
%% Repeat at higher sampling using last iteration as initial condition

in_smpl=randperm(60646,5000);

xlon = ocean{1}.lon(in_smpl);
ylat = ocean{1}.lat(in_smpl)

rgb2 = tsne(full(xxx{i}(in_smpl,1:94)),...
        'Algorithm','barneshut',... % 'exact' or 'barneshut'
        'Distance', @bray_curtis,...
        'NumDimensions',3,...
        'Options',opts,...
        'InitialY',Vq(in_smpl,:),...
        'Verbose',1,'NumPrint',100);
    
rgb2 = normalize(rgb2,'range');


ax=axesm ('eqaazim','frame','on','FlineWidth',0.5,'Origin',[-46.283333,-86.666665]);
scatterm(ocean{1}.lat(in_smpl),ocean{1}.lon(in_smpl),10,rgb2,'filled');
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
axis off
drawnow



[X,Y,Z] = sph2cart(deg2rad(ocean{1}.lon),deg2rad(ocean{1}.lat),1);

clear F2
for ic=1:3
    disp(ic)
    F2{ic} = scatteredInterpolant(X(in_smpl),Y(in_smpl),Z(in_smpl),rgb2(:,ic),'natural');
end
clear Vq2
for ic=1:3
    disp(ic)
    Vq2(:,ic) = F2{ic}(X,Y,Z);
end

% subplot(312)
clf
[ax] = plot_image(Vq2,mygrid,ocean{1})
[x,y,u,v] = get_circ_vectors(ocean{i}); % second input is coarse graining resolution
iplot=randsample(60646,1e4);
h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
axis off
set(gcf,'Color','w')
drawnow

Vq=Vq2;

sname=['Global_rgb3.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')


return










