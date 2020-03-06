clear
clc
addpath ~/GitHub/EPMD/Diagnostic_functions/
load('Basin8.mat');

save_figs = 'false';

load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);



Biodiversity_type = 'Metagenomic'; % 'Metagenomic or OTU
for Size_class        = 2 % 1 to 6
    disp(Size_class)
    
    fnames = dir('~/GitHub/EPMD/Richter_2020/*dissimilarity*');
    
    if startsWith(Biodiversity_type,'Metagenomic')
        fnames = dir('~/GitHub/EPMD/Richter_2020/*metagenomic_dissimilarity*');
    elseif startsWith(Biodiversity_type,'OTU')
        fnames = dir('~/GitHub/EPMD/Richter_2020/*OTU_dissimilarity*');
    else
        error('Biodiversity_type must be ''Metagenomic'' or ''OTU''.')
    end
    
    Community_filename  = fnames(Size_class).name;
    Env_filename        = 'supplementary_table_02.environmental_parameters.mat';
    Geographic_filename = 'supplementary_table_16.geographic_distance.mat';
    
    pathname='~/GitHub/EPMD_description/';
    addpath(pathname)
    %% Environmental data
    load(Env_filename);
    Env_data      = Data;
    Env_stations  = char(Env_data.Station);
    Env_stationID = str2num(Env_stations(:,6:end));
    Env_layer     = Env_data.Depth;
    Env_Temp      = normalize(Env_data.TemperatureC);
    Env_Nitr      = normalize(Env_data.NO2NO3molL);
    Env_Phos      = normalize(Env_data.PhosphatemolL);
    Env_Iron      = normalize(Env_data.Ironmmolm3);
    Env_Sili      = normalize(Env_data.SilicateWOA13molL);
    Env_Chlo      = normalize(Env_data.Chlorophyllmgm3);
    OceanBasin    = Env_data.Oceanandsearegion;
    Biome         = Env_data.Marinebiome;
    Latitude      = Env_data.Latitude;
    Longitude     = Env_data.Longitude;
    
    EnvV=[Env_Temp Env_Nitr Env_Iron Env_Sili Env_Chlo];
    EnvV(isnan(EnvV))=0;
    
    % isSO=find(startsWith(cellstr(OceanBasin),'Southern Ocean'));
    % EnvV(isSO,:)=0;
    
    %% Community data
    load(Community_filename);
    Community_data      = Data;
    Community_distance  = Community_data{:,2:end};
    Community_stations  = char(Community_data.VarName1);
    Community_stationID = str2num(Community_stations(:,1:3));
    Community_layer     = cellstr(Community_stations(:,5:end));
    
    %% Geographic data
    load(Geographic_filename);
    Geographic_data      = Data;
    Geographic_distance  = Geographic_data{:,2:end};
    Geographic_distance(isnan(Geographic_distance))=0;
    Geographic_stations  = char(Geographic_data.VarName1);
    Geographic_stationID = str2num(Geographic_stations(:,1:3));
    insurface            = find(startsWith(Community_layer,'SRF'));
    
    %% Collate distance matrices
    % find intersect of environmental and community data
    [C,ia,ib] = intersect(Env_stationID,Community_stationID(insurface));
    
    OBind = grp2idx(OceanBasin(ia));
    EnvVars=EnvV(ia,:);
    
    Edist = squareform(pdist(EnvVars));
    Mdist = Community_distance(ib,ib);
    Gdist = Geographic_distance(C,C);
        
    % classical MDS
    [Geo eigen1] = cmdscale(Gdist,3);
    [Env eigen2] = cmdscale(Edist,3);
    [Com eigen3] = cmdscale(Mdist,3);
    
    opts = statset('MaxIter',1e5);
    Com  = tsne(Mdist,'Algorithm','exact',...
                      'NumDimensions',2,...
                      'Distance',@dist2dist,...
                      'Perplexity',20,...
                      'Options',opts);
    Com3 = tsne(Mdist,'Algorithm','exact',...
                      'NumDimensions',3,...
                      'Distance',@dist2dist,...
                      'Perplexity',20,...
                      'Options',opts);
%     Com3=Com;
%     Com=Com(:,1:2);
        
    % normalise distances
    Geo  = normalize(Geo,'range');
    Env  = normalize(Env,'range');
%     Com  = normalize(Com,'range');
%     Com3 = normalize(Com3,'range');
    Com = (Com - min(Com(:))) ./ (max(Com(:)) - min(Com(:)));
    Com3 = (Com3 - min(Com3(:))) ./ (max(Com3(:)) - min(Com3(:)));
    
    
    %% plot Martiny clusters
    fsize=20;
    
    f1=figure(1);
    f1.Position = [73 536 913 809];
    clf
    
    dotsz=100;
    
    sh1=subplot(221);
    scatter(Com(:,1),Com(:,2),dotsz,grp2idx(OceanBasin(ia)),'filled')
    hold on
    scatter(Com(:,1),Com(:,2),dotsz,'k')
    sh1.Position(1)=sh1.Position(1)+0.1;
    shpos=sh1.Position;
    axis square,box on
    title(['(a) Ocean basin'],'FontSize',fsize)
    axis([-0.05 1.05 -0.05 1.05])
    set(gca,'XTick',[],'YTick',[])
    sh1.Colormap=bsnclr;
    ch=colorbar;
    ch.Location='westoutside';
    ch.Ticks=linspace(1.4375,7.5625,8);
    ch.TickLabels=cellstr(unique(OceanBasin));
    ch.FontSize=14;
    caxis([1 8])
    sh1.Position=shpos;
    
    sh2=subplot(222);
    scatter(Com(:,1),Com(:,2),dotsz,Com3,'filled')
    hold on
    scatter(Com(:,1),Com(:,2),dotsz,'k')
    title(['(b) ' Biodiversity_type ' clusters'],'FontSize',fsize)
    axis square,box on
    axis([-0.05 1.05 -0.05 1.05])
    set(gca,'XTick',[],'YTick',[])
    sh2.Position(3:4)=shpos(3:4);
    sh2.Position(1)=sh2.Position(1)+0.1;
    sh2.Position(1)=sh2.Position(1)-0.1;
    
    sh3=subplot(223);
    scatter(Com(:,1),Com(:,2),dotsz,Env(:,1:3),'filled')
    hold on
    axis([-0.05 1.05 -0.05 1.05])
    scatter(Com(:,1),Com(:,2),dotsz,'k')
    rho   = corr(Edist(:),Mdist(:),'Type','Spearman');
    rho_p = partialcorr(Mdist(:),Edist(:),Gdist(:),'Type','Spearman');
    title(['(c) Environmental clusters'],'FontSize',fsize)
    axis square,box on
    set(gca,'XTick',[],'YTick',[])
    sh3.Position(1)=sh3.Position(1)+0.1;
    sh3.Position(3:4)=shpos(3:4);
    sh3.Position(2)=sh3.Position(2)+0.09;
    
    sh4=subplot(224);
    scatter(Com(:,1),Com(:,2),dotsz,Geo(:,1:3),'filled')
    hold on
    axis([-0.05 1.05 -0.05 1.05])
    scatter(Com(:,1),Com(:,2),dotsz,'k')
    rho   = corr(Gdist(:),Mdist(:),'Type','Spearman');
    rho_p = partialcorr(Mdist(:),Gdist(:),Edist(:),'Type','Spearman');
    title(['(d) Geographic distance clusters'],'FontSize',fsize)
    axis square,box on
    set(gca,'XTick',[],'YTick',[])
    sh4.Position(1)=sh4.Position(1)+0.1;
    sh4.Position(3:4)=shpos(3:4);
    sh4.Position(2)=sh4.Position(2)+0.09;
    sh4.Position(1)=sh4.Position(1)-0.1;
    
    if save_figs
        sname=[pathname 'Tara_Figures//Martiny_Clusters_' Biodiversity_type '_' num2str(Size_class) '.png'];
        set(gcf,'Color','w')
        export_fig(sname,'-r300')
    end
    %% plot dendrogram and metadata
    f2 = figure(2);
    f2.Position = [537 8 873 1337];
    clf
    
    % plot dendrogram
    subplot(1,5,[1 2])
    % test different clustering methods
    method={'single','complete','average','weighted','centroid','median','ward'};
    for imeth=1:7
        tree = linkage(squareform(Mdist),method{imeth});
        cnet(imeth)=cophenet(tree,squareform(Mdist));
    end
    [~,iord] = sort(cnet,'descend');
    disp(['UPGMA ranked #' num2str(find(iord==3)) '.'])
    disp(['(''' method{iord(1)} ''' ranked #1.)'])
    
    % then use UPGMA
    tree = linkage(squareform(Mdist),'average');
    leafOrder = optimalleaforder(tree,Mdist);
    [H,T,outperm] = dendrogram(tree,0,...
        ...'Reorder',leafOrder,...
        'Orientation','left');    
    title({'Community','Distance','Clusters','(a)'},'FontSize',fsize)
    set(gca,'XTick',[],'YTick',[])
    box on
    % get the current tick labels
    set(gca,'YTickLabel',cellstr(OceanBasin(ia(outperm))),'XTick',[])
    % get the current tick labeks
    ticklabels = get(gca,'YTickLabel');
        
    ncol=14;
    
    sh=subplot(1,ncol,ncol-2);
    clear clr
    clr(1,:,:)=Com3(outperm,:);
    clr=permute(clr,[2 1 3]);
    imagesc(clr)
    axis xy
    ylim([0 size(Com3,1)+1])
    set(gca,'XTick',1,...
        'XTickLabel',{'Community'},...
        'FontSize',fsize)
    xtickangle(90)
    set(gca,'YTick',[])
    shpos=sh.Position;
    sh.Position(1)=shpos(1)-0.32;
    
    sh=subplot(1,ncol,ncol-1);
    clear clr
    clr(1,:,:)=Geo(outperm,1:3);
    clr=permute(clr,[2 1 3]);
    imagesc(clr)
    axis xy
    ylim([0 size(Com,1)+1])
    set(gca,'XTick',1,...
        'XTickLabel',{'Geographic'},...
        'FontSize',fsize)
    xtickangle(90)
    set(gca,'YTick',[])
    title({'RGB','clusters',' ','(b)     (c)     (d)'},'FontSize',fsize)
    shpos=sh.Position;
    sh.Position(1)=shpos(1)-0.32;
    
    sh=subplot(1,ncol,ncol);
    clear clr
    clr(1,:,:)=Env(outperm,1:3);
    clr=permute(clr,[2 1 3]);
    imagesc(clr)
    axis xy
    ylim([0 size(Com,1)+1])
    set(gca,'XTick',1,...
        'XTickLabel',{'Environmental'},...
        'FontSize',fsize)
    xtickangle(90)
    set(gca,'YTick',[])
    shpos=sh.Position;
    sh.Position(1)=shpos(1)-0.32;
    
    
    % plot environmental variables
    sh=subplot(1,ncol,[ncol-4 ncol-2]);
    imagesc(EnvVars(outperm,:))
    sh.Colormap=flipud(redblue(64));
    axis xy
    ylim([0 size(Com,1)+1])
    set(gca,'XTick',1:size(EnvV,2),...
        'XTickLabel',{'Temperature','Nitrate','Iron','Silicate','Chlorophyll a'},...
        'FontSize',fsize)
    xtickangle(90)
    set(gca,'YTick',1:numel(ia))
    title({'Normalised','Environmental','Variables','(e)'},'FontSize',fsize)
    shpos=sh.Position;
    sh.Position(1)=shpos(1)-0.04;

    % prepend a color for each tick label
    ticklabels_new = cell(size(ticklabels));
    Tara_Basin = OceanBasin(ia);
    set(gca,'YTickLabel','','TickLength',[0 0])
    for ii = 1:length(ticklabels)
        jj=outperm(ii);
        r=bsnclr(grp2idx(Tara_Basin(jj)),1);
        g=bsnclr(grp2idx(Tara_Basin(jj)),2);
        b=bsnclr(grp2idx(Tara_Basin(jj)),3);
        ticklabels_new{ii} = ['\color[rgb]{' num2str(r) ',' num2str(g) ',' num2str(b) '} ' ticklabels{ii}];
        text(5.5,ii,ticklabels_new{ii},'FontSize',12)
    end
    % set the tick labels
%     set(gca, 'YTickLabel', ticklabels_new);
    for i=1:size(H)
        H(i).Color='k';
        H(i).LineWidth=1;
    end
    box on
    
    
    if save_figs
        sname=[pathname 'Tara_Figures//Dendrogram_' Biodiversity_type '_' num2str(Size_class) '.png'];
        set(gcf,'Color','w')
        export_fig(sname,'-r300')
    end
    %% plot maps
    
    dotsz=75;
    f3 = figure(3);
    f3.Position = [1000 875 860 400];
    clf
    
    ax=axesm('MapProjection','mollweid','frame','on');
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
    scatterm(Latitude(ia),Longitude(ia),dotsz,Com3,'filled')
    scatterm(Latitude(ia),Longitude(ia),dotsz,'k')
    title({'Community clustering'},'FontSize',fsize)
    axis off
    
    if save_figs
        sname=[pathname 'Tara_Figures//Map_' Biodiversity_type '_Dist_' num2str(Size_class) '.png'];
        set(gcf,'Color','w')
        export_fig(sname,'-r300')
    end
    
    f4 = figure(4);
    f4.Position = [1000 1 860 800];
    clf
    
    subplot(211)
    ax=axesm('MapProjection','mollweid','frame','on');
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
    scatterm(Latitude(ia),Longitude(ia),dotsz,Geo(:,1:3),'filled')
    scatterm(Latitude(ia),Longitude(ia),dotsz,'k')
    title('Geographic clustering','FontSize',fsize)
    axis off
    
    sh=subplot(212);
    ax=axesm('MapProjection','mollweid','frame','on');
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
    scatterm(Latitude(ia),Longitude(ia),dotsz,Env(:,1:3),'filled')
    scatterm(Latitude(ia),Longitude(ia),dotsz,'k')
    title('Environmental clustering','FontSize',fsize)
    axis off
    sh.Position(2)=sh.Position(2)+0.1;
    
    if save_figs
        sname=[pathname 'Tara_Figures//Maps_Env_Geo_' Biodiversity_type '_' num2str(Size_class) '.png'];
        set(gcf,'Color','w')
        export_fig(sname,'-r300')
    end
    
    %%
    %% Maps
f5=figure(5);
f5.Position=[1253 918 1307 427];
clf
dotsz=75;

n_eg=8;
clear shpos
    
sh=subplot(1,3,[1 2]);
ax=axesm('MapProjection','mollweid','frame','on');
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7])
scatterm(Latitude(ia),Longitude(ia),dotsz,Com3,'filled')
scatterm(Latitude(ia),Longitude(ia),dotsz,'k')
title({'(a)'},'FontSize',fsize)
axis off

sh=subplot(133);
scatter(Com(:,1),Com(:,2),150,grp2idx(OceanBasin(ia)),'filled')
hold on
scatter(Com(:,1),Com(:,2),150,'k')
shpos=sh.Position;
sh.Position(1)=shpos(1)-0.05;
shpos=sh.Position;
set(gca,'XTick',[],'YTick',[]),axis square
ch=colorbar;
ch.Location='eastoutside';
ch.Ticks=linspace(1.4375,7.5625,8);
ch.TickLabels=cellstr(unique(Tara_Basin));
    ch.FontSize=14;
caxis([1 8])
title(['(b)'],'FontSize',fsize)
box on
axis([-0.05 1.05 -0.05 1.05])
sh.Position=shpos;
sh.Colormap=bsnclr;

sname=[pathname 'Tara_Figures//Richter_' Biodiversity_type '_cluster_' num2str(Size_class) '.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
end