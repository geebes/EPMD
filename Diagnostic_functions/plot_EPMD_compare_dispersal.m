clear
clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;

input_filename = {   'neutral_stochastic_static_GUD_X01_surface_transport',...
                   'neutral_stochastic_seasonal_GUD_X01_surface_transport',...
                     'neutral_stochastic_static_GUD_X01_weighted_transport',...
                   'neutral_stochastic_seasonal_GUD_X01_weighted_transport',...
       'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport',...
     'nonadaptive_dispersal_stochastic_seasonal_GUD_X01_weighted_transport',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01',...
       'selective_dispersal_stochastic_seasonal_GUD_X01_weighted_transport_m0.01',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1',...
       'selective_dispersal_stochastic_seasonal_GUD_X01_weighted_transport_m0.1',...
...
                     'neutral_stochastic_static_GUD_X17_surface_transport',...
                   'neutral_stochastic_seasonal_GUD_X17_surface_transport',...
                     'neutral_stochastic_static_GUD_X17_weighted_transport',...
                   'neutral_stochastic_seasonal_GUD_X17_weighted_transport',... 
       'nonadaptive_dispersal_stochastic_static_GUD_X17_weighted_transport',...
     'nonadaptive_dispersal_stochastic_seasonal_GUD_X17_weighted_transport',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.01',...
       'selective_dispersal_stochastic_seasonal_GUD_X17_weighted_transport_m0.01',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.1',...
       'selective_dispersal_stochastic_seasonal_GUD_X17_weighted_transport_m0.1'}; 
  
lgnd = {'Neutral, surface',...
        'Neutral, depth-integrated',...
        'Selective, depth-integrated',...
        'Adaptive, depth-integrated, 1% mutation',...
        'Adaptive, depth-integrated, 10% mutation'};

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

for i=1:numel(input_filename)
    pathname   = '~/GitHub/EPMD/Output/';
    matObj  = matfile([pathname input_filename{i} '.mat']);

    clear ocean t_occupied
    ocean       = matObj.ocean;
    t_occupied  = matObj.t_occupied;

    i_lastyr{i} = matObj.yrs_saved;
    disp([num2str(i_lastyr{i}) ' years evaluated.'])
 

    T = ocean.forcing_temp;
    t_occ{i}=full(t_occupied(:,1:numel(ocean.sample_points)));
    
    ctime=t_occ{i}(:);
    ctime(isnan(ctime))=inf;
    
    disp(input_filename{i})
    disp(['ocean is ' num2str(100.*nnz(t_occ{i})/numel(t_occ{i}),'%2.0f') '% connected.'])
    disp(':::::::::::::::::::::::::::')

    t_occ{i}(~t_occ{i}(:)) = NaN;
    
    K{i}  = mean(ocean.forcing_PCapacity,2);
    xx{i} = cell2mat(matObj.x(i_lastyr{i},1));
    
    run_options{i} = matObj.run_options;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4 = figure(4);
f4.Position=[89 692 1160 653];
clf

set(gcf,'defaultAxesColorOrder',[023 063 095;
                                 032 099 154;
                                 060 173 162;
                                 245 212 092;
                                 236 085 059;
                                 000 000 000]./255)
xscale='log';

for i=1:numel(input_filename)
    if contains(input_filename{i},'X01')
        ax=subplot(2,1,1);
        t=title(['(a) 0.6 ' char(181) 'm \it{Prochlorococcus}'],'HorizontalAlignment','Left');
    else
        ax=subplot(2,1,2);
        t=title(['(b) 6 ' char(181) 'm Diatom'],'HorizontalAlignment','Left');        
    end
    t.Position(1)=7/365;
    
    hold on
    clear N edges
    [N,edges] = histcounts(t_occ{i}(:),0:(1/96):i_lastyr{i});
    
    
    y=[NaN cumsum(N)./numel(t_occ{i}) NaN];
    if contains(input_filename{i},'static')
        hp=plot([0 edges],y,'LineW',2);
    else
        ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
        hp=plot([0 edges],y,':','LineW',2,'HandleVisibility','off');
    end
    imed=find(y<=0.5,1,'last');
    disp(['Median time = ' num2str(hp.XData(imed)) ' years.'])
    
    ylim([0 1])
    box on
    ylabel('Connectance') 
    
    switch xscale
        case 'log'
            set(gca,'XScale','log')
            xlim([7/365 100])
            set(gca,'XTick',[7/365 1/12 1 10 100],'YTick',0:0.1:1);
            set(gca,'XTickLabels',{'week','month','year','decade','century'});
            grid on
            ax.XMinorTick = 'off';
            ax.XAxis.MinorTickValues = [1/365 (1:11)./12 1 2 3 4 5 6 7 8 9 10:10:100];
            ax.XMinorGrid = 'on';
        case 'lin'
            set(gca,'XScale','lin')
            xlim([0 30])
            set(gca,'XTick',0:30,'YTick',0:0.1:1);
            grid on           
    end
    
%     if contains(input_filename{i},'GUD_X17_surface_transport')
%         delete(hp)
%     end
    set(gca,'FontSize',14)
end
subplot(211);
switch xscale
    case 'log'
        legend(lgnd,'Location','NorthWest')
        sname=['../Figures/cumulative_connections.png'];
    case 'lin'
        legend(lgnd,'Location','East')
        sname=['../Figures/cumulative_connections_linear.png'];
end
subplot(212);
xlabel('Time')

set(gcf,'Color','w');
export_fig(sname,'-r300')

return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
set(gcf,'Position',[123 1 1879 1344])
clf

seed_ID=8;

% get abundance data

pp=0;
rp=0;
for i=[3 5 7 9]
    rp=rp+1;
    
    x = xx{i} .* K{i};
%     n = 100; % sample volume
%     x=floor(x .* n ./ ocean.volume); % Veil based on n litre sample size
    
    % sum across seed populations
    switch run_options{i}.seed_dist
        case 'nonadaptive_dispersal'
            x2 = reshape(full(x),60646,run_options{i}.nphen,run_options{i}.nlineages);
            x  = squeeze(sum(x2,2));
    end
    
    x_i=x(:,seed_ID);
    
    pp=pp+1;
    sh1=subplot(4,3,pp);
    [ax] = plot_vector(x_i,'lin',mygrid,ocean);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
%     caxis([0 1e20])
    hold on
    title(['(' char(96+pp) ')'])
    scatterm(ocean.lat(ocean.sample_points(seed_ID)),ocean.lon(ocean.sample_points(seed_ID)),25,'m')
    ch=colorbar;
    ch.Location = 'westoutside';
%     ch.TickLabels={'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'};
    sh1.Position(1) = sh1.Position(1) + 0.02;
    sh1.Position(2) = sh1.Position(2) + (rp-1).*0.04;
    set(gca,'FontSize',14)
   
    pp=pp+1;
    sh2=subplot(4,3,pp);    
    [n,edges,bin] = histcounts(ocean.ann_theta,linspace(-2,36,77));
    for ib=1:numel(edges)
        freq(ib)=sum(x_i(bin==ib));
    end
    bar(edges,freq,1)
    hold on
%     in9999=cumsum(freq)>sum(freq).*0.01 & cumsum(freq)<sum(freq).*0.99;
%     bar(edges,freq.*in9999,1)
    mask=edges==round(2*ocean.ann_theta(ocean.sample_points(seed_ID)))./2;
    bar(edges,freq.*mask,1)
    set(gca,'YScale','lin')
    if i==9
        xlabel('Temperature (^\circC)');
    else
        set(gca,'XTickLabel','')
    end    
    ylabel('Abundance')
    title(['(' char(96+pp) ')'])
%     legend('All','1% to 99%','T = T_{seed}','Location','NorthWest')
    xlim([-2.5 36.5])
%     set(gca,'YTick',10.^(0:5:25),...
%             'YTickLabel',{'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'})
    sh2.Position(2) = sh2.Position(2) + (rp-1).*0.04 +0.02;
    set(gca,'FontSize',14)
   
    pp=pp+1;
    sh3=subplot(4,3,pp);
    ab=sort(x(ocean.sample_points,1:94),2,'descend');
    plot(1:94,ab','Color',[0 0 0 0.25],'LineW',1);
    hold on
%     plot(1:94,ab(seed_ID,:)','r-','LineW',1);
    axis([1 100 1 1e25])
    set(gca,'XScale','Log',...
            'XTick',[1 2 5 10 20 50 100],...
            'YScale','Log',...
            'YTick',10.^(0:5:25),...
            'YTickLabel',{'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'})
    if i==9
        xlabel('Rank'); 
    else
        set(gca,'XTickLabel','')
    end
    ylabel('Abundance')
    title(['(' char(96+pp) ')'])
    sh3.Position(1) = sh3.Position(1) - 0.03;
    sh3.Position(2) = sh3.Position(2) + (rp-1).*0.04 +0.02;
    set(gca,'FontSize',14)
    
    drawnow
    
end
sname=['../Figures/RankAbundance.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(67)
set(gcf,'Position',[223 1 1879 1344])
clf

pp=0;
rp=0;
for i=[3 5 7 9]
    rp=rp+1;
    
    x = xx{i};
    % sum across seed populations
    switch run_options{i}.seed_dist
        case 'nonadaptive_dispersal'
            x2 = reshape(full(x),60646,run_options{i}.nphen,run_options{i}.nlineages);
            x  = squeeze(sum(x2,2));
    end
    
    % species number
    nspecies = sum(x(:,1:94)>0,2);
    % shannon index
    lnx = log(x);
    lnx(isinf(lnx))=0;
    shannon = -sum(x.*lnx,2)
    
  
    pp=pp+1;
    sh1=subplot(4,2,pp);
    [ax] = plot_vector(nspecies,'lin',mygrid,ocean);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
%     caxis([0 1e20])
    hold on
    title(['(' char(96+pp) ')'])
    scatterm(ocean.lat(ocean.sample_points(seed_ID)),ocean.lon(ocean.sample_points(seed_ID)),25,'m')
    ch=colorbar;
    ch.Location = 'westoutside';
%     ch.TickLabels={'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'};
    set(gca,'FontSize',14)
    caxis([0 94])
  
    pp=pp+1;
    sh1=subplot(4,2,pp);
    [ax] = plot_vector(shannon,'lin',mygrid,ocean);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
%     caxis([0 1e20])
    hold on
    title(['(' char(96+pp) ')'])
    scatterm(ocean.lat(ocean.sample_points(seed_ID)),ocean.lon(ocean.sample_points(seed_ID)),25,'m')
    ch=colorbar;
    ch.Location = 'westoutside';
%     ch.TickLabels={'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'};
    set(gca,'FontSize',14)
    caxis([0 0.3])
    
    drawnow
    
end
% sname=['../Figures/RankAbundance.png'];
% set(gcf,'Color','w')
% export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1=figure(1);
clf

prc=95;
t_immigration = prctile(t_occ{1},prc,2)./prctile(t_occ{2},prc,2);
t_emmigration = prctile(t_occ{1},prc,1)./prctile(t_occ{2},prc,1);

[ax] = plot_vector(t_immigration,'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
% th=title(['Prochlorococcus']);
ch=colorbar;
colormap(flipud(redblue))
caxis(log10([1/10 10]));
ch.Ticks=log10([1/10 1/5 1/3 1/2 1 2 3 5 10]);
ch.TickLabels={'1/10','1/5','1/3','1/2','1','2','3','5','10'};

scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),25,log10(t_emmigration),'filled')
scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),25,'k','LineWidth',0.1)

sname=['../Figures/surface_vs_depth_map.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3=figure(3);
f3.Position = [75 130 560 715];
clf

D           = t_occ{1}(ocean.sample_points,:)./ t_occ{2}(ocean.sample_points,:);
D(isnan(D)) = inf;

basins=ocean.basins;
basins(find(ocean.lat<-60))=5; % southern ocean
basins(find(ocean.basins==1 & ocean.lon<180 & ocean.lat<50  & ocean.lat>30))=6; % mediterranean

basins(basins==1)=8;
basins(basins==2)=10;
basins(basins==3)=11;
basins(basins==4)=7;
basins(basins==5)=12;
basins(basins==6)=9;
basin_names{7}='Arctic';
basin_names{8}='Atlantic';
basin_names{9}='Med. Sea';
basin_names{10}='Indian';
basin_names{11}='Pacific';
basin_names{12}='S. Ocean';


[ibasins,ind] = sort(basins(ocean.sample_points));
imagesc(log10(D(ind,ind)))
hold on
isdv=find(diff(ibasins))';
axlim=xlim';
% ones(2,numel(isdv)).*axlim;
plot([axlim],repmat(isdv,2,1),'k-','LineW',1.5)
plot(repmat(isdv,2,1),[axlim],'k-','LineW',1.5)

for i=unique(ibasins)'
    ii=find(ibasins==i);
    text(mean(ii),axlim(2).*1.01,basin_names(i),...
         'HorizontalAlignment','l',...
         'Rotation',90)
    text(axlim(2).*1.01,mean(ii),basin_names(i),...
         'HorizontalAlignment','l')
     set(gca,'XTick',[],'YTick',[])
end
xlabel('source')
ylabel('destination')
axis xy

ch=colorbar('Location','SouthOutside');
ch.Label.String={' ';'Connectance Time'};
colormap(flipud(redblue))
caxis(log10([1/10 10]));
ch.Ticks=log10([1/10 1/5 1/3 1/2 1 2 3 5 10]);
ch.TickLabels={'1/10','1/5','1/3','1/2','1','2','3','5','10'};
ch.FontSize=11;
axis square

sname=['../Figures/connection_matrix.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%f5=figure(5);
f5=figure(5);
clf

clrs={'r','b'};
for i=1:2
    
    t_immigration = prctile(t_occ{i},prc,2);
    t_emmigration = prctile(t_occ{i},prc,1);
    
    if i==1
        clr=abs(ocean.lat(ocean.sample_points));
        sz=25;
    else
        clr='k';
        sz=1;
    end
    scatter(t_emmigration,t_immigration(ocean.sample_points),sz+2,'k')
    hold on
    scatter(t_emmigration,t_immigration(ocean.sample_points),sz,clr,'filled')
    caxis([0 90])
end
colormap(redblue)
colorbar
box on
xlabel('Emigration time (years)')
ylabel('Immigration time (years)')
axis square
axis([0 80 0 80])
hold on
plot(xlim,ylim,'k-')
set(gcf,'Color','w')

sname=['../Figures/imm_vs_em.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
