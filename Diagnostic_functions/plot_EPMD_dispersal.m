clear
% clc
addpath(genpath('~/GitHub/EPMD'))
diag_fcns = diagnostics;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_filename = {   'neutral_stochastic_static_GUD_X01_surface_transport',...
                     'neutral_stochastic_static_GUD_X01_surface_transport_slow',...
                     'neutral_stochastic_static_GUD_X01_weighted_transport',...
                     'neutral_stochastic_static_GUD_X01_weighted_transport_slow',...
       'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport',...
       'nonadaptive_dispersal_stochastic_static_GUD_X01_weighted_transport_slow',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.01_slow',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1',...
         'selective_dispersal_stochastic_static_GUD_X01_weighted_transport_m0.1_slow',...
    ...
                     'neutral_stochastic_static_GUD_X17_surface_transport',...
                     'neutral_stochastic_static_GUD_X17_surface_transport_slow',...
                     'neutral_stochastic_static_GUD_X17_weighted_transport',...
                     'neutral_stochastic_static_GUD_X17_weighted_transport_slow',...
       'nonadaptive_dispersal_stochastic_static_GUD_X17_weighted_transport',...
       'nonadaptive_dispersal_stochastic_static_GUD_X17_weighted_transport_slow',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.01',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.01_slow',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.1',...
         'selective_dispersal_stochastic_static_GUD_X17_weighted_transport_m0.1_slow'};
  
input_filename = input_filename{13};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathname   = '~/GitHub/EPMD/Output/';
matObj  = matfile([pathname input_filename '.mat']);

ocean       = matObj.ocean;
run_options = matObj.run_options;
t_occupied  = matObj.t_occupied;

i_lastyr    = matObj.yrs_saved;
disp([num2str(i_lastyr) ' years evaluated.']) 

grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
gcmfaces_global
load coastlines
land = shaperead('landareas', 'UseGeoCoords', true);

ocean.ann_abundance = mean(ocean.forcing_PCapacity,2);

K = ocean.forcing_PCapacity(:,end);
T = ocean.forcing_temp;
t_occ=full(t_occupied(:,1:numel(ocean.sample_points)));
   
disp(['ocean is ' num2str(100.*nnz(t_occ)/numel(t_occ),'%2.0f') '% connected.'])

t_occ(~t_occ(:))=NaN;


if exist([pathname input_filename])~=7
    mkdir([pathname input_filename]);
end

% return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=figure(1);
clf

mean_abundance = mean(ocean.forcing_PCapacity,2);
[ax] = plot_vector(mean_abundance,'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
caxis([0 25]);

% transport vectors
[x,y,u,v] = get_circ_vectors(ocean); % second input is coarse graining resolution
iplot=randsample(60646,1e4);
h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');

% th=title(['Prochlorococcus']);
colormap(parula)
ch=colorbar('Location','SouthOutside');
ch.TickLabels={'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'};
drawnow

sname=[pathname input_filename '/abundance.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f11=figure(11);
clf

mean_temp = mean(ocean.forcing_temp,2);
[ax] = plot_vector(mean_temp,'lin',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
caxis([-2 36]);

% transport vectors
[x,y,u,v] = get_circ_vectors(ocean); % second input is coarse graining resolution
iplot=randsample(60646,1e4);
h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');
% h(1).Color=[1 1 1].*0.1;
% h(2).Color=[1 1 1].*0.1;

th=title(['Mean Annual SST (^\circC)']);
colormap(parula)
ch=colorbar('Location','SouthOutside');
ch.Ticks=-2:2:36;
drawnow

sname=[pathname input_filename '/temp_and_circulation.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2=figure(2);
f2.Position = [209 369 930 685];
clf
t_occ(isnan(t_occ))=100;
prc=95;
t_immigration = prctile(t_occ,prc,2);
t_emmigration = prctile(t_occ,prc,1);

[ax] = plot_vector(t_immigration,'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!

% transport vectors
[x,y,u,v] = get_circ_vectors(ocean); % second input is coarse graining resolution
iplot=randsample(60646,1e4);
h=quiverm(y(iplot),x(iplot),v(iplot),u(iplot),'k');

% th=title(['Prochlorococcus']);
% ch=colorbar;
colormap(flipud(turbo)) 
caxis(log10([2 100]));
ch.Ticks=log10([1 2 5 10 20 50 100 200]);
ch.TickLabels={'1','2','5','10','20','50','100','200'};

scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),50,log10(t_emmigration),'filled')
scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),50,'k','LineWidth',0.1)

sname=[pathname input_filename '/connection_times_map.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2=figure(22);
f22.Position = [209 369 930 685];
clf

ch=colorbar('Location','SouthOutside');

colormap(flipud(turbo(64)))
caxis(log10([1/365 100]));
ch.Ticks=log10([1/365 7/365 1/12 1 10 100]);
ch.TickLabels={'day','week','month','year','decade','century'};
axis off

sname=[pathname input_filename '/colorbar.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3=figure(3);
f3.Position = [75 130 560 715];
clf

D           = t_occ(ocean.sample_points,:);
% D(isnan(D)) = inf;
indx        = find(speye(size(D)));
D(indx)     = 0; 

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
pclr=log10(D(ind,ind));
pclr=padarray(pclr,[1 1],'post'); % pad array so no data lost by pcolor
pcolor(0:size(D,1),0:size(D,1),pclr)
shading flat
hold on
isdv=find(diff(ibasins))';
axlim=xlim';
ones(2,numel(isdv)).*axlim;
plot([axlim],repmat(isdv,2,1),'k-','LineW',1.5)
plot(repmat(isdv,2,1),[axlim],'k-','LineW',1.5)
set(gca,'LineW',1.5,'layer','top')

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
colormap(flipud(turbo))
caxis(log10([1/365 100]));
ch.Ticks=log10([1/365 7/365 1/12 1 10 100]);
ch.TickLabels={'day','week','month','year','decade','century'};
ch.FontSize=11;
axis square


sname=[pathname input_filename '/connection_matrix.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5=figure(5);
clf

h=scatter(t_emmigration,t_immigration(ocean.sample_points),...
        25,'k','Marker','.');
h=scatter(t_emmigration,t_immigration(ocean.sample_points),...
        25,abs(ocean.lat(ocean.sample_points)),'filled','Marker','s');
colormap(redblue)
colorbar

box on
xlabel('Emigration time')
ylabel('Immigration time')
axis square
axis([0 50 0 50])
hold on
plot([0 100],[0 100],'k-')

sname=[pathname input_filename '/imm_vs_em.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4 = figure(4);
clf
hold on

[N,edges] = histcounts(t_occ(:),0:(1/96):i_lastyr);

clr=rand(1,3);
plot([0 edges],[NaN cumsum(N)./numel(t_occ) NaN],'Color',clr,'LineW',2)
ylim([0 1])
box on
xlabel('Time (years)')
ylabel('Connectance')
set(gca,'XScale','log')

sname=[pathname input_filename '/cumulative_connections.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

Cseed=t_occ(ocean.sample_points,:);
Cdiag=find(speye(size(Cseed)));
Cseed(Cdiag)=NaN;
disp(['Median connection time is ' num2str(nanmedian(Cseed(:))) ' years'])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f8 = figure(8);
clf

iseed=1; % 261,279
[ax] = plot_vector(t_occ(:,iseed),'log',mygrid,ocean);
geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
% th=title(['Prochlorococcus']);

ch=colorbar('Location','SouthOutside');
ch.Label.String={'Connectance Time'};
colormap(flipud(turbo(16)))
caxis(log10([1/365 100]));
ch.Ticks=log10([1/365 7/365 1/12 1 10 100]);
ch.TickLabels={'day','week','month','year','decade','century'};
ch.FontSize=11;

scatterm(ocean.lat(ocean.sample_points(iseed)),ocean.lon(ocean.sample_points(iseed)),25,'m')

sname=[pathname input_filename '/dispersal_map.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(6)

seed_ID=8;

% get abundance data

for iyr = i_lastyr
    clf
    x  = cell2mat(matObj.x(iyr,1)) .* ocean.ann_abundance;
    
    [~,cols] = find(x);
    cols=unique(cols);
    x_i=x(:,cols(seed_ID));
    
    subplot(3,1,1:2)
    [ax] = plot_vector(x_i,'lin',mygrid,ocean);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
%     caxis([0 25])
    hold on
    title('(a)')
    scatterm(ocean.lat(ocean.sample_points(seed_ID)),ocean.lon(ocean.sample_points(seed_ID)),25,'m')
    ch=colorbar;
%     ch.TickLabels={'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'};
    
    sh2=subplot(313);
    [n,edges,bin] = histcounts(ocean.ann_theta,unique(run_options.T_opt));
    for i=1:numel(edges)
        freq(i)=sum(x_i(bin==i));
    end
    bar(edges,freq,1)
    hold on
    in9999=cumsum(freq)>sum(freq).*0.01 & cumsum(freq)<sum(freq).*0.99;
    bar(edges,freq.*in9999,1)
    mask=edges==run_options.T_opt(cols(seed_ID));
    bar(edges,freq.*mask,1)
    set(gca,'YScale','log')
%     ylim([1 1e25])
    xlabel('Temperature (^\circC)')
    ylabel('Abundance')
    title('(b)')
    sh2.Position(2)=sh2.Position(2)+0.1;
    legend('All','1% to 99%','T = T_{opt}')
    set(gca,'YTick',10.^(0:5:25),...
            'YTickLabel',{'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'})
    
    drawnow
    
    sname=[pathname input_filename '/Seed_' num2str(seed_ID,'%03i') '_Year_'  num2str(iyr,'%03i') '.png'];
    set(gcf,'Color','w')
    export_fig(sname,'-r300')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(9)
cmap=lhsdesign(94,3);

% plot dominance patterns

for iyr = 1:i_lastyr
    clf
    
    x  = cell2mat(matObj.x(iyr,1));
    % sum across seed populations
    switch run_options.seed_dist
        case 'nonadaptive_dispersal'
            x2 = reshape(full(x),60646,run_options.nphen,run_options.nlineages);
            x  = squeeze(sum(x2,2));
    end
    [~,I] = max(x(:,1:94),[],2);
    
    [ax] = plot_vector(I,'lin',mygrid,ocean);
    colormap(cmap)
    scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),50,cmap,'filled')
    scatterm(ocean.lat(ocean.sample_points),ocean.lon(ocean.sample_points),50,'k')
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!
    
    title(['Year ' num2str(iyr)])
    
    drawnow
    
    sname=[pathname input_filename '/Dominance_Year_'  num2str(iyr,'%03i') '.png'];
    set(gcf,'Color','w')
    export_fig(sname,'-r300')
end


%%
figure(99)
    clf

% plot dominance patterns
cmap=turbo(100);
for iyr = 1:i_lastyr
    
    x  = cell2mat(matObj.x(iyr,1));
    x = x .* K;
    
    
    % sum across seed populations
    switch run_options.seed_dist
        case 'nonadaptive_dispersal'
            x2 = reshape(full(x),60646,run_options.nphen,run_options.nlineages);
            x  = squeeze(sum(x2,2));
    end
    
    x=x(ocean.sample_points(8),:);
    
    ab=sort(x(1:94),2,'descend');
    plot(1:94,ab','Color',cmap(iyr,:),'LineW',1);
    hold on
    %     plot(1:94,ab(seed_ID,:)','r-','LineW',1);
    axis([1 100 1 1e25])
    set(gca,'XScale','Log',...
        'XTick',[1 2 5 10 20 50 100],...
        'YScale','Log',...
        'YTick',10.^(0:5:25),...
        'YTickLabel',{'10^{0}','10^{5}','10^{10}','10^{15}','10^{20}','10^{25}'})
    xlabel('Rank');
    set(gca,'XTickLabel','')
    ylabel('Abundance')
    title(iyr)
    
    drawnow

end  
    
return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
clf
D           = t_occ(ocean.sample_points,:);
D(isnan(D)) = inf;
indx        = find(speye(size(D)));
D(indx)     = 0; 

triuD=triu(D);

tree = linkage(triuD);

Dv=squareform(triuD);

leafOrder = optimalleaforder(tree,Dv,'Criteria','group','Transformation','inverse');

[H, nodes, outperm] = dendrogram(tree,341,'ReOrder',leafOrder);
























