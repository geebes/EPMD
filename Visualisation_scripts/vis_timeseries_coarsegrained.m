% EPMD global connectivity and homogenisation through time
% Connectivity uses first occupation across native ECCO cells.
% Area-weighted homogenisation uses diagnostics when available.
clear
clc

%% Settings
% outputDir = '/Users/baw1d17/GitHub/EPMD/Output/neutral_stochastic_seasonal_GUD_X01_upper50m_buoyant_transport_nobio_years';
outputDir = '/Users/baw1d17/GitHub/EPMD/Output/neutral_stochastic_seasonal_GUD_X01_upper50m_buoyant_transport_years';

presentationMode = true; % true = black presentation, false = white paper

backgroundColour = [1 1 1];
foregroundColour = [0 0 0];
frameColour = [0 0 0];
seriesColours = get(groot,'factoryAxesColorOrder');
if presentationMode
    backgroundColour = [0 0 0];
    foregroundColour = [1 1 1];
    frameColour = [0.55 0.55 0.55];
    seriesColours = 0.7*seriesColours + 0.3;
end

%% Initial conditions

initialFile = fullfile(outputDir,'initial.mat');
if ~isfile(initialFile)
    error('Initial conditions not found: %s',initialFile);
end
load(initialFile,'x0','run_options');
daysPerYear = run_options.nday;
clear run_options
[nLocations,nPopulations] = size(x0);
nSeedPopulations = nPopulations-1; % Final population is the resident

% Identify the founding cell for each seeded population.
[seedCell,seedPop] = find(x0(:,1:nSeedPopulations));
assert(numel(seedCell)==nSeedPopulations, ...
    'Expected one initial cell per seeded population.');
[seedPop,order] = sort(seedPop);
seedCell = seedCell(order);
assert(isequal(seedPop,(1:nSeedPopulations)'), ...
    'Unexpected initial population distribution.');

%% Completed annual outputs
files = dir(fullfile(outputDir,'year_*.mat'));
completedYears = [];
for k = 1:numel(files)
    token = regexp(files(k).name,'^year_(\d{6})\.mat$','tokens','once');
    if ~isempty(token)
        completedYears(end+1,1) = str2double(token{1}); %#ok<SAGROW>
    end
end
completedYears = sort(completedYears);
if isempty(completedYears)
    error('No completed annual population files found in %s.',outputDir);
end
latestYear = completedYears(end);

%% Global connectivity from first-occupation times
latestFile = fullfile(outputDir,sprintf('year_%06d.mat',latestYear));
data = load(latestFile,'t_occupied');
t_occupied = data.t_occupied;
assert(isequal(size(t_occupied),[nLocations,nSeedPopulations]), ...
    'Unexpected dimensions for t_occupied in %s.',latestFile);

% Exclude each population's founding cell. Every other cell/population pair
% is a possible connection, counted when first occupied.
homeLinear = sub2ind(size(t_occupied),seedCell,seedPop);
t_occupied(homeLinear) = 0;
firstArrivals = nonzeros(t_occupied);
firstArrivals = firstArrivals(firstArrivals>0);
nPossibleConnections = nLocations*nSeedPopulations-nSeedPopulations;

if isempty(firstArrivals)
    connectionTime = [0; latestYear];
    connectivity = [0; 0];
else
    [arrivalTime,~,arrivalIndex] = unique(firstArrivals);
    newConnections = accumarray(arrivalIndex,1);
    connectionTime = [0; arrivalTime];
    connectivity = [0; cumsum(newConnections)/nPossibleConnections];
    if connectionTime(end)<latestYear
        connectionTime(end+1,1) = latestYear;
        connectivity(end+1,1) = connectivity(end);
    end
end
clear t_occupied firstArrivals data

%% Locate completed diagnostics
files = dir(fullfile(outputDir,'diagnostics_*.mat'));
diagnosticFiles = {};
diagnosticYears = [];
for k = 1:numel(files)
    token = regexp(files(k).name,'^diagnostics_(\d{6})\.mat$','tokens','once');
    if ~isempty(token)
        diagnosticYears(end+1,1) = str2double(token{1}); %#ok<SAGROW>
        diagnosticFiles{end+1,1} = fullfile(files(k).folder,files(k).name); %#ok<SAGROW>
    end
end
if isempty(diagnosticFiles)
    error('No completed diagnostic files found in %s.',outputDir);
end
[diagnosticYears,order] = sort(diagnosticYears);
diagnosticFiles = diagnosticFiles(order);

%% Figure 1. Global connectivity and homogenisation
[~,simName] = fileparts(outputDir);
weightedTime = [];
weightedHomogenisation = [];
for k = 1:numel(diagnosticFiles)
    diagnostic = matfile(diagnosticFiles{k});
    if ~ismember('homogenisationMetric',who(diagnostic)) || ...
            ~strcmp(diagnostic.homogenisationMetric,'area_weighted_mean_pairwise_native')
        continue
    end
    time = diagnostic.time;
    values = diagnostic.homogenisation;
    assert(numel(time)==numel(values), ...
        'Unexpected diagnostic dimensions in %s.',diagnosticFiles{k});
    time = time(:);
    values = values(:);
    positive = time>0;
    weightedTime = [weightedTime; time(positive)]; %#ok<AGROW>
    weightedHomogenisation = [weightedHomogenisation; values(positive)]; %#ok<AGROW>
end
useAnnualHomogenisation = isempty(weightedTime);
if useAnnualHomogenisation
    homogenisationTime = completedYears(:);
    homogenisation = nan(size(homogenisationTime));
else
    [homogenisationTime,keep] = unique(weightedTime);
    homogenisation = weightedHomogenisation(keep);
end

fig = figure(1);
set(fig,'Color',backgroundColour,'Name','EPMD connectivity and homogenisation');
tl = findobj(fig,'Type','tiledlayout');
if isempty(tl)
    clf(fig);
    tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
end
ax1 = nexttile(tl,1);
ax2 = nexttile(tl,2);
colororder(ax1,seriesColours);

% Use the output directory name as the simulation label.
h1 = findobj(ax1,'Type','stair','DisplayName',simName);
h2 = findobj(ax2,'Type','line','DisplayName',simName);

positiveTime = connectionTime>0;

if isempty(h1)
    hold(ax1,'on');
    h1 = stairs(ax1,connectionTime(positiveTime), ...
        connectivity(positiveTime),'LineWidth',1.5, ...
        'DisplayName',simName);
else
    delete(h1(2:end));
    h1 = h1(1);
    set(h1,'XData',connectionTime(positiveTime), ...
        'YData',connectivity(positiveTime));
end

if isempty(h2)
    hold(ax2,'on');
    h2 = plot(ax2,homogenisationTime,homogenisation,'LineWidth',1.2, ...
        'Color',h1.Color,'DisplayName',simName);
else
    delete(h2(2:end));
    h2 = h2(1);
    set(h2,'XData',homogenisationTime,'YData',homogenisation,'Color',h1.Color);
end
legend(ax2,'off');
% Keep paired series in sync when restyling an existing figure.
connectivitySeries = findobj(ax1,'Type','stair');
for k = 1:numel(connectivitySeries)
    series = connectivitySeries(k);
    paired = findobj(ax2,'Type','line','DisplayName',series.DisplayName);
    set(paired,'Color',series.Color);
end
lgd = legend(ax1,'show','Interpreter','none','Location','best');
set(lgd,'Color',backgroundColour,'TextColor',foregroundColour,'EdgeColor',foregroundColour);

set([ax1 ax2],'XScale','log','XGrid','on','YGrid','on','Box','on', ...
    'Color',backgroundColour,'XColor',foregroundColour,'YColor',foregroundColour, ...
    'GridColor',foregroundColour,'MinorGridColor',foregroundColour);
axis(ax1,'tight');
axis(ax2,'tight');

% Include every series in both tiles, using only points on the log axis.
curves = findobj([ax1 ax2],'-property','YData');
xMax = latestYear;
for k = 1:numel(curves)
    x = curves(k).XData;
    y = curves(k).YData;
    plotted = isfinite(x) & x>0 & isfinite(y);
    xMax = max([xMax; reshape(x(plotted),[],1)]);
end
set([ax1 ax2],'YLim',[0 1]);
ylabel(ax1,'Global connectivity','Color',foregroundColour);
ylabel(ax2,'Area-weighted homogenisation','Color',foregroundColour);
xlabel(ax2,'Time (model years)','Color',foregroundColour);
linkaxes([ax1 ax2],'x');

% Show subannual milestones followed by powers of ten in model years.
yearTicks = 10.^(0:floor(log10(xMax)));
timeTicks = [[1 7 30]/daysPerYear yearTicks];
timeLabels = ["day" "week" "month" compose("%.0f",yearTicks)];
set([ax1 ax2],'XLim',[timeTicks(1) xMax], ...
    'XTick',timeTicks,'XTickLabel',timeLabels);

% Place minor ticks at subannual milestones and between decade marks.
minorDays = [1:6 7:7:21 30:30:330];
minorYears = (1:9)' .* yearTicks;
minorTicks = [minorDays/daysPerYear minorYears(:)'];
minorTicks = setdiff(minorTicks(minorTicks<=xMax),timeTicks);
for timeAx = [ax1 ax2]
    timeAx.XAxis.MinorTickValues = minorTicks;
    timeAx.XMinorTick = 'on';
end
drawnow

%% Figure 2. Annual native-cell Bray–Curtis maps
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoDir,'Diagnostic_funcs'),fullfile(repoDir,'EPMD_functions'));
if isempty(which('gcmfaces_global'))
    addpath(fullfile(fileparts(repoDir),'gcmfaces'));
end
gcmfacesDir = fileparts(which('gcmfaces_global'));
addpath(fullfile(gcmfacesDir,'gcmfaces_IO'));
grid_load([fullfile(repoDir,'nctiles_grid') filesep],5,'nctiles');
gcmfaces_global

gridFields = {'iface','ix','iy','iz'};
for k = 1:numel(gridFields)
    field = gridFields{k};
    indices = h5read(initialFile,['/ocean/' field]);
    assert(numel(indices)==nLocations, ...
        'Native grid indices do not match the population rows.');
    oceanGrid.(field) = indices(:);
end
cellArea = native_ecco_cell_area(oceanGrid);
areaWeights = cellArea/sum(cellArea);

land = shaperead('landareas','UseGeoCoords',true);
landColour = [0.75 0.75 0.75];
landLatParts = arrayfun(@(shape) [shape.Lat(:); NaN],land,'UniformOutput',false);
landLonParts = arrayfun(@(shape) [shape.Lon(:); NaN],land,'UniformOutput',false);
landLat = vertcat(landLatParts{:});
landLon = vertcat(landLonParts{:});

mapFig = figure(2);
clf(mapFig);
set(mapFig,'Color',backgroundColour,'Name','Native-cell Bray–Curtis similarity');

% Match the grid to the figure and the 2:1 aspect ratio of each map.
nMaps = numel(completedYears);
candidateCols = 1:nMaps;
candidateRows = ceil(nMaps./candidateCols);
figureAspect = mapFig.Position(3)/mapFig.Position(4);
[~,bestLayout] = min(abs(log(2*candidateCols./candidateRows/figureAspect)));
ncols = candidateCols(bestLayout);
nrows = candidateRows(bestLayout);
mapLayout = tiledlayout(mapFig,nrows,ncols,'TileSpacing','tight','Padding','tight');
title(mapLayout,simName,'Interpreter','none','Color',foregroundColour);
colormap(mapFig,turbo);

annualDataSeconds = 0;
annualPlotSeconds = 0;
% Convert population fractions to concentration using the matching day's
% carrying capacity and the volume represented by each surface cell.
cellVolume = h5read(initialFile,'/ocean/volume');
cellVolume = cellVolume(:);
assert(numel(cellVolume)==nLocations && all(isfinite(cellVolume) & cellVolume>0), ...
    'Expected a positive model volume for each native cell.');
initialCapacity = h5read(initialFile,'/ocean/forcing_PCapacity',[1 1],[nLocations 1]);
yearEndCapacity = h5read(initialFile,'/ocean/forcing_PCapacity', ...
    [1 daysPerYear],[nLocations 1]);
yearEndConcentrationScale = yearEndCapacity(:)./cellVolume;
concentration = zeros(nLocations,nMaps+1);
concentration(:,1) = full(sum(x0,2)).*initialCapacity(:)./cellVolume;

annualConnectivity = interp1(connectionTime,connectivity,completedYears,'previous');
fprintf('\nGlobal metrics for %s\n',simName);
fprintf('%14s  %12s  %s\n', ...
    'Model year','Connectivity','Area-weighted homogenisation');
for k = 1:nMaps
    stageTimer = tic;
    yr = completedYears(k);
    annualFile = fullfile(outputDir,sprintf('year_%06d.mat',yr));
    annual = load(annualFile,'x');
    assert(isequal(size(annual.x),[nLocations nPopulations]), ...
        'Unexpected population dimensions in %s.',annualFile);
    [similarity,annualHomogenisation] = ...
        meanBCToOtherLocations(annual.x(:,1:nSeedPopulations),areaWeights);
    concentration(:,k+1) = full(sum(annual.x,2)).*yearEndConcentrationScale;
    clear annual
    annualDataSeconds = annualDataSeconds + toc(stageTimer);

    stageTimer = tic;
    % plot_vector sets Position, so build the map before assigning its tile.
    axes('Parent',mapFig);
    ax = plot_vector(similarity,'lin',mygrid,oceanGrid,'mollweid');
    set(findobj(ax,'Type','patch','Tag','Frame'),'EdgeColor',frameColour);
    hold(ax,'on');
    if k==1
        % Render the multipart land geometry as one fill patch.
        landOverlay = geoshow(ax,landLat,landLon,'DisplayType','polygon', ...
            'FaceVertexForm',true,'FaceColor',landColour, ...
            'EdgeColor','k','LineWidth',0.5);
    else
        copyobj(landOverlay,ax);
    end
    clim(ax,[0 1]);
    axis(ax,'off');
    title(ax,sprintf('Year %d',yr),'Color',foregroundColour);
    ax.Parent = mapLayout;
    ax.Layout.Tile = k;
    annualPlotSeconds = annualPlotSeconds + toc(stageTimer);

    if useAnnualHomogenisation
        homogenisation(k) = annualHomogenisation;
        h2.YData = homogenisation;
    end
    fprintf('%14d  %12.5f  %.5f\n',yr,annualConnectivity(k),annualHomogenisation);
    drawnow
end

stageTimer = tic;
cb = styledColorbar(ax,foregroundColour);
cb.Layout.Tile = 'east';
cb.Label.String = 'Area-weighted Bray–Curtis similarity to other locations';
cb.Ticks = 0:0.2:1;
drawnow
annualPlotSeconds = annualPlotSeconds + toc(stageTimer);
fprintf('Yearly maps: %.1f s reading and calculating, %.1f s plotting.\n', ...
    annualDataSeconds,annualPlotSeconds);

%% Figure 3. Effective divergence and transport vectors
transport = load(initialFile,'ocean');
B = transport.ocean.B;
lat = double(transport.ocean.lat(:));
lon = mod(double(transport.ocean.lon(:))+180,360)-180;
capacity = mean(transport.ocean.forcing_PCapacity,2);
stepsPerDay = (86400/(4*transport.ocean.dt_sec))*transport.ocean.nsubtime;
clear transport

% The transport of a uniform field measures effective convergence.
uniform = ones(nLocations,1);
divergence = (uniform-B*uniform)*stepsPerDay;
divergenceColour = nthroot(divergence,3);

% Convert concentration transfers to mass transfers and sum displacements
% at each source cell. Wrap longitude differences across the date line.
[dest,src,weight] = find(B);
valid = dest~=src & weight>0 & capacity(src)>0 & ...
    isfinite(capacity(src)) & isfinite(capacity(dest));
dest = dest(valid);
src = src(valid);
weight = weight(valid).*capacity(dest)./capacity(src);
dLon = mod(lon(dest)-lon(src)+180,360)-180;
dLat = lat(dest)-lat(src);
u = accumarray(src,weight.*dLon,[nLocations 1]);
v = accumarray(src,weight.*dLat,[nLocations 1]);

% Average vectors in 3-degree bins, using the reference arrow scaling.
[~,~,group] = unique(floor([lat lon]/3),'rows');
latq = accumarray(group,lat,[],@mean);
lonq = accumarray(group,lon,[],@mean);
uq = 4*accumarray(group,u,[],@mean);
vq = 4*accumarray(group,v,[],@mean);

transportFig = figure(3);
clf(transportFig);
set(transportFig,'Color',backgroundColour,'Name','Effective divergence and transport');
axes('Parent',transportFig);
transportAx = plot_vector(divergenceColour,'lin',mygrid,oceanGrid,'mollweid');
set(findobj(transportAx,'Type','patch','Tag','Frame'),'EdgeColor',frameColour);
hold(transportAx,'on');
copyobj(landOverlay,transportAx);
if presentationMode
    % Outline dark vectors so they read over dark parts of the map.
    vectorOutline = quiverm(latq,lonq,vq,uq,'w',4);
    set(vectorOutline,'LineWidth',1.5);
end
quiverm(latq,lonq,vq,uq,'k',4);
axis(transportAx,'off');
limit = max(abs(divergenceColour));
if limit==0
    limit = 1;
end
clim(transportAx,[-limit limit]);
colormap(transportAx,interp1([-1; 0; 1], ...
    [0.15 0.35 0.75; 1 1 1; 0.75 0.15 0.15],linspace(-1,1,257)'));
cb = styledColorbar(transportAx,foregroundColour);
cb.Label.String = 'Cube root of effective divergence ((day^{-1})^{1/3})';
title(transportAx,'Effective divergence and transport vectors','Color',foregroundColour);



%% Figure 4. Diagnostic Bray–Curtis similarity movie
animationFile = fullfile(outputDir,'spatial_BC_similarity.mp4');
movieFig = figure(4);
clf(movieFig);
set(movieFig,'Color',backgroundColour,'Name','Diagnostic Bray–Curtis similarity');
axes('Parent',movieFig);
[movieAx,movieSurfaces] = plot_vector((1:nLocations)','lin',mygrid,oceanGrid,'mollweid');

% Reuse the projected mesh. Its initial colours identify the native rows.
nativeRows = get(movieSurfaces,'CData');
mapped = cellfun(@isfinite,nativeRows,'UniformOutput',false);
frameColours = cellfun(@(rows) nan(size(rows)),nativeRows,'UniformOutput',false);
set(findobj(movieAx,'Type','patch','Tag','Frame'),'EdgeColor',frameColour);
hold(movieAx,'on');
copyobj(landOverlay,movieAx);
[seedX,seedY] = projfwd(getm(movieAx),lat(seedCell),lon(seedCell));
seedMarkers = scatter(movieAx,seedX,seedY,10,'o', ...
    'MarkerFaceColor','w','MarkerEdgeColor','k');
clim(movieAx,[0 1]);
colormap(movieAx,turbo);
cb = styledColorbar(movieAx,foregroundColour);
cb.Label.String = 'Diagnostic Bray–Curtis similarity';
movieTitle = title(movieAx,'','Color',foregroundColour);

fprintf('Writing diagnostic movie to %s\n',animationFile);
movieTimer = tic;
movieWriter = VideoWriter(animationFile,'MPEG-4');
movieWriter.FrameRate = 10;
open(movieWriter);
frameNumber = 0;
try
    for k = 1:numel(diagnosticFiles)
        data = load(diagnosticFiles{k},'time','bcMap');
        if ismember('bcMapMetric',who('-file',diagnosticFiles{k}))
            metadata = load(diagnosticFiles{k},'bcMapMetric');
            switch metadata.bcMapMetric
                case 'area_weighted_mean_pairwise_native'
                    cb.Label.String = 'Area-weighted Bray–Curtis similarity to other locations';
                case 'mean_pairwise_native'
                    cb.Label.String = 'Equal-cell mean Bray–Curtis similarity to other locations';
                otherwise
                    error('Unknown Bray–Curtis metric in %s.',diagnosticFiles{k});
            end
        else
            cb.Label.String = 'Bray–Curtis similarity to global mean';
        end
        assert(isequal(size(data.bcMap),[nLocations numel(data.time)]), ...
            'Unexpected diagnostic map dimensions in %s.',diagnosticFiles{k});
        for j = 1:numel(data.time)
            for face = 1:numel(movieSurfaces)
                rows = nativeRows{face}(mapped{face});
                frameColours{face}(mapped{face}) = double(data.bcMap(rows,j));
                movieSurfaces(face).CData = frameColours{face};
            end
            movieTitle.String = sprintf('Year %.5g',data.time(j));
            seedAlpha = max(0,1-data.time(j)*daysPerYear/30);
            set(seedMarkers,'MarkerFaceAlpha',seedAlpha,'MarkerEdgeAlpha',seedAlpha);
            drawnow
            writeVideo(movieWriter,getframe(movieFig));
            frameNumber = frameNumber+1;
        end
    end
catch err
    close(movieWriter);
    rethrow(err)
end
close(movieWriter);
fprintf('Saved %d diagnostic frames in %.1f s to %s\n', ...
    frameNumber,toc(movieTimer),animationFile);

%% Figure 5. Total local concentration movie
positiveConcentration = isfinite(concentration) & concentration>0;
concentration(~positiveConcentration) = NaN;

concentrationFile = fullfile(outputDir,'total_local_concentration.mp4');
concentrationFig = figure(5);
clf(concentrationFig);
set(concentrationFig,'Color',backgroundColour,'Name','Total local concentration');
axes('Parent',concentrationFig);
[concentrationAx,concentrationSurfaces] = ...
    plot_vector(concentration(:,1),'lin',mygrid,oceanGrid,'mollweid');
set(findobj(concentrationAx,'Type','patch','Tag','Frame'),'EdgeColor',frameColour);
hold(concentrationAx,'on');
copyobj(landOverlay,concentrationAx);
colormap(concentrationAx,turbo);
concentrationAx.ColorScale = 'log';
clim(concentrationAx,[1e3 1e15]);
cb = styledColorbar(concentrationAx,foregroundColour);
cb.Ticks = 10.^(3:15);
cb.TickLabels = compose('10^{%d}',3:15);
cb.TickLabelInterpreter = 'tex';
cb.Label.String = 'Total concentration (cells m^{-3})';
concentrationTitle = title(concentrationAx,'','Color',foregroundColour);
concentrationYears = [0; completedYears(:)];

fprintf('Writing %d concentration frames to %s\n',numel(concentrationYears),concentrationFile);
writeConcentrationMovie(concentrationFile,concentrationFig,concentrationTitle, ...
    concentrationSurfaces,nativeRows,mapped,concentration,concentrationYears);

function cb = styledColorbar(ax,foregroundColour)
cb = colorbar(ax);
cb.Color = foregroundColour;
cb.Label.Color = foregroundColour;
end

function writeConcentrationMovie(filename,fig,titleHandle,surfaces,nativeRows,mapped,values,years)
writer = VideoWriter(filename,'MPEG-4');
writer.FrameRate = 2;
open(writer);
cleanup = onCleanup(@() close(writer));
colours = cellfun(@(rows) nan(size(rows)),nativeRows,'UniformOutput',false);
for k = 1:numel(years)
    for face = 1:numel(surfaces)
        rows = nativeRows{face}(mapped{face});
        colours{face}(mapped{face}) = values(rows,k);
        surfaces(face).CData = colours{face};
    end
    titleHandle.String = sprintf('Year %d',years(k));
    drawnow
    writeVideo(writer,getframe(fig));
end
end

function [similarity,homogenisation] = meanBCToOtherLocations(abundance,areaWeights)
totals = full(sum(abundance,2));
occupied = totals>0;
distanceSum = [];
if nnz(occupied)>1
    distanceSum = native_pairwise_distance_sum(abundance,totals,areaWeights);
end
[similarity,homogenisation] = native_pairwise_similarity(distanceSum,occupied,areaWeights);
end
