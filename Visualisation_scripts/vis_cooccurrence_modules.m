% Milke-style lineage co-occurrence modules in saved neutral EPMD output.
% Fit modules at year 1 and track those fixed lineage groups through time.
% A dot's colour represents one draw from the cell's local module mixture.

%% Analysis settings


outputDir = 'neutral_stochastic_seasonal_GUD_X01_upper50m_buoyant_transport_nobio_years';
% outputDir = 'neutral_stochastic_seasonal_GUD_X01_upper50m_buoyant_transport_years';

repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoDir,'EPMD_functions'));
outputPath = fullfile(repoDir,'Output',outputDir);
fitCells = 384;  % Area-weighted sample of native cells
nPermutations = 999;
minCorrelation = 0.5;
maxAdjustedP = 0.05;
fractionFloor = 1e-8;   % For log ratios of continuous model fractions

%% Plot settings
presentationMode = true; % true = black presentation, false = white paper
plotAnnualMaps = true;
exportMovie = true;
dotsPerCell = 15; % Average across wet cells. Each cell's allocation scales with its area.
annualMarkerSize = 0.4;
seedMarkerSize = 5;
movieMarkerSize = 1.4;
frameRate = 4;
paletteSeed = 123458;

%% Fit modules and aggregate abundances
% Run once for the selected output. Keep the results for the plot sections.
[years,module,ocean,abundanceByYear] = analyseModules(outputPath,fitCells, ...
    nPermutations,minCorrelation,maxAdjustedP,fractionFloor);

%% Prepare plots
% Rerun this section after changing the dot density or palette.
if plotAnnualMaps || exportMovie
    nMaps = numel(years);
    palette = modulePalette(max(module),paletteSeed);
    % Keep both dot locations and colour quantiles fixed throughout the movie.
    stream = RandStream('mt19937ar','Seed',123457);
    cellArea = native_ecco_cell_area(ocean);
    nCells = numel(cellArea);
    nDots = nCells*dotsPerCell;
    % Systematic rounding keeps the dot budget exact and the expected
    % allocation proportional to area, even for cells smaller than one dot.
    boundaries = [0; cumsum(cellArea)/sum(cellArea)*nDots];
    boundaries(end) = nDots;
    offsets = floor(boundaries+rand(stream,1,1));
    counts = diff(offsets);
    dotCell = repelem((1:nCells)',counts);
    withinCell = (0:nDots-1)'-offsets(dotCell);
    quantiles = (withinCell+rand(stream,nDots,1))./counts(dotCell);
    [dotLat,dotLon] = nativeDotCoordinates(ocean,repoDir,dotCell,stream);
    dotIDs = cell(nMaps+1,1);
    for k = 1:numel(dotIDs)
        moduleAbundance = abundanceByYear{k};
        totals = sum(moduleAbundance,2);
        ids = ones(nDots,1,'uint16');
        cumulative = zeros(nCells,1);
        dotTotals = totals(dotCell);
        draws = quantiles.*dotTotals;
        for j = 1:size(moduleAbundance,2)-1
            cumulative = cumulative + moduleAbundance(:,j);
            ids = ids + uint16(draws>=cumulative(dotCell));
        end
        ids(dotTotals<=0) = 0;
        dotIDs{k} = ids;
    end
    land = shaperead('landareas','UseGeoCoords',true);
    latParts = arrayfun(@(shape) [shape.Lat(:); NaN],land,'UniformOutput',false);
    lonParts = arrayfun(@(shape) [shape.Lon(:); NaN],land,'UniformOutput',false);
    landLat = vertcat(latParts{:});
    landLon = vertcat(lonParts{:});
    clear landOverlay
end

%% Figure appearance
% Rerun this section and the figures after switching presentationMode.
if plotAnnualMaps || exportMovie
    backgroundColour = [1 1 1];
    foregroundColour = [0 0 0];
    frameColour = [0 0 0];
    plotPalette = palette;
    unassignedLabel = 'Black = unassigned';
    if presentationMode
        backgroundColour = [0 0 0];
        foregroundColour = [1 1 1];
        frameColour = [0.55 0.55 0.55];
        plotPalette(1,:) = [0.55 0.55 0.55];
        plotPalette(2:end,:) = 0.75*palette(2:end,:) + 0.25;
        unassignedLabel = 'Grey = unassigned';
    end
end

%% Figure 2 equivalent. Annual maps
if plotAnnualMaps
    mapFig = figure(2);
    clf(mapFig);
    set(mapFig,'Color',backgroundColour,'Name','Annual co-occurrence modules');
    cols = 1:nMaps;
    rows = ceil(nMaps./cols);
    figureAspect = mapFig.Position(3)/mapFig.Position(4);
    [~,layoutIndex] = min(abs(log(2*cols./rows/figureAspect)));
    mapRows = max(1,rows(layoutIndex)-2);
    mapCols = max(1,cols(layoutIndex)-1);
    layout = tiledlayout(mapFig,mapRows,mapCols, ...
        'TileSpacing','tight','Padding','tight');
    title(layout,['Year-1 co-occurrence modules | ' unassignedLabel ' lineages'], ...
        'Color',foregroundColour);
    % Spread the available tiles across the saved snapshots.
    mapIndices = round(linspace(1,nMaps,min(nMaps,mapRows*mapCols)));
    for tile = 1:numel(mapIndices)
        k = mapIndices(tile);
        ax = moduleMapAxes(mapFig,frameColour);
        if tile==1
            [dotX,dotY] = projectDots(getm(ax),dotLat(:),dotLon(:));
            dotX = single(dotX);
            dotY = single(dotY);
        end
        ids = dotIDs{k+1};
        use = ids>0;
        scatter(ax,dotX(use),dotY(use),annualMarkerSize,plotPalette(ids(use),:), ...
            'filled','MarkerEdgeColor','none');
        if tile==1
            landOverlay = plotModuleLand(ax,landLat,landLon);
        else
            copyobj(landOverlay,ax);
        end
        title(ax,sprintf('Year %d',years(k)),'Color',foregroundColour);
        ax.Parent = layout;
        ax.Layout.Tile = tile;
    end
end

%% Figure 4 equivalent. Movie
if exportMovie
    movieFig = figure(4);
    clf(movieFig);
    set(movieFig,'Color',backgroundColour,'Name','Co-occurrence module movie');
    ax = moduleMapAxes(movieFig,frameColour);
    reuseMap = exist('landOverlay','var') && isgraphics(landOverlay);
    if ~reuseMap
        [dotX,dotY] = projectDots(getm(ax),dotLat(:),dotLon(:));
        dotX = single(dotX);
        dotY = single(dotY);
    end
    ids = dotIDs{1};
    use = ids>0;
    dots = scatter(ax,dotX(use),dotY(use),seedMarkerSize,plotPalette(ids(use),:), ...
        'filled','MarkerEdgeColor','none');
    if reuseMap
        copyobj(landOverlay,ax);
    else
        plotModuleLand(ax,landLat,landLon);
    end
    frameTitle = title(ax,'Day 0 (seed locations)','Color',foregroundColour);
    writeModuleMovie(fullfile(outputPath,'cooccurrence_modules_year1.mp4'),movieFig, ...
        dots,frameTitle,years,dotIDs(2:end),dotX,dotY,plotPalette,frameRate, ...
        movieMarkerSize,unassignedLabel);
end

%%

function [years,module,ocean,abundanceByYear] = analyseModules(outputDir,fitCells, ...
    nPermutations,minCorrelation,maxAdjustedP,fractionFloor)
initialFile = fullfile(outputDir,'initial.mat');
initial = load(initialFile,'x0');
[nCells,nPopulations] = size(initial.x0);
nLineages = nPopulations-1; % Exclude the resident throughout.
fields = {'iface','ix','iy'};
ocean = struct();
for k = 1:numel(fields)
    field = fields{k};
    values = h5read(initialFile,['/ocean/' field]);
    assert(numel(values)==nCells,'Grid and population rows do not match.');
    ocean.(field) = values(:);
end
files = dir(fullfile(outputDir,'year_*.mat'));
files = files(~cellfun(@isempty,regexp({files.name},'^year_\d{6}\.mat$','once')));
[~,order] = sort({files.name});
files = files(order);
assert(~isempty(files),'No completed population snapshots in %s.',outputDir);
years = arrayfun(@(file) str2double(file.name(6:11)),files);
nMaps = numel(years);

referenceIndex = find(years==1,1);
assert(~isempty(referenceIndex),'Year 1 snapshot (year_000001.mat) is required in %s.',outputDir);

% Fit the network from occupied year-1 cells only.
reference = loadLineages(files(referenceIndex),nCells,nPopulations);
cellArea = native_ecco_cell_area(ocean);
stream = RandStream('mt19937ar','Seed',123456);
[~,areaOrder] = sort(-log(rand(stream,nCells,1))./cellArea);
totals = full(sum(reference,2));
rows = areaOrder(totals(areaOrder)>0);
rows = rows(1:min(fitCells,numel(rows)));
fractions = full(reference(rows,:))./totals(rows);
profiles = log(max(fractions,fractionFloor));
if ~isempty(rows)
    profiles = profiles-mean(profiles,1);
end
keep = sum(fractions>fractionFloor,1)>=5 & sum(profiles.^2,1)>eps;
fprintf('Year 1: fitting co-occurrence on %d lineages from %d cells.\n', ...
    nnz(keep),numel(rows));
module = zeros(nLineages,1);
if nnz(keep)>=4
    fitStream = RandStream('mt19937ar','Seed',123457);
    module(keep) = inferModules(profiles(:,keep), ...
        nPermutations,minCorrelation,maxAdjustedP,fitStream);
end
nModules = max(module);
fprintf('Year 1: found %d modules containing %d lineages.\n',nModules,nnz(module));
if nModules==0
    warning('No year-1 modules passed the filters. All seeded lineages are unassigned.');
end

% Use the same lineage-to-module mapping for day zero and every saved year.
membership = sparse((1:nLineages)',module+1,ones(nLineages,1), ...
    nLineages,nModules+1);
abundanceByYear = cell(nMaps+1,1);
abundanceByYear{1} = full(initial.x0(:,1:nLineages)*membership);
abundanceByYear{referenceIndex+1} = full(reference*membership);
clear initial reference
for k = setdiff(1:nMaps,referenceIndex)
    A = loadLineages(files(k),nCells,nPopulations);
    abundanceByYear{k+1} = full(A*membership);
    fprintf('Aggregated year-1 modules for year %d (%d/%d).\n',years(k),k,nMaps);
end
end

function palette = modulePalette(nModules,seed)
previous = rng;
cleanup = onCleanup(@() rng(previous));
rng(seed,'twister');
palette = [0 0 0];
if nModules>0
    palette = [palette; lhsdesign(nModules,3)];
end
palette = single(palette);
end

function overlay = plotModuleLand(ax,lat,lon)
overlay = geoshow(ax,lat,lon,'DisplayType','polygon','FaceVertexForm',true, ...
    'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0.35 0.35 0.35],'LineWidth',0.35);
end

function A = loadLineages(file,nCells,nPopulations)
filename = fullfile(file.folder,file.name);
data = load(filename,'x');
assert(isequal(size(data.x),[nCells nPopulations]), ...
    'Expected native rows and a final resident column in %s.',filename);
A = data.x(:,1:end-1);
end

function module = inferModules(profiles,nPermutations,minCorrelation,maxAdjustedP,stream)
% SparCC basis-correlation estimator on fixed-floor continuous fractions.
% Shuffle each lineage independently across sampled cells for the null.
[nRows,nLineages] = size(profiles);
offsets = (0:nLineages-1)*nRows;
rho = sparccCorrelation(profiles'*profiles);
[i,j] = find(triu(true(nLineages),1));
pairs = sub2ind(size(rho),i,j);
observed = rho(pairs);
module = zeros(nLineages,1);
if ~any(observed>=minCorrelation)
    return
end
exceedances = zeros(size(observed));
for b = 1:nPermutations
    [~,order] = sort(rand(stream,nRows,nLineages),1);
    shuffled = profiles(order+offsets);
    null = sparccCorrelation(shuffled'*shuffled);
    exceedances = exceedances+(abs(null(pairs))>=abs(observed));
    if mod(b,100)==0 || b==nPermutations
        fprintf('Co-occurrence permutations: %d/%d\n',b,nPermutations);
    end
end
p = (1+exceedances)/(1+nPermutations);
sortedP = sort(p);
last = find(sortedP<=maxAdjustedP*(1:numel(p))'/numel(p),1,'last');
pass = false(size(p));
if ~isempty(last)
    pass = p<=sortedP(last) & observed>=minCorrelation;
end
G = graph(i(pass),j(pass),[],nLineages);
active = degree(G)>0;
fprintf('Partitioning %d connected lineages with %d edges\n',nnz(active),numedges(G));
if any(active)
    module(active) = edgeBetweennessPartition(subgraph(G,active));
end
end

function rho = sparccCorrelation(crossProduct)
% The common sample-variance divisor cancels in correlations.
variance = diag(crossProduct);
logRatioVariance = max(0,variance+variance'-2*crossProduct);
n = numel(variance);
allowed = ~eye(n);
for iteration = 0:10
    basisVariance = (diag(sum(allowed,2))+allowed)\sum(logRatioVariance.*allowed,2);
    basisVariance = max(basisVariance,eps(max(variance)));
    scale = sqrt(basisVariance);
    rho = (basisVariance+basisVariance'-logRatioVariance)./(2*(scale*scale'));
    [strongest,index] = max(triu(abs(rho).*allowed,1),[],'all','linear');
    if strongest<=0.1 || iteration==10
        break
    end
    [a,b] = ind2sub([n n],index);
    allowed(a,b) = false;
    allowed(b,a) = false;
    if any(sum(allowed,2)<2)
        break
    end
end
rho = min(1,max(-1,rho));
end

function best = edgeBetweennessPartition(original)
% Girvan–Newman, retaining the cut with maximum original-graph modularity.
G = original;
ends = original.Edges.EndNodes;
m = numedges(original);
degrees = degree(original);
bestScore = -Inf;
while true
    groups = conncomp(G)';
    inside = nnz(groups(ends(:,1))==groups(ends(:,2)))/m;
    if inside<=bestScore
        break
    end
    expected = accumarray(groups,degrees)/(2*m);
    score = inside-sum(expected.^2);
    if score>bestScore
        best = groups;
        bestScore = score;
    end
    if numedges(G)==0
        break
    end
    centrality = edgeBetweenness(G);
    mostCentral = max(centrality);
    G = rmedge(G,find(centrality>=mostCentral-max(1,mostCentral)*1e-12));
end
end

function scores = edgeBetweenness(G)
% Brandes accumulation on shortest paths in the unweighted graph.
n = numnodes(G);
ends = G.Edges.EndNodes;
m = size(ends,1);
scores = zeros(m,1);
edgeID = sparse([ends(:,1); ends(:,2)],[ends(:,2); ends(:,1)], ...
    [(1:m)'; (1:m)'],n,n);
neighbours = cell(n,1);
for v = 1:n
    neighbours{v} = find(edgeID(:,v));
end
for s = 1:n
    depth = -ones(n,1);
    depth(s) = 0;
    paths = zeros(n,1);
    paths(s) = 1;
    queue = zeros(n,1);
    queue(1) = s;
    tail = 1;
    head = 1;
    while head<=tail
        v = queue(head);
        adjacent = neighbours{v};
        unseen = adjacent(depth(adjacent)<0);
        depth(unseen) = depth(v)+1;
        queue(tail+(1:numel(unseen))) = unseen;
        tail = tail+numel(unseen);
        next = adjacent(depth(adjacent)==depth(v)+1);
        paths(next) = paths(next)+paths(v);
        head = head+1;
    end
    dependencies = zeros(n,1);
    for k = tail:-1:2
        v = queue(k);
        adjacent = neighbours{v};
        prev = adjacent(depth(adjacent)==depth(v)-1);
        credit = paths(prev)/paths(v)*(1+dependencies(v));
        dependencies(prev) = dependencies(prev)+credit;
        edges = full(edgeID(prev,v));
        scores(edges) = scores(edges)+credit;
    end
end
end

function [lat,lon] = nativeDotCoordinates(ocean,repoDir,dotCell,stream)
% Bilinear interpolation of the four corner unit vectors avoids longitude
% wrap at the antimeridian and handles the rotated Arctic face.
faceSizes = [90 270; 90 270; 90 90; 270 90; 270 90];
nDots = numel(dotCell);
u = rand(stream,nDots,1);
v = rand(stream,nDots,1);
lat = nan(nDots,1);
lon = lat;
dotFace = ocean.iface(dotCell);
for face = 1:size(faceSizes,1)
    selected = dotFace==face;
    dims = faceSizes(face,:)+1;
    filename = fullfile(repoDir,'nctiles_grid',sprintf('tile%03d.mitgrid',face));
    [file,message] = fopen(filename,'r','ieee-be');
    assert(file>=0,'Cannot read %s: %s',filename,message);
    cleanup = onCleanup(@() fclose(file));
    assert(fseek(file,5*prod(dims)*8,'bof')==0,'Cannot seek to ECCO grid corners.');
    [corners,count] = fread(file,[prod(dims) 2],'double=>double');
    assert(count==2*prod(dims),'Incomplete ECCO grid corners.');
    clear cleanup
    [x,y,z] = sph2cart(deg2rad(corners(:,1)),deg2rad(corners(:,2)),1);
    cells = dotCell(selected);
    first = sub2ind(dims,ocean.ix(cells),ocean.iy(cells));
    indices = first+[0 1 dims(1) dims(1)+1];
    a = u(selected,:);
    b = v(selected,:);
    weights = [(1-a).*(1-b),a.*(1-b),(1-a).*b,a.*b];
    x = sum(weights.*reshape(x(indices),size(indices)),2);
    y = sum(weights.*reshape(y(indices),size(indices)),2);
    z = sum(weights.*reshape(z(indices),size(indices)),2);
    [lambda,phi] = cart2sph(x,y,z);
    lon(selected,:) = rad2deg(lambda);
    lat(selected,:) = rad2deg(phi);
end
end

function ax = moduleMapAxes(fig,frameColour)
ax = axes('Parent',fig);
axes(ax);
ax = axesm('mollweid','MapLonLimit',[-245 115],'MapLatLimit',[-90 90], ...
    'Frame','on','FLineWidth',0.5,'FEdgeColor',frameColour);
axis(ax,'off');
tightmap;
hold(ax,'on');
end

function [x,y] = projectDots(projection,lat,lon)
centre = projection.origin(2);
lon = mod(lon-centre+180,360)-180+centre;
[x,y] = projfwd(projection,lat,lon);
end

function writeModuleMovie(filename,fig,dots,frameTitle,years,ids,x,y,palette, ...
    frameRate,markerSize,unassignedLabel)
fprintf('Writing %d frames to %s\n',numel(years)+1,filename);
writer = VideoWriter(filename,'MPEG-4');
writer.FrameRate = frameRate;
open(writer);
cleanup = onCleanup(@() close(writer));
drawnow
writeVideo(writer,getframe(fig));
for k = 1:numel(years)
    use = ids{k}>0;
    set(dots,'XData',x(use),'YData',y(use),'CData',palette(ids{k}(use),:),'SizeData',markerSize);
    frameTitle.String = sprintf('Year %d | %s',years(k),unassignedLabel);
    drawnow
    writeVideo(writer,getframe(fig));
end
clear cleanup
fprintf('Saved module movie to %s\n',filename);
end
