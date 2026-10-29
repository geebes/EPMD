%% EPMD neutral connectivity and homogenisation
% Load EPMD output and calculate connectivity and homogenisation.
% The final population represents the widespread background population.
% Exclude it from both connectivity and homogenisation calculations.

clearvars

%% Settings

outputFile = '/Users/baw1d17/GitHub/EPMD/Output/neutral_stochastic_static_GUD_X01_surface_transport.mat';

nSample = 1000;              % Number of native locations to compare
detectionThreshold = 0;      % Minimum relative abundance for detection
rng(12345);                  % Reproducible sampling

%% Load output

fprintf('Loading %s\n',outputFile);

load(outputFile,'x','x_restart','yrs_saved');

% Check that population time slices are available.
if ~exist('x','var') || isempty(x)
    error('No population time slices found in output file.');
end

% Count completed years, ignoring empty time slices.
nTimes = find(~cellfun('isempty',x),1,'last');

% Number of native ECCO locations.
nLocations = size(x_restart,1);

% Number of explicitly tracked founding lineages.
nLineages = size(x_restart,2)-1;

fprintf('Loaded %d time slices, %d locations and %d tracked lineages.\n', ...
    nTimes,nLocations,nLineages);

%% Select locations

% Sample native ECCO locations without replacement.
sampleLocations = randperm(nLocations,min(nSample,nLocations));

% Number of unordered pairs of sampled locations.
nPairs = numel(sampleLocations)*(numel(sampleLocations)-1)/2;

% Precompute upper-triangular pair indices.
upperPairs = triu(true(numel(sampleLocations)),1);

%% Calculate diagnostics

years = (1:nTimes)';
sharedConnectivity = zeros(nTimes,1);
homogenisation = zeros(nTimes,1);

for k = 1:nTimes

    % Extract the abundance matrix at this time slice.
    % Rows = locations; columns = lineages.
    abundance = full(x{k}(sampleLocations,1:nLineages));

    % Normalise by total abundance, INCLUDING background population.
    % This preserves the actual relative frequency of each tracked lineage.
    allPopulations = full(x{k}(sampleLocations,:));
    total = sum(allPopulations,2);
    total(total==0) = 1;

    relative = abundance./total;

    % Apply detection threshold.
    present = relative>detectionThreshold;

    % Fraction of paired locations sharing at least one lineage.
    overlap = double(present)*double(present)';

    sharedConnectivity(k) = nnz(overlap(upperPairs)>0)/nPairs;

    % Homogenisation: mean pairwise Bray-Curtis similarity.
    % Renormalise tracked lineages, excluding the background population.
    relative(~present) = 0;
    trackedTotal = sum(relative,2);
    occupied = trackedTotal>0;

    if nnz(occupied)>1

        community = relative(occupied,:)./trackedTotal(occupied);

        pairSimilarity = 1-0.5*pdist(community,'cityblock');

        % Empty communities contribute zero similarity.
        homogenisation(k) = sum(pairSimilarity)/nPairs;

    end

    fprintf('Year %3d: connectivity %.5f, homogenisation %.5f\n', ...
        years(k),sharedConnectivity(k),homogenisation(k));

end

%% Plot

figure(1); clf
set(gcf,'Color','w','Name','EPMD neutral connectivity and homogenisation')

tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

ax1 = nexttile;
plot([0;years],[0;sharedConnectivity],'o-','LineWidth',1.5)
ylabel('Connectivity')
ylim([0 1])
grid on
title(sprintf('Fraction of paired locations sharing at least one lineage > %g', ...
    detectionThreshold))

ax2 = nexttile;
plot([0;years],[0;homogenisation],'o-','LineWidth',1.5)
ylabel('Homogenisation')
xlabel('Time (years)')
ylim([0 1])
grid on
title('Mean fraction of relative community composition shared between paired locations')

linkaxes([ax1 ax2],'x')