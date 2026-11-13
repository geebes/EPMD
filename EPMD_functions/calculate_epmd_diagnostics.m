function [H,BC] = calculate_epmd_diagnostics( ...
    x,localIDs,nSeedPopulations,detectionThreshold,areaWeights)
% Calculate Bray-Curtis diagnostics from population blocks distributed
% across SPMD workers, without reconstructing the global matrix.

% Identify seeded populations in this worker's column block.
seedColumns = localIDs<=nSeedPopulations;

% Calculate total abundance including the resident.
totalNative = gplus(full(sum(x,2)));

% Extract seeded populations from this worker.
seededNative = full(x(:,seedColumns));

% Apply the same detection threshold.
seededNative( ...
    seededNative<=detectionThreshold.*totalNative) = 0;

% Calculate total seeded abundance across all workers.
seededTotalNative = gplus(sum(seededNative,2));

% Identify cells containing seeded populations.
occupiedNative = seededTotalNative>0;

totalDistance = [];
if nnz(occupiedNative)>1
    % Each worker contributes distances for its own population columns.
    localDistance = native_pairwise_distance_sum( ...
        seededNative,seededTotalNative,areaWeights);
    totalDistance = gplus(localDistance);
end
[similarity,H] = native_pairwise_similarity(totalDistance,occupiedNative,areaWeights);
BC = single(similarity);
end
