function distanceSum = native_pairwise_distance_sum(abundance,totals,areaWeights)
% Sum area-weighted L1 distances from each occupied cell across these population columns.
% Totals must include all seeded populations, including other workers.
occupied = find(totals>0);
occupiedWeights = areaWeights(occupied);
totalOccupiedWeight = sum(occupiedWeights);
inverseTotals = 1./totals(occupied);
distanceCorrection = zeros(size(occupied));
distanceFromZeros = 0;
for j = 1:size(abundance,2)
    [rows,~,values] = find(abundance(occupied,j));
    [values,order] = sort(values.*inverseTotals(rows));
    rows = rows(order);

    weights = occupiedWeights(rows);
    weightedValues = weights.*values;

    % Start at the distance from zero, then correct the positive rows.
    distanceFromZeros = distanceFromZeros + sum(weightedValues);
    zeroWeight = totalOccupiedWeight-sum(weights);
    prefixWeight = zeroWeight+cumsum(weights);
    distanceCorrection(rows) = distanceCorrection(rows) + ...
        (2*prefixWeight-totalOccupiedWeight).*values-2*cumsum(weightedValues);
end
distanceSum = distanceFromZeros + distanceCorrection;
end
