function [similarity,homogenisation] = native_pairwise_similarity(distanceSum,occupied,areaWeights)
% Convert area-weighted L1 distance sums to focal-cell and global similarity.

similarity = zeros(size(areaWeights));
homogenisation = 0;
if nnz(occupied)<2
    return
end

occupiedWeights = areaWeights(occupied);
comparisonWeight = sum(areaWeights)-areaWeights;
% Exclude the focal cell and give pairs involving empty cells zero similarity.
similarity(occupied) = (sum(occupiedWeights)-occupiedWeights-0.5*distanceSum) ...
    ./comparisonWeight(occupied);
similarity = min(max(similarity,0),1);

% Each distinct pair has weight proportional to the product of its areas.
pairWeights = areaWeights.*comparisonWeight;
homogenisation = sum(pairWeights.*similarity)/sum(pairWeights);
end
