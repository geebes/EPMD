function write_diagnostic_output(yr,time,homogenisation,bcMap, ...
    detectionThreshold,cmat)
% Publish one completed diagnostic file for this model year.

bcMapMetric = 'area_weighted_mean_pairwise_native';
homogenisationMetric = bcMapMetric;

% Write atomically so readers see only completed files.
temporaryFile = fullfile(cmat.outputDir, ...
    sprintf('diagnostics_%06d.partial.mat',yr));

finalFile = fullfile(cmat.outputDir, ...
    sprintf('diagnostics_%06d.mat',yr));

save(temporaryFile, ...
    'time','homogenisation','bcMap','bcMapMetric', ...
    'homogenisationMetric','detectionThreshold','-v7.3');

[success,message] = movefile(temporaryFile,finalFile);

if ~success
    error('Unable to finalise diagnostics: %s',message);
end
end
