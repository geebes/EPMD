function run_options = default_output_options(run_options)
% Fill output and diagnostic options for both scripts and direct callers.

defaults.diagnostics_per_decade = 30;
defaults.snapshots_per_decade = 9;
defaults.diagnostic_detection_threshold = 0;

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(run_options,name)
        run_options.(name) = defaults.(name);
    end
end
end
