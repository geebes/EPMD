function [cmat,run_options] = initialise_output(run_options,ocean)
% Initialise EPMD output directory and metadata.
% Schedule diagnostics and full population snapshots separately.

%% Identify time-series locations

tseries_lon = [-064 -158 -145 +062 -140 -180 -019 +068 -170 ];
tseries_lat = [+032 +023 +050 +016 +000 -076 +047 -051 -61.5];

run_options.n_tseries_loc = numel(tseries_lon);

for i = 1:run_options.n_tseries_loc

    % Find the nearest ocean grid cell.
    distances = sum( ...
        ([tseries_lon(i) tseries_lat(i)]-[ocean.lon ocean.lat]).^2,2);

    [~,run_options.site_indices(i)] = min(distances);

end

%% Initialise output

if ~run_options.save_data
    cmat = [];
    return
end

% Preserve the existing experiment naming convention.
if run_options.dispersal && strcmp(run_options.TM_scheme,'surface_transport')
    transport_scheme = ['GUD_' run_options.DARWIN_pop '_surface_transport'];
elseif run_options.dispersal
    transport_scheme = run_options.TM_scheme;
else
    transport_scheme = ['GUD_' run_options.DARWIN_pop '_no_dispersal'];
end

runName = [run_options.seed_dist '_' ...
           run_options.trajectory '_' ...
           run_options.annual_cycle '_' ...
           transport_scheme ...
           run_options.suffix];

% Create a separate directory for each experiment.
outputDir = fullfile(getenv('HOME'),'GitHub','EPMD','Output', ...
    [runName '_years']);

% Prevent accidental overwriting of an earlier simulation.
if isfolder(outputDir)
    error('Output directory already exists: %s',outputDir);
end

[success,message] = mkdir(outputDir);

if ~success
    error('Cannot create output directory: %s',message);
end

% Retain metadata needed for the initial output.
cmat.outputDir = outputDir;
cmat.tseries_lon = tseries_lon;
cmat.tseries_lat = tseries_lat;
cmat.tseries_ind = run_options.site_indices;

fprintf('Output directory: %s\n',outputDir);

%% Schedule output
% Diagnostics can be saved on any day. Full snapshots are saved at year end.
% Round to the available resolution and save the final day and year.
lastYear = run_options.nyear;
nday = run_options.nday;
cmat.saveYears = [];
cmat.diagnosticDaysByYear = cell(lastYear,1);
if lastYear>0
    lastDay = lastYear*nday;

    diagnosticSteps = unique(round(logspace(0,log10(lastDay), ...
        1+ceil(run_options.diagnostics_per_decade*log10(lastDay)))));
    cmat.saveYears = unique(round(logspace(0,log10(lastYear), ...
        1+ceil(run_options.snapshots_per_decade*log10(lastYear)))));

    for step = diagnosticSteps
        yearIndex = floor((step-1)/nday)+1;
        dayInYear = mod(step-1,nday)+1;
        cmat.diagnosticDaysByYear{yearIndex}(end+1) = dayInYear;
    end
    cmat.diagnosticDaysByYear{1} = [0 cmat.diagnosticDaysByYear{1}];
end

end
