function write_initial_output(x,ocean,run_options,cmat)
% Save initial conditions and metadata at t=0.

% Initial population distribution.
x0 = sparse(x);

% Preserve the original output metadata.
Temperature = ocean.forcing_temp;
CarryingCapacity = ocean.forcing_PCapacity;
T_optima = run_options.T_opt;

tseries_lon = cmat.tseries_lon;
tseries_lat = cmat.tseries_lat;
tseries_ind = cmat.tseries_ind;

yrs_saved = 0;

% Save to temporary file.
temporaryFile = fullfile(cmat.outputDir,'initial.partial.mat');
finalFile = fullfile(cmat.outputDir,'initial.mat');

save(temporaryFile, ...
    'x0','ocean','run_options', ...
    'Temperature','CarryingCapacity','T_optima', ...
    'tseries_lon','tseries_lat','tseries_ind', ...
    'yrs_saved','-v7.3');

% Publish completed initial conditions.
[success,message] = movefile(temporaryFile,finalFile);

if ~success
    error('Unable to finalise initial conditions: %s',message);
end
end
