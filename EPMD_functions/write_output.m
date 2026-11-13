function cmat = write_output(yr,x,tseries_x,daily_x,run_options,cmat)
% Save one complete annual EPMD time slice.
% Write a temporary MAT-file and rename it after saving.

%% Prepare output filenames

if ~run_options.save_data
    return
end

finalFile = fullfile(cmat.outputDir, ...
    sprintf('year_%06d.mat',yr));

temporaryFile = fullfile(cmat.outputDir, ...
    sprintf('year_%06d.partial.mat',yr));

% Never overwrite a completed annual output.
if isfile(finalFile)
    error('Annual output already exists: %s',finalFile);
end

% Do not reuse an incomplete file from a failed write.
if isfile(temporaryFile)
    delete(temporaryFile);
end

%% Prepare annual population distributions

% Record the completed year.
yrs_saved = yr;

% Preserve the full population distribution for restarting.
% For all experiments except selective_dispersal, x itself
% is already the complete restart state.
if strcmp(run_options.seed_dist,'selective_dispersal')

    x_restart = x;

    % Aggregate phenotypes within each lineage.
    % Exclude the resident population as in the original code.
    xx = reshape(full(x),[],run_options.nphen,run_options.nlineages);
    xx = squeeze(sum(xx(:,:,1:end-1),2));

    x = sparse(xx);

end

%% Prepare first-occupation times

t_occupied = [];

switch run_options.seed_dist

    case {'neutral','selective_dispersal','nonadaptive_dispersal'}

        % Obtain first-occupation times from the simulation.
        tt = reshape(full(run_options.t_occupied), ...
            [],run_options.nphen,run_options.nlineages);

        % Exclude zero times and the resident population.
        tt(tt==0) = NaN;
        tt = squeeze(min(tt(:,:,1:end-1),[],2));

        % Restore zeros for locations never occupied.
        tt(isnan(tt)) = 0;

        t_occupied = sparse(tt);

end

%% Prepare daily outputs

% Preserve the original rule for time-series saving.
switch run_options.seed_dist

    case {'selective_dispersal','nonadaptive_dispersal'}
        tseries_x = [];

end

% Save full daily distributions only in the final year.
if ~run_options.save_daily_data || yr~=run_options.nyear
    daily_x = [];
end

%% Write annual MAT-file

fprintf('Saving year %d to %s\n',yr,temporaryFile);

% Selective dispersal requires a separate full restart state.
if strcmp(run_options.seed_dist,'selective_dispersal')

    save(temporaryFile, ...
        'yr','yrs_saved','x','x_restart', ...
        't_occupied','tseries_x','daily_x','-v7.3');

else

    save(temporaryFile, ...
        'yr','yrs_saved','x', ...
        't_occupied','tseries_x','daily_x','-v7.3');

end

% Publish the completed file with an atomic rename.
[success,message] = movefile(temporaryFile,finalFile);

if ~success
    error('Unable to finalise year %d: %s',yr,message);
end

fprintf('Completed year %d: %s\n',yr,finalFile);

end