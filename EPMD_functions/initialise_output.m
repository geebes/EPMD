function [cmat,run_options] = initialise_output(run_options,ocean)
% Initialise Matfile output files...
% 
%  Syntax:
%    [cmat,run_options] = initialise_output(run_options,ocean)
%
%  (run_options is structural array defined in run_EPMD)
%  (ocean is structural array defined in allocate_ocean)

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

%%
    % Coordinates of JGOFS sites for saving time series output
    tseries_lon = [-064 -158 -145 +062 -140 -180 -019 +068 -170 ];
    tseries_lat = [+032 +023 +050 +016 +000 -076 +047 -051 -61.5];
    run_options.n_tseries_loc=numel(tseries_lon);  
    
    for i=1:run_options.n_tseries_loc
        % find distances from JGOFS sites to all ocean grid coordinates
        distances = sum(([tseries_lon(i) tseries_lat(i)]-[ocean.lon ocean.lat]).^2,2);
        % get index of grid coordinates closest to JGOFS sites
        [~,run_options.site_indices(i)] = min(distances);
    end
    
    % generate output filename
    if run_options.save_data
        if strmatch(run_options.TM_scheme,'surface_transport')    
            fname = [run_options.seed_dist '_' ...
                     run_options.trajectory '_' ...
                     run_options.annual_cycle '_GUD_' ...
                     run_options.DARWIN_pop '_surface_transport'...
                     run_options.suffix '.mat'];
        else
            fname = [run_options.seed_dist '_' ...
                     run_options.trajectory '_' ...
                     run_options.annual_cycle '_' ...
                     run_options.TM_scheme ...
                     run_options.suffix '.mat'];
        end

        % create Matfile object
        warning('off')
        filename = ['~/GitHub/EPMD/Output/' fname];
        delete(filename);
        warning('on')
        cmat  = matfile(filename,'Writable', true);

        % Forcing (Climatological)
        cmat.Temperature        =ocean.forcing_temp;
        cmat.CarryingCapacity	=ocean.forcing_PCapacity;
        cmat.T_optima           =run_options.T_opt;

        % yearly snapshots (1st Jan) 
        cmat.x(run_options.nyear,1)        ={[]}; 
        cmat.x_restart                     =[];      
        cmat.t_occupied                    =[];

        % Timeseries (Daily)
        cmat.tseries_x(run_options.nyear,1)      ={[]};  
        cmat.tseries_lon =tseries_lon;
        cmat.tseries_lat =tseries_lat;
        cmat.tseries_ind =run_options.site_indices;
        cmat.run_options =run_options;
        cmat.ocean       =ocean;        
        cmat.yrs_saved   = 0;
    else
        cmat=[];
    end
end