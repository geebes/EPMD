function [cmat] = write_output(yr,x,tseries_x,run_options,cmat)

    if run_options.save_data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE GLOBAL SNAPSHOT DATA (YEARLY ON 31 DEC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Saving time-slice ...')
        cmat.yrs_saved         = yr;
        
        switch run_options.seed_dist
            case 'selective_dispersal'
                % Do not write x for huge selective_dispersal runs
            otherwise
                cmat.x(yr,1)           = {[x]};
        end
        
        switch run_options.seed_dist
            case {'neutral','selective_dispersal'}
                % Write dispersal times
                cmat.t_occupied   = run_options.t_occupied;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE DAILY TIMESERIES DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch run_options.seed_dist
            case 'selective_dispersal'
                % Do not write x time-series for huge selective_dispersal runs
            otherwise
                disp('Saving time-series ...')
                cmat.tseries_x(yr,1)       = {[tseries_x]};
        end
    end

end