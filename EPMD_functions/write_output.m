function [cmat] = write_output(yr,x,tseries_x,run_options,cmat)

    if run_options.save_data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE GLOBAL SNAPSHOT DATA (YEARLY ON 31 DEC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Saving time-slice ...')
        cmat.yrs_saved         = yr;
        cmat.x(yr,1)           = {[x]};
        
        switch run_options.seed_dist
            case 'neutral'
                cmat.t_occupied   = run_options.t_occupied;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE DAILY TIMESERIES DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~run_options.resident
            disp('Saving time-series ...')
            cmat.tseries_x(yr,1)       = {[tseries_x]};
        end
    end

end