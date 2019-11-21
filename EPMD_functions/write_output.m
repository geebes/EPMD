function [cmat] = write_output(yr,x,tseries_x,run_options,cmat)

    if run_options.save_data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE GLOBAL SNAPSHOT DATA (YEARLY ON 31 DEC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Saving time-slice ...')
        cmat.yrs_saved         = yr;
        
        switch run_options.seed_dist
            case 'selective_dispersal'
                % Integrate cell numbers across all phenotypes in each lineage
                % (N.B. exclude resident)
                xx=reshape(full(x),[],run_options.nphen,run_options.nlineages);
                xx=squeeze(sum(xx(:,:,1:end-1),2));
                xx=sparse(xx);
            otherwise
                xx=x;
        end
        cmat.x(yr,1)           = {[xx]};
        
        switch run_options.seed_dist
            case {'neutral','selective_dispersal'}
                % find minimum occupancy time across all phenotypes in each lineage
                % (N.B. exclude resident)
                t_occ = run_options.t_occupied;
                % reshape into 3D array
                tt=reshape(full(t_occ),[],run_options.nphen,run_options.nlineages);
                % set zeros to NaN
                tt(tt==0)=NaN;
                % find minimum time for each lineage/location (excluding resident)
                tt=squeeze(min(tt(:,:,1:end-1),[],2));
                % put back zeros
                tt(isnan(tt))=0;
                % Write dispersal times
                cmat.t_occupied   = sparse(tt);
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