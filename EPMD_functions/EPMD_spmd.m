% Run EPMD_spmd
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate mutation matrix
[run_options] = make_mutation_matrix(run_options);

% Initialise ocean structural aaray
[ocean] = allocate_ocean(run_options,3,'quad');

% initialise the global metacommunity
[x,run_options,ocean] = seed_metacommunity(run_options,ocean);

% initialise output files
[cmat,run_options] = initialise_output(run_options,ocean);

% setup Serial or Parallel 
switch run_options.solver
    case 'parallel'
        myCluster               = parcluster('local');
        run_options.parp_size   = myCluster.NumWorkers;
        run('setup_SPMD'); % (Single Program Multiple Data)
    case 'serial'
        run_options.parp_size   = 1;
        run('setup_SPMD'); % (Single Program Multiple Data)
        nxc                     = run_options.nxc;
        indx                    = 1:nxc;
end

%%

for yr=1:run_options.nyear

    spmd(run_options.parp_size)
        tseries_x  = zeros(run_options.n_tseries_loc,nxc,run_options.nday);

        for dy = 1:run_options.nday
            tday=tic;
            % DAILY-RESOLVED CARRYING CAPACITY
            N=ocean.forcing_PCapacity(:,dy);

            if run_options.selection
                % CALCULATE SELECTION COEFFICIENT (s) AS FUNCTION OF TEMPERATURE AND PREFERENCE (s<=1)
                s = exp(-((ocean.forcing_temp(:,dy)-T_opt)./run_options.w).^2);  % Seasonal Temperature limitation
            else
                s = 1;
            end
            
            % GENERATION time steps in each day
            for t=1:(3600/(ocean.dt_sec*4))*24 

                if run_options.mutation
                    % TRAIT DIFFUSION
                    dxdt = x * mutmat;   % Redistribute mutants
                    x    = x + dxdt;
                end
                
                % PHYSICAL TRANSPORT
                for st=1:ocean.nsubtime   % nsubtime transport timesteps per generation
                    x=ocean.B*x;          % calculate probability of each population in each grid cell
                end

                % SELECTION (abundance and fitness weighted or just abundance weighted)
                % calculate abundance and fitness weighted probability of
                % selection in next generation, normalising so sum(x)=1
                popn_selected = s.*x;
                comn_selected = sum(popn_selected,2);
                global_comn_selected = gplus(comn_selected);
                
                p = popn_selected ./ global_comn_selected;
                
                ind         = find(p & N);  % work only with locations/genotypes where nonzero probability of cells
                [loc,trc]   = ind2sub(size(N),ind);
                p_l         = full(p(ind));
                N_l         = full(N(loc));       % carrying capacities for extant populations
                
                % STOCHASTIC OR DETERMINISTIC POPULATION DYNAMICS
                switch run_options.trajectory
                    case 'stochastic'
                        mu_x   =N_l.*p_l;
                        sigma_x=sqrt(N_l.*p_l.*(1-p_l));
                        X_l=normrnd(mu_x,sigma_x); % sample population
                        % Set abundance to integer value
                        X_l=floor(X_l);
                    case 'deterministic'
                        X_l=N_l.*p_l;
                end 

                % Set abundance to positive value
                X_l(X_l<0)=0;

                % Calculate as fraction of local carrying capacity
                x_l = X_l./N_l;
  
                % convert x and X from vector to matrix form
                x=sparse(loc,trc,x_l,run_options.nxr,nxc);
                X=sparse(loc,trc,X_l,run_options.nxr,nxc);
            end

            switch run_options.seed_dist
                case 'neutral'
                    occdate=(yr-1 + dy./run_options.nday);
                    isoccupied = find(x & t_occ==0);
                    t_occ(isoccupied) = occdate;
            end
%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % COLLATE TIMESERIES OUTPUT (DAILY) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if run_options.save_data & ~run_options.resident
                tser_x(:,:,dy) = full(x(run_options.site_indices,:));
            end
            
            if labindex==1
                disp(['Year ' num2str(yr,'%03i') ', Day ' num2str(dy,'%03i') ' (' num2str(toc(tday)) ' seconds).']);
            end

        end % end day loops
   
        % put local data back in global codistributed array
        xD(:,indx)           = x;
        t_occD(:,indx)       = t_occ;
        tseries_xD(:,indx,:) = tser_x;
        
    end  % end spmd block (exiting to write data)
    
    % Gather x data
    xG                       = gather(xD);
    tseries_xG               = gather(tseries_xD);
    run_options.t_occupied   = gather(t_occD);
    
    disp('----------------------------------')
    if run_options.save_data
        [cmat] = write_output(yr,xG,tseries_xG,run_options,cmat);
        disp('----------------------------------')
    end
end % end year loop

%%
