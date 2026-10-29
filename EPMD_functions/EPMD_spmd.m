function EPMD_spmd(run_options)
% Wrapper function for EPMD
%
%  Syntax:
%    EPMD_spmd(run_options) 
%
%  (run_options is structural array defined in run_EPMD)

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate mutation matrix
disp('Create mutation matrix for lineages')
[run_options] = make_mutation_matrix(run_options);

% Initialise ocean structural array
disp('Allocate ocean metadata and generate seeding points')
[ocean] = allocate_ocean(run_options,run_options.seedseed,'quad');

% initialise the global metacommunity
disp('Set up initial seed populations')
[x,run_options,ocean] = seed_metacommunity(run_options,ocean);

% initialise output files
disp('Initialise output files')
[cmat,run_options] = initialise_output(run_options,ocean);

% Set the number of parallel workers.
switch run_options.solver
    case 'parallel'
        pool = gcp('nocreate');
        if isempty(pool)
            run_options.parp_size = parcluster('local').NumWorkers;
        else
            run_options.parp_size = pool.NumWorkers;
        end
    case 'serial'
        run_options.parp_size = 1;
end

% Distribute lineage blocks without codistributed arrays.
[x,T_opt,mutmat,indx,nxc,t_occ,run_options] = ...
    setup_SPMD(run_options,x);

%% Extract data from structural arrays for use inside SPMD block
nday                = run_options.nday;
n_tseries_loc       = run_options.n_tseries_loc;
selection           = run_options.selection;
w                   = run_options.w;
mutation            = run_options.mutation;
nxr                 = run_options.nxr;
trajectory          = run_options.trajectory;
seed_dist           = run_options.seed_dist;
save_data           = run_options.save_data;
save_daily_data     = run_options.save_daily_data;
site_indices        = run_options.site_indices;

forcing_PCapacity   = ocean.forcing_PCapacity;
forcing_temp        = ocean.forcing_temp;
dt_sec              = ocean.dt_sec;
nsubtime            = ocean.nsubtime;
B                   = ocean.B;
%%
for yr=1:run_options.nyear
    forcing_temp        = ocean.forcing_temp ...
                        + yr.*run_options.warming_rate;

    disp('Opening SPMD block')
    tyear = tic;
    spmd(run_options.parp_size)
        tser_x  = zeros(n_tseries_loc,nxc,nday);
        if save_daily_data 
            daily_x = zeros(nxr,nxc,nday);
        end
        
        for dy = 1:nday
            tday=tic;
            % DAILY-RESOLVED CARRYING CAPACITY
            N=forcing_PCapacity(:,dy);

            if selection
                % CALCULATE SELECTION COEFFICIENT (s) AS FUNCTION OF TEMPERATURE AND PREFERENCE (s<=1)
                s = exp(-((mean(forcing_temp(:,dy),2)-T_opt)./w).^2);  % Seasonal Temperature limitation
            else
                s = 1;
            end
            
            % GENERATION time steps in each day
            for t=1:(3600/(dt_sec*4))*24 

                if mutation
                    % TRAIT DIFFUSION
                    dxdt = (x) * mutmat;   % Redistribute mutants
                    x    = x + dxdt;
                end

                % PHYSICAL TRANSPORT
                for st=1:nsubtime   % nsubtime transport timesteps per generation
                    x=B*x;          % calculate probability of each population in each grid cell
                end
                
                % set negative s to zero
                x = max(x,0);
                
                % SELECTION (abundance and fitness weighted or just abundance weighted)
                % calculate abundance and fitness weighted probability of
                % selection in next generation, normalising so sum(x)=1                
                popn_selected = s.*x;
                comn_selected = sum(popn_selected,2);
                global_comn_selected = gplus(comn_selected);
                
                p = popn_selected ./ global_comn_selected;
                p(global_comn_selected<=0,:)=0; % check divide by zero
                
                mu_x    = N.*p;   % mean = unselected + selected part
                sigma_x = sqrt(N.*p.*(1-p)); % variance based only on selected part
                X       = normrnd(mu_x,sigma_x); % sample population
                    
                % Set abundance to integer value
                X = floor(X);
                % Set abundance to positive value
                X(X<0)=0;
                % Calculate as fraction of local carrying capacity
                x= X./N;
                % crude check for div 0
                x(isnan(x))=0;
            end

            switch seed_dist
                case {'neutral','selective_dispersal','nonadaptive_dispersal'}
                    
                    % get current date
                    occdate=(yr-1 + dy./nday);
                    % find all pppulations and locations occupied for the first time
                    % (i.e. occupied but no date of first occupation)
                    [i,j] = find(x & t_occ==0);
                    % assign current date to newly occupied 
                    t_occ = t_occ + sparse(i,j,occdate,size(t_occ,1),size(t_occ,2));
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % COLLATE TIMESERIES OUTPUT (DAILY) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if save_data 
                tser_x(:,:,dy) = full(x(site_indices,:));
            end
            if save_daily_data && yr==run_options.nyear
                daily_x(:,:,dy) = full(x);
            end
            
            if labindex==1
                disp(['Year ' num2str(yr,'%03i') ', Day ' num2str(dy,'%03i') ' (' num2str(toc(tday)) ' seconds).']);
            end

        end % end day loops

        
    end  % end spmd block (exiting to write data)
    disp('----------------------------------')
    disp('Closed SPMD block to write output data')
    disp(['Year ' num2str(yr) ' evaluated in ' num2str(toc(tyear)) ' seconds.'])

    % Reconstruct global arrays from worker-owned lineage blocks.
    xG = zeros(nxr,run_options.nxc);

    tseries_xG = zeros(n_tseries_loc,run_options.nxc,nday);

    run_options.t_occupied = zeros(nxr,run_options.nxc);

    if save_daily_data && yr==run_options.nyear
        daily_xG = zeros(nxr,run_options.nxc,nday);
    else
        daily_xG = [];
    end

    % Collect each worker's populations.
    for worker = 1:run_options.parp_size

        ids = indx{worker};

        xG(:,ids) = x{worker};

        run_options.t_occupied(:,ids) = t_occ{worker};

        if save_data
            tseries_xG(:,ids,:) = tser_x{worker};
        end

        if save_daily_data && yr==run_options.nyear
            daily_xG(:,ids,:) = daily_x{worker};
        end

    end
    
    if run_options.save_data
        [cmat] = write_output(yr,xG,tseries_xG,daily_xG,run_options,cmat);
        disp('----------------------------------')
    end
    
    disp('Clear large output arrays')
    clear xG tseries_xG daily_xG
    run_options.t_occupied = [];
end % end year loop

%%
