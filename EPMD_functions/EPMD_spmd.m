% Run EPMD_spmd
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate mutation matrix
disp('Create mutation matrix for lineages')
[run_options] = make_mutation_matrix(run_options);

% Initialise ocean structural array
disp('Allocate ocean metadata and generate seeding points')
[ocean] = allocate_ocean(run_options,3,'quad');

% initialise the global metacommunity
disp('Set up initial seed populations')
[x,run_options,ocean] = seed_metacommunity(run_options,ocean);

% initialise output files
disp('Initialise output files')
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
site_indices        = run_options.site_indices;

forcing_PCapacity   = ocean.forcing_PCapacity;
forcing_temp        = ocean.forcing_temp;
dt_sec              = ocean.dt_sec;
nsubtime            = ocean.nsubtime;
B                   = ocean.B;

for yr=1:run_options.nyear

    disp('Opening SPMD block')
    spmd(run_options.parp_size)
        tser_x = zeros(n_tseries_loc,nxc,nday);
        
        for dy = 1:nday
            tday=tic;
            % DAILY-RESOLVED CARRYING CAPACITY
            N=forcing_PCapacity(:,dy);

            if selection
                % CALCULATE SELECTION COEFFICIENT (s) AS FUNCTION OF TEMPERATURE AND PREFERENCE (s<=1)
                s = exp(-((forcing_temp(:,dy)-T_opt)./w).^2);  % Seasonal Temperature limitation
            else
                s = 1;
            end
            
            % GENERATION time steps in each day
            for t=1:(3600/(dt_sec*4))*24 

                if mutation
                    % TRAIT DIFFUSION
                    dxdt = x * mutmat;   % Redistribute mutants
                    x    = x + dxdt;
                end
                
                % PHYSICAL TRANSPORT
                for st=1:nsubtime   % nsubtime transport timesteps per generation
                    x=B*x;          % calculate probability of each population in each grid cell
                end

                % SELECTION (abundance and fitness weighted or just abundance weighted)
                % calculate abundance and fitness weighted probability of
                % selection in next generation, normalising so sum(x)=1
                popn_selected = s.*x;
                comn_selected = sum(popn_selected,2);
                global_comn_selected = gplus(comn_selected);
                
                p = popn_selected ./ global_comn_selected;
                
                if run_options.xsparse
                    ind         = find(p & N);  % work only with locations/genotypes where nonzero probability of cells
                    [loc,trc]   = ind2sub(size(N),ind);
                    p_l         = full(p(ind));
                    N_l         = full(N(loc));       % carrying capacities for extant populations

                    % STOCHASTIC OR DETERMINISTIC POPULATION DYNAMICS
                    switch trajectory
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
                    x=sparse(loc,trc,x_l,nxr,nxc);
                    X=sparse(loc,trc,X_l,nxr,nxc);
                else
                    mu_x    = N.*p;
                    sigma_x = sqrt(N.*p.*(1-p));
                    X_l     = normrnd(mu_x,sigma_x); % sample population
                    % Set abundance to integer value
                    X = floor(X);
                     % Set abundance to positive value
                    X(X<0)=0;
                    % Calculate as fraction of local carrying capacity
                    x=X./N;
                end
            end

            switch seed_dist
                case {'neutral','selective_dispersal'}
                    % get current date
                    occdate=(yr-1 + dy./nday);
                    % find all pppulations and locations occupied for the first time
                    % (i.e. occupied but no date of first occupation)
                    [i,j] = find(x & t_occ==0);
                    % assign current date to newly occupied 
                    t_occ = t_occ + sparse(i,j,occdate,size(t_occ,1),size(t_occ,2));
            end
%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % COLLATE TIMESERIES OUTPUT (DAILY) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if save_data 
                tser_x(:,:,dy) = full(x(site_indices,:));
            end
            
            if labindex==1
                disp(['Year ' num2str(yr,'%03i') ', Day ' num2str(dy,'%03i') ' (' num2str(toc(tday)) ' seconds).']);
            end

        end % end day loops
   
        % put local data back in global codistributed array
        xD(:,indx)           = x;
        tseries_xD(:,indx,:) = tser_x;
        t_occD(:,indx)       = t_occ;
        
    end  % end spmd block (exiting to write data)
    disp('----------------------------------')
    disp('Closed SPMD block to write output data')

    %%
    % Gather x data
    xG                       = gather(xD);
    tseries_xG               = gather(tseries_xD);
    run_options.t_occupied   = gather(t_occD);
    
    if run_options.save_data
        [cmat] = write_output(yr,xG,tseries_xG,run_options,cmat);
        disp('----------------------------------')
    end
end % end year loop

%%
