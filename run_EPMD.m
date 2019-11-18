clear
clc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load preallocated surface transport matrix and associated metadata
load('TM_data/surface_transport.mat');
addpath EPMD_functions
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set options

run_options.solver          = 'parallel';	% serial or parallel
run_options.seed_dist       = 'neutral';    % 'preadapted', 'equal', 'singleancestor', 'neutral', 'mutation_dispersal'
run_options.trajectory      = 'stochastic'; % 'stochastic' or 'deterministic'
run_options.annual_cycle    = 'static';     % 'static' or 'seasonal'
run_options.DARWIN_pop      = 'X_01';

run_options.mutation        = false;      	% allow mutations
run_options.selection       = false;
run_options.save_data       = true;

run_options.nyear           = 100;          
run_options.nday        	= 365;          
run_options.parp_size       = 18;           % n processors to use in parallel
run_options.nlevel          = 1;            % number of levels inhabite

% N.B. not applicable for neutral model
run_options.nancestral      = 1;           % number of distinct ancestral lineages
run_options.npopn           = 77;           % number of phenotypes
run_options.w               = 6;           	% Niche breadth
run_options.sigma_m         = 0.1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate evenly-spaced sample points
[ocean] = sample_grid(ocean,3,'quad');

% generate mutation matrix
[run_options] = make_mutation_matrix(run_options);

% initialise the global metacommunity
[x,run_options,ocean] = seed_metacommunity(ocean,run_options);

% use SPMD if comparing multiple independent lineages
if run_options.nancestral>1
    setup_SPMD(x,ocean,run_options)
end

[cmat,run_options] = initialise_output(run_options,ocean);


%%

for yr=1:run_options.nyear
    tseries_x  = zeros(run_options.n_tseries_loc,run_options.nxc,run_options.nday);

    

        
    for dy = 1:run_options.nday
        tday=tic;

        % DAILY RESOLVED CARRYING CAPACITY
        N=ocean.forcing_PCapacity(:,dy);


        if run_options.selection
            % CALCULATE SELECTION COEFFICIENT (s) AS FUNCTION OF TEMPERATURE AND PREFERENCE (s<=1)
            s = exp(-((ocean.forcing_temp(:,dy)-run_options.T_opt)./run_options.w).^2);  % Seasonal Temperature limitation
        else
            s = 1;
        end

        % GENERATION time steps in each day
        for t=1:(3600/(ocean.dt_sec*4))*24 

            if run_options.mutation
                % TRAIT DIFFUSION
                dxdt = x * run_options.mutmat;   % Redistribute mutants
                x = x + dxdt;
            end

            % PHYSICAL TRANSPORT
            for st=1:ocean.nsubtime   % nsubtime transport timesteps per generation
                x=ocean.B*x;          % calculate probability of each population in each grid cell
            end

            % SELECTION (abundance and fitness weighted or just abundance weighted)
            % calculate abundance and fitness weighted probability of
            % selection in next generation, normalising so sum(x)=1
            
            p = s.*x ./ sum(s.*x,2);

            

            [loc,trc,p_l]=find(p);  % work only with locations/genotypes where nonzero probability of cells
            N_l=N(loc);             % carrying capacities for extant populations

            % STOCHASTIC OR DETERMINISTIC POPULATION DYNAMICS
            switch run_options.trajectory
                case 'stochastic'
                    mu_x   =N_l.*p_l;
                    sigma_x=sqrt(N_l.*p_l.*(1-p_l));
                    X=normrnd(mu_x,sigma_x); % sample population
                    % multinomial (~normal)
                case 'deterministic'
                    X=N_l.*p_l;
            end

            % Set abundance to positive integer value
            X=floor(X);
            X(X<0)=0;

            % Calculate as fraction of local carrying capacity
            x = X./N_l;

            % convert x and X from vector to matrix form
            x=sparse(loc,trc,x,run_options.nxr,run_options.nxc);
            X=sparse(loc,trc,X,run_options.nxr,run_options.nxc);

%             % Check to see if metacommunity matrix is actually sparse
%             if nnz(x)./numel(x) > 0.25
%                 disp('full!')
%                 x=full(x);
%                 X=full(X);
%             end

        end

        switch run_options.seed_dist
            case 'neutral'
                isoccupied = find(x & run_options.t_occupied==0);
                date=(yr-1 + dy./run_options.nday);
                run_options.t_occupied(isoccupied) = date;
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COLLATE TIMESERIES OUTPUT (DAILY) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if run_options.save_data & ~run_options.resident
            tseries_x(:,:,dy) = x(run_options.site_indices,:);
        end


        disp(['Year ' num2str(yr,'%03i') ', Day ' num2str(dy,'%03i') ' (' num2str(toc(tday)) ' seconds).']);


    end % end day loops
    
    disp('----------------------------------')
    [cmat] = write_output(yr,x,tseries_x,run_options,cmat);
    disp('----------------------------------')
end % end year loop

%%
