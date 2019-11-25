
function [x,run_options,ocean] = seed_metacommunity(run_options,ocean)

    % Load carrying capacity
    if strmatch(run_options.TM_scheme,'surface_transport')
        % Can be selected if using surface transport
        load(['./GUD_forcing/GUD_' run_options.DARWIN_pop '_abundance.mat']);
        run_options.PCapacity = abundance(ocean.Ib,:) .* ocean.volume;
    else
        % Predefined if using depth-integrated biomass-weighted transport
        run_options.PCapacity = ocean.abundance;
    end
    
    % initialise populations adapted to annual temperature
    temp=ocean.ann_theta;
        
    switch run_options.annual_cycle
        case 'static'
            ocean.forcing_temp       = repmat(mean(ocean.theta,2),1,run_options.nday);
            ocean.forcing_PCapacity  = repmat(mean(run_options.PCapacity,2),1,run_options.nday);
        case 'seasonal'
            ocean.forcing_temp        = interp1(0:12,[ocean.theta ocean.theta(:,1)]',linspace(0,12,run_options.nday))';
            ocean.forcing_PCapacity   = interp1(0:12,[run_options.PCapacity run_options.PCapacity(:,1)]',linspace(0,12,run_options.nday))';
    end
    
    T_opt       = run_options.T_opt;
    delta_Topt  = run_options.delta_Topt;
    nphen       = run_options.nphen;
    K           = ocean.sample_points;
    B           = ocean.B;
    
    run_options.nlineages	= 1;
    run_options.t_occupied  = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set initial abundance distributions
    switch run_options.seed_dist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'equal'
            % share initial abundance among all types equally
            x=ones(length(B),nphen)./nphen;
            
            run_options.solver    = 'serial';	% serial or parallel
            run_options.mutation  = true;	
            run_options.selection = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'preadapted' % set one completely dominant type preadapted to each location
            % bin  temperature according to T_opt range
            [~,~,bin] = histcounts(temp,[-inf T_opt(1:end-1)+delta_Topt/2 inf]);
            % initialise EiE matrix
            x=zeros(length(B),nphen);
            % find matching seed locations in EiE matrix
            iseed=sub2ind(size(x),ocean.Ib,bin);
            % set type abundance of best adapted type to 1 in each location
            x(iseed)=1;
            x=sparse(x);
            
            run_options.solver    = 'serial';	% serial or parallel
            run_options.mutation  = true;	
            run_options.selection = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'neutral'

            nlineages=size(ocean.sample_points,1);
            
            x=sparse(ocean.sample_points,1:nlineages,1,length(B),nlineages);
            % add global resident population
            x(:,end+1)=1;
            x(ocean.sample_points,end)=0;
            
            run_options.nlineages = nlineages+1;
            run_options.nphen     = 1;
            run_options.T_opt     = zeros(1,run_options.nlineages);
            run_options.selection = false;
            run_options.mutation  = false;
            run_options.rel_s     = 1;
            run_options.solver    = 'parallel';	% serial or parallel
            
            % initialise array for connectivity times
            run_options.t_occupied=zeros(size(x));
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'neutral_lineages' % set locally addapted phenotype as x=1 in each location
            % bin  temperature according to T_opt range
            [~,~,bin] = histcounts(temp,[-inf T_opt(1:end-1)+delta_Topt/2 inf]);
            % initialise EiE matrix
            x1=sparse(length(B),nphen);
            % find matching seed locations in EiE matrix
            iseed=sub2ind(size(x1),ocean.Ib,bin);
            % set type abundance of best adapted type to 1 in each location
            x1(iseed)=1;
            
            nlineages = nphen;
            
            ilin    = (1:nlineages);
            linindx = (ilin-1).*nphen + ilin;
            
            x=sparse(length(B),nphen*nlineages);
            x(:,linindx) = x1;
            
            run_options.T_opt     = repmat(T_opt,1,nlineages); % replicate T_opt for ancestral species
            run_options.nlineages = nlineages;
            run_options.mutation  = true;
            run_options.selection = false;
            run_options.solver    = 'parallel';	% serial or parallel           
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'lineages' % set locally addapted phenotype as x=1 in each location
            % bin  temperature according to T_opt range
            [~,~,bin] = histcounts(temp,[-inf T_opt(1:end-1)+delta_Topt/2 inf]);
            % initialise EiE matrix
            x1=sparse(length(B),nphen);
            % find matching seed locations in EiE matrix
            iseed=sub2ind(size(x1),ocean.Ib,bin);
            % set type abundance of best adapted type to 1 in each location
            x1(iseed)=1;
            
            nlineages = nphen;
            
            ilin    = (1:nlineages);
            linindx = (ilin-1).*nphen + ilin;
            
            x=sparse(length(B),nphen*nlineages);
            x(:,linindx) = x1;
            
            run_options.T_opt     = repmat(T_opt,1,nlineages); % replicate T_opt for ancestral species
            run_options.nlineages = nlineages;
            run_options.mutation  = true;
            run_options.selection = true;
            run_options.solver    = 'parallel';	% serial or parallel
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'selective_dispersal'

            nlineages   = size(ocean.sample_points,1);
            
            x1=sparse(ocean.sample_points,1:nlineages,1,length(B),nlineages);
            
            % find temperature of water at seed locations
            Tsample=ocean.ann_theta(ocean.sample_points);
            diffT=abs(Tsample-run_options.T_opt);
            [~,isample]=min(diffT,[],2);
            
            ilin    = (1:nlineages);
            linindx = (ilin-1).*nphen + isample';
            
            x=sparse(length(B),nlineages*nphen);
            x(:,linindx) = x1;
            
            % add global resident population
            % find best adapted phenotype for each location
            diffT=abs(ocean.ann_theta-run_options.T_opt);
            [~,isample]=min(diffT,[],2);
            % add nphen columns for resident
            x(:,end+1:end+nphen)=0;
            % get index of best adapted phenotype in each location
            resind = sub2ind(size(x),(1:size(x,1))',nlineages*nphen+isample);
            % set resident to 1
            x(resind)=1;
            % reset sample sites to zero
            x(resind(ocean.sample_points))=0; 
                        
            nlineages = nlineages + 1;
            
            run_options.T_opt     = repmat(T_opt,1,nlineages); % replicate T_opt for ancestral species
            run_options.nlineages = nlineages;
            run_options.npopn     = nphen;
            run_options.selection = true;
            run_options.mutation  = true;
            run_options.rel_s     = 1;
            run_options.solver    = 'parallel';	% serial or parallel
            
            % initialise array for connectivity times
            run_options.t_occupied=sparse(size(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    run_options.nxr =size(x,1); % n grid cells
    run_options.nxc =size(x,2); % n phenotypes
    

end











