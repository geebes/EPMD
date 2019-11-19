
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
            ocean.forcing_temp        = ...
                interp1(0:12,[ocean.theta           ocean.theta(:,1)          ]',linspace(0,12,run_options.nday))';
            ocean.forcing_PCapacity   = ...
                interp1(0:12,[run_options.PCapacity run_options.PCapacity(:,1)]',linspace(0,12,run_options.nday))';
    end
    
    T_opt       = run_options.T_opt;
    delta_Topt  = run_options.delta_Topt;
    npopn       = run_options.npopn;
    nancestral  = run_options.nancestral;
    K           = ocean.sample_points;
    B           = ocean.B;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set initial abundance distributions
    switch run_options.seed_dist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'preadapted' % set one completely dominant type preadapted to each location
            % bin  temperature according to T_opt range
            [~,~,bin] = histcounts(temp,[-inf T_opt(1:end-1)+delta_Topt/2 inf]);
            % initialise EiE matrix
            x=zeros(length(B),npopn.*nancestral);
            % find matching seed locations in EiE matrix
            iseed=sub2ind(size(x),ocean.Ib,bin);
            % set type abundance of best adapted type to 1 in each location
            x(iseed)=1;
            x=sparse(x);
            
            run_options.resident  = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'equal'
            % share initial abundance among all types equally
            x=ones(length(B),npopn)./npopn;
            
            run_options.resident  = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'neutral'

            npopn=size(ocean.sample_points,1);
            
            x=sparse(ocean.sample_points,1:npopn,1,length(B),npopn);
            % add global resident population
            x(:,end+1)=1;
            x(ocean.sample_points,end)=0;
            
            run_options.npopn     = npopn+1;
            run_options.T_opt     = zeros(1,run_options.npopn);
            run_options.selection = false;
            run_options.mutation  = false;
            run_options.rel_s     = 1;
            
            run_options.resident  = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         case 'lineages'
%             T_opt2=repmat(T_opt,1,nxc); % replicate T_opt for ancestral species
%             T_opt2=T_opt2(:)'; % reshape to row vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    run_options.nxr =size(x,1); % n grid cells
    run_options.nxc =size(x,2); % n phenotypes
    
    % initialise array for connectivity times
    run_options.t_occupied=zeros(size(x));

end











