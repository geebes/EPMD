clear
clc
cd ~/GitHub/EPMD
addpath(genpath('~/GitHub/EPMD'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath EPMD_functions
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set options
run_options.TM_scheme       = 'GUD_X01_weighted_transport'; % 'surface_transport' or 'GUD_X01_weighted_transport', or similar
run_options.seed_dist       = 'lineages';    % 'preadapted', 'equal', 'lineages', 'neutral'
run_options.trajectory      = 'stochastic'; % 'stochastic' or 'deterministic'
run_options.annual_cycle    = 'static';     % 'static' or 'seasonal'

% note depth-integrated TMs are each associated with a particular abundance distribution
% (this will be overwritten if using depth-integrated biomass-weighted transport)
run_options.DARWIN_pop      = 'X01'; 

run_options.save_data       = true;

run_options.nyear           = 100;          
run_options.nday        	= 365;          

% N.B. not applicable for neutral model
run_options.npopn           = 77;           % number of phenotypes
run_options.w               = 6;           	% Niche breadth
run_options.sigma_m         = 0.1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('../EPMD_functions/EPMD_spmd')