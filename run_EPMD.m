% Script to run EPMD

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>


clear
clc

cd ~/GitHub/EPMD
addpath(genpath('~/GitHub/EPMD'))
addpath EPMD_functions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set options
run_options.TM_scheme       = 'GUD_X01_weighted_transport'; % 'surface_transport', 'GUD_X01_surface_transport' or 'GUD_X01_weighted_transport', or similar
run_options.seed_dist       = 'selective_dispersal';    % 'preadapted', 'equal', 'lineages', 'neutral', 'selective_dispersal'
run_options.trajectory      = 'stochastic'; % 'stochastic' or 'deterministic'
run_options.annual_cycle    = 'static';     % 'static' or 'seasonal'
run_options.seedseed        = 2;            % seed for global seeding sites

% optional string to add to output filenames
run_options.suffix          = '_try';

% note depth-integrated TMs are each associated with a particular abundance distribution
% (this will be overwritten if using depth-integrated biomass-weighted transport)
run_options.DARWIN_pop      = 'X01'; 

run_options.save_data       = true;

run_options.nyear           = 200;          
run_options.nday        	= 365;          

% N.B. not applicable for neutral model
run_options.nphen           = 77;           % number of phenotypes
run_options.w               = 6;           	% Niche breadth
run_options.sigma_m         = 0.1;          % Mutation size

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run EPMD_spmd
EPMD_spmd(run_options)