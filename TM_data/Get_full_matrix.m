clear
% clc
cd ~/GitHub/EPMD/TM_data
addpath(genpath('~/GitHub/EPMD'))

%%
depth_scheme = 'alldepths'; % 'surface_transport' or 'alldepths'
specID       = 'GUD_X17';

%%


% Initialise ocean metadata
disp('Initialising ocean metadata')

% load transport matrix data
TM_boxes        = load('~/Transport_Matrices/MITgcm_ECCO_v4/Matrix13/Data/boxes.mat');
TM_matrices     = load('~/Transport_Matrices/MITgcm_ECCO_v4/Matrix13/TMs/matrix_nocorrection_annualmean.mat','Aexpms');
TM_grid         = load('~/Transport_Matrices/MITgcm_ECCO_v4/grid.mat');
TM_basin_mask   = load('~/Transport_Matrices/MITgcm_ECCO_v4/GCM/basin_mask.mat');
gcmfaces_global
grid_load('../nctiles_grid/',5,'nctiles');

coastlines      = load('coastlines');
ocean.land      = shaperead('landareas','UseGeoCoords',true);
ocean.dt_sec    = 60*60*6; % 6 hours
ocean.nsubtime  = 4;
ocean.dt_rat    = ocean.dt_sec/TM_grid.deltaT;

% extract grid metadata
ocean.Ib        = find(TM_boxes.izBox==1); % surface boundary
ocean.Ii        = find(TM_boxes.izBox~=1); % interior
nb              = numel(TM_boxes.izBox);
I               = speye(nb,nb);
volb            = TM_boxes.volb;
% extract TM data
Aexpms          = TM_matrices.Aexpms;

% get surface grid parameters
ocean.lon       = TM_boxes.Xboxnom;
ocean.lat       = TM_boxes.Yboxnom;
ocean.z         = TM_boxes.Zboxnom;
ocean.ix        = TM_boxes.ixBox;
ocean.iy        = TM_boxes.iyBox;
ocean.iz        = TM_boxes.izBox;
ocean.volume    = TM_boxes.volb;
ocean.iface     = TM_boxes.boxfacenum;
% temperatures are surface temperatures
load('~/GitHub/EPMD/GUD_forcing/Theta.mat');
ocean.theta         = theta;
ocean.ann_theta     = mean(theta,2);

data                = TM_basin_mask.basin_mask;
ocean.basins        = gcmfaces2vector(data,TM_boxes);
ocean.basin_names	= TM_basin_mask.basin_names;

% Pre-process TMs
disp('Correcting TM for mass conservation')
A_conc = no_negatives(Aexpms,volb);

ocean.A = A_conc;


save('full_matrix.mat','ocean', '-v7.3');
















