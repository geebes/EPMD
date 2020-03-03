function [ocean] = allocate_ocean(run_options,varargin)
% Load pre-allocated ocean grid and generate seed locations
%
%  Examples:
%    [ocean] = allocate_ocean(run_options) 
%  generates seed locations with default n=3 subdivisions on 'quad' grid
%
%    [ocean] = allocate_ocean(run_options,n,mesh) 
%  generates seed locations with n subdivisions on 'mesh' grid
%
%  (run_options is structural array defined in run_EPMD)

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

%% Load pre-allocated ocean grid
load(['TM_data/pre-rolled/' run_options.TM_scheme '.mat']);

%% generate seed locations using SubdivideSphericalMesh
% (equally spaced points on surface of a sphere)
if nargin == 1
    n    = 3;
elseif nargin==2
    if isnumeric(varargin{1})
        mesh = 'quad';
        n    = varargin{1};
    else
        mesh = varargin{1};
        n    = 3;
    end
elseif nargin==3
    if isnumeric(varargin{1})
        n    = varargin{1};
        mesh = varargin{2};
    else
        mesh = varargin{1};
        n    = varargin{2};
    end
end

switch mesh
    case 'tri'
        TR=SubdivideSphericalMesh(IcosahedronMesh,n);
        xyz=TR.X;
    case 'quad'
        fv=SubdivideSphericalMesh(QuadCubeMesh,n);
        xyz=fv.vertices;
end

%%

% convert to spherical coordinates
[thetaR,phi,~]=cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));

% convert to degrees
xs=rad2deg(thetaR);
xs(xs<0)=xs(xs<0)+360;
ys=rad2deg(phi);

% find sample points nearest to each model gridpoint
Kall = dsearchn([xs ys],[ocean.lon ocean.lat]);

% remove duplicates
[K,~,Kall] = unique(Kall);

% xs(K) & ys(K) are grid coordinates nearest to sample grid
% find index of those coordinates in model grid
K = dsearchn([ocean.lon ocean.lat],[xs(K),ys(K)]);

% place in structural array
ocean.sample_points = K;

disp(['Placed ' num2str(numel(K)) ' seed sites'])