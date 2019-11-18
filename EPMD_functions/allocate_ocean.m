function [ocean] = allocate_ocean(run_options,varargin)

load(['TM_data/pre-rolled/' run_options.TM_scheme '.mat']);

% generate equally spaced points on surface of a sphere

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

% convert to spherical coordinates
[thetaR,phi,~]=cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));

% convert to degrees
xs=rad2deg(thetaR);
xs(xs<0)=xs(xs<0)+360;
ys=rad2deg(phi);

lon = ocean.lon;
lat = ocean.lat;

% find sample points nearest to each model gridpoint
Kall = dsearchn([xs ys],[ocean.lon ocean.lat]);

% remove duplicates
[K,~,Kall] = unique(Kall);

% xs(K) & ys(K) are grid coordinates nearest to sample grid
% find index of those coordinates in model grid
K = dsearchn([ocean.lon ocean.lat],[xs(K),ys(K)]);

ocean.sample_points = K;