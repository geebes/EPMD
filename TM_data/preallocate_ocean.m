% Script to generate pre-rolled ocean data files
% These include the transport matrix and associated metadata

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

% clear
% clc
cd ~/GitHub/EPMD/TM_data
addpath(genpath('~/GitHub/EPMD'))

%%
depth_scheme = 'alldepths'; % 'surface_transport', 'surface_depth_integ' or 'alldepths'
specID       = 'GUD_X17';

%%

load(['~/GitHub/EPMD/GUD_forcing/' specID '_abundance.mat']);

if strmatch(depth_scheme,'alldepths')
    TM_filename = ['pre-rolled/' specID '_weighted_transport.mat'];
elseif strmatch(depth_scheme,'surface_depth_integ')
    TM_filename = ['pre-rolled/' specID '_surface_transport.mat'];
    
elseif strmatch(depth_scheme,'surface_transport')
    TM_filename = ['pre-rolled/surface_transport.mat'];
end

disp(['Preparing ' TM_filename])

[ocean] = initialise_ocean(depth_scheme,abundance);

save(TM_filename,'ocean')
%%

function [ocean] = initialise_ocean(depth_scheme,cell_conc)

    % Initialise ocean metadata
    disp('Initialising ocean metadata')

    % load transport matrix data
    TM_boxes        = load('~/TM_data/MITgcm_ECCO_v4/Matrix13/Data/boxes.mat');
    TM_matrices     = load('~/TM_data/MITgcm_ECCO_v4/Matrix13/TMs/matrix_nocorrection_annualmean.mat','Aexpms','Aimpms');
    TM_grid         = load('~/TM_data/MITgcm_ECCO_v4/grid.mat');
    TM_basin_mask   = load('~/TM_data/MITgcm_ECCO_v4/GCM/basin_mask.mat');
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
    
    % get surface grid parameters
    ocean.lon       = TM_boxes.Xboxnom(ocean.Ib);
    ocean.lat       = TM_boxes.Yboxnom(ocean.Ib);
    ocean.ix        = TM_boxes.ixBox(ocean.Ib);
    ocean.iy        = TM_boxes.iyBox(ocean.Ib);
    ocean.iz        = TM_boxes.izBox(ocean.Ib);
    ocean.volume    = TM_boxes.volb(ocean.Ib);
    ocean.iface     = TM_boxes.boxfacenum(ocean.Ib);
    % temperatures are surface temperatures
    load('~/GitHub/EPMD/GUD_forcing/Theta.mat');
    ocean.theta         = theta(ocean.Ib,:);
    ocean.ann_theta     = mean(theta(ocean.Ib,:),2);
    
    data                = TM_basin_mask.basin_mask;
    ocean.basins        = gcmfaces2vector(data,TM_boxes);
    ocean.basin_names	= TM_basin_mask.basin_names;

    % gather transport matrices (explicit and implicit)
    Aexpms = TM_matrices.Aexpms;
    
%%

    
    if strmatch(depth_scheme,{'surface_transport','surface_depth_integ'})
        disp('Generating surface-only transport matrix')
        disp(' - fluxes originating in interior are set to zero')
        disp(' - fluxes from surface to interior are kept at source (i.e. put on diagonal)')
        disp(' - implicit matrix can be ignored as it contains no horizontal fluxes')

        % convert from concentration flux to mass flux (* sink volume) [nb x nb]
        [isnk,isrc,val]=find(Aexpms);
        val=val.*volb(isnk);
        A_mass=sparse(isnk,isrc,val,size(Aexpms,1),size(Aexpms,2));
        
        % Fluxes within surface layer
        B_mass  = A_mass(ocean.Ib,ocean.Ib);
        % Downwelling fluxes (Ib to Ii) should loop back into surface (i.e. trapped particles)
        Bdn     = A_mass(ocean.Ii,ocean.Ib);
        % sum of downwelling fluxes out of each surface cell (prevented from leaving)
        Fdn     = sum(Bdn,1);
        
        I       = speye(size(B_mass));
        dindx   = find(I); % get index for diagonal
        
        % Generate new surface-only matrix
        divergent=true;
        if divergent % divergent version
            B_mass(dindx) = B_mass(dindx) + Fdn';
        else % mass-conservative version
            B_mass(dindx) = 0 - (sum(B_mass,2) - B_mass(dindx));
        end
       
        % get volume of surface grid boxes to convert back to concentration flux
        ocean.volume=volb(ocean.Ib);
        
        % convert from concentration flux to mass flux (* sink volume) [nb x nb]
        [isnk,isrc,val]=find(B_mass);
        val=val./ocean.volume(isnk);
        B=sparse(isnk,isrc,val,size(B_mass,1),size(B_mass,2));

        if strmatch(depth_scheme,'surface_depth_integ')
            disp('Also save depth-integrated abundances')
            
            % get index for unique [x,y] locations
            [~,ia,ic] = unique([TM_boxes.Xboxnom TM_boxes.Yboxnom],'rows','stable');
            % Create sparse mapping matrix [60646,2406992]
            % This maps all fluxes with a non-surface component to within the surface
            % Vertical fluxes with no horizontal component go on diagonal
            % Any horizontal fluxes are mapped to off diagonal
            % ones in each row correspond to spatial indices for one location
            % e.g. location n has 50 depth levels, corresponding to index n1:n50
            % so row n has ones in columns n1:n50
            M1=sparse(ic,1:size(ic,1),1,size(ia,1),size(ic,1));
            M2=M1';
            % Need to weight layers by fraction of column total biomass
            abundance       = max(cell_conc.*volb,1);   % get cell abundance in grid boxes
            ann_abundance   = mean(abundance,2);        % get annual mean abundance in grid boxes
            % get volume of water columns to convert back to concentration flux
            ocean.abundance     = M1*abundance;         % save water-column abundances
            ocean.ann_abundance = M1*ann_abundance;     % save annual water-column abundances
        end
    elseif strmatch(depth_scheme,'alldepths')
        disp('Generating depth-integrated transport matrix')
        disp(' - all horizontal fluxes are mapped into a single layer')
        disp(' - implicit matrix can be ignored as it contains no horizontal fluxes')
        % get index for unique [x,y] locations
        [~,ia,ic] = unique([TM_boxes.Xboxnom TM_boxes.Yboxnom],'rows','stable');

        % Create sparse mapping matrix [60646,2406992]
        % This maps all fluxes with a non-surface component to within the surface
        % Vertical fluxes with no horizontal component go on diagonal
        % Any horizontal fluxes are mapped to off diagonal
        % ones in each row correspond to spatial indices for one location
        % e.g. location n has 50 depth levels, corresponding to index n1:n50
        % so row n has ones in columns n1:n50
        M1=sparse(ic,1:size(ic,1),1,size(ia,1),size(ic,1));
        M2=M1';
    
        % Need to weight layers by fraction of column total biomass
        abundance       = max(cell_conc.*volb,1);   % get cell abundance in grid boxes
        ann_abundance   = mean(abundance,2);        % get annual mean abundance in grid boxes
        int_abundance   = M1*ann_abundance;
        
        % convert from concentration flux to cell flux (* sink abundance) [nb x nb]
        [isnk,isrc,val]=find(Aexpms);
        val=val.*ann_abundance(isnk); % scaled by volume and cell conc
        A_cell=sparse(isnk,isrc,val,size(Aexpms,1),size(Aexpms,2));
        
        % map from full grid to single layer
        B_cell = M1*A_cell*M2;
        
        % convert from cell flux to mass flux (* sink volume) [nb x nb]
        [isnk,isrc,val]=find(B_cell);
        val=val./int_abundance(isnk);
        B=sparse(isnk,isrc,val,size(B_cell,1),size(B_cell,2));
            
%         % correct mass conservation
%         B = conserve_mass(B,ocean.volume);
        
        % get volume of water columns to convert back to concentration flux 
        ocean.volume        = M1*TM_boxes.volb;
        ocean.abundance     = M1*abundance;         % save water-column abundances
        ocean.ann_abundance = M1*ann_abundance;     % save annual water-column abundances

    end
    
    % correct negatives
    B = no_negatives(B,ocean.volume);
    
    % Timestep
    B = B.*ocean.dt_sec;
    
    % Add identity matrix to get discrete time form
    ocean.B = B + speye(size(B));

end





