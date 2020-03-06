function [coastlat coastlon mygrid land] = load_geo_grid(nctiles_path)
%LOAD_GEO_GRID Load ocean grid data
%   
    grid_load(nctiles_path,5,'nctiles')
    gcmfaces_global
    load coastlines
    land = shaperead('landareas', 'UseGeoCoords', true);

end

