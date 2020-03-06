function [] = load_EPMD_output(fname)

    grid_load('~/GitHub/EPMD/nctiles_grid/',5,'nctiles')
    gcmfaces_global
    load coastlines
    land = shaperead('landareas', 'UseGeoCoords', true);