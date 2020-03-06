function [ocean,run_options,t_occ,x,i_lastyr] = load_EPMD_output(input_filename,pathname)
%load_EPMD_output Load EPMD .mat output file
    
    % load MatObj file
    matObj  = matfile([pathname '/' input_filename '.mat']);
    
    % Load data structures
    ocean       = matObj.ocean;
    run_options = matObj.run_options;
    
    % load abundance data
    x = matObj.x;
    
    % load dispersal times
    t_occupied  = matObj.t_occupied;
    t_occ       = full(t_occupied(:,:));
    t_occ(~t_occ(:)) = NaN;
    
    % identify last year of simulation
    i_lastyr    = matObj.yrs_saved;
    
    disp('::::::::::::::::::::::::::::::::::')
    disp(input_filename)
    disp([num2str(i_lastyr) ' years evaluated.'])
    disp(['ocean is ' num2str(100.*nnz(t_occ)/numel(t_occ),'%2.0f') '% connected.'])
    disp('::::::::::::::::::::::::::::::::::')
    
end