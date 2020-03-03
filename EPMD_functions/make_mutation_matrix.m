function [run_options] = make_mutation_matrix(run_options);
% Initialises mutation matrix for adaptive simulations
%
%  Syntax:
%    [run_options] = make_mutation_matrix(run_options)
%
%  (run_options is structural array defined in run_EPMD)

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

    % generate vector of thermal optima for n phenotypes
    run_options.T_opt = linspace(-2,36,run_options.nphen);
    % get non-structured copy
    T_opt=reshape(run_options.T_opt',1,[]);

    % Generate Mutation Matrix
    % calculate delta temperature between phenotypes
    run_options.delta_Topt = range(T_opt')./(numel(unique(T_opt))-1);
    % calculate mutation rate
    run_options.delta_m = run_options.sigma_m.^2 ./ (3.*run_options.delta_Topt.^2);

    % find diagonals in mutation matrix
    d1=find(abs(T_opt-T_opt')==0); 
    % find index of adjacent phenotypes
    d2=find(abs(T_opt-T_opt')==run_options.delta_Topt); 

    % allocate sparse mutation matrix
    mutmat=spalloc(run_options.nphen,run_options.nphen,3*run_options.nphen);
    % off-diagonals set to mutation rate
    mutmat(d2)=mutmat(d2)+run_options.delta_m;
    % diagonals set to 1 - sum of rows
    mutmat(d1)=0-sum(mutmat,1);

    % add identity matrix for mutmat_I
    mutmat_I=mutmat+speye(size(mutmat));

    % save in structural array
    run_options.mutmat_I = mutmat_I;
    run_options.mutmat   = mutmat;

end