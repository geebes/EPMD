function [run_options] = make_mutation_matrix(run_options);

    run_options.T_opt = linspace(-2,36,run_options.npopn);

    % Mutation matrix

    T_opt=reshape(run_options.T_opt',1,[]);

    run_options.delta_Topt = range(T_opt')./(numel(unique(T_opt))-1);
    run_options.delta_m = run_options.sigma_m.^2 ./ (3.*run_options.delta_Topt.^2);

    d1=find(abs(T_opt-T_opt')==0); % find diagonals
    d2=find(abs(T_opt-T_opt')==run_options.delta_Topt); % find adjacents

    mutmat=spalloc(run_options.npopn,run_options.npopn,3*run_options.npopn);
    mutmat(d2)=mutmat(d2)+run_options.delta_m;
    mutmat(d1)=0-sum(mutmat,1);

    mutmat_I=mutmat+speye(size(mutmat));

    run_options.mutmat_I = mutmat_I;
    run_options.mutmat   = mutmat;

end