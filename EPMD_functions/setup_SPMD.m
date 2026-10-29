function [x,T_opt,mutmat,indx,nxc,t_occ,run_options] = setup_SPMD(run_options,xInitial)
% Initialise SPMD variables...
% 
%  Syntax:
%    [...] = setup_SPMD(run_options,x)
%
%  (run_options is structural array defined in run_EPMD)
%  (output is full list of SPMD variables)

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>


disp('Setting up SPMD variables')

% Do not use more processsors than model has independent lineages
run_options.parp_size = min(run_options.parp_size,run_options.nlineages);

% open parallel pool
open_parallel(run_options.parp_size);

parp_size       = run_options.parp_size;
nphen           = run_options.nphen;
nlineages       = run_options.nlineages;
nxc             = run_options.nxc;
nxr             = run_options.nxr;
T_opt           = run_options.T_opt;
mutmat1         = run_options.mutmat;
nday            = run_options.nday;
n_tseries_loc   = run_options.n_tseries_loc;
t_occupied      = run_options.t_occupied;

ii=round(linspace(0,nlineages,parp_size+1)).*nphen; % divide populations among n proc
idiv = sort(diff(ii),'descend');       % get size of blocks on each proc
ii=[0 cumsum(idiv)];                   % sort ii as for idiv

%%

mutation = run_options.mutation;

% Retain the global initial conditions before entering SPMD.
ToptInitial = T_opt;
toccInitial = t_occupied;

spmd(run_options.parp_size)

    % Assign complete lineages to each worker.
    indx = (ii(spmdIndex)+1):ii(spmdIndex+1);

    % Number of populations assigned to this worker.
    nxc = numel(indx);

    % Extract the worker's population arrays.
    x = full(xInitial(:,indx));
    T_opt = full(ToptInitial(:,indx));
    t_occ = full(toccInitial(:,indx));

    % Construct the mutation matrix for local populations.
    if mutation
        nrep = nxc/nphen;
        mutmat = kron(eye(nrep),mutmat1);
    else
        mutmat = [];
    end

end