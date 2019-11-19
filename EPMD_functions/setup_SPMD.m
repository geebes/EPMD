
% open parallel pool
open_parallel(run_options.parp_size);

parp_size       = run_options.parp_size;
npopn           = run_options.npopn;
nlineages       = run_options.nlineages;
nxc             = run_options.nxc;
nxr             = run_options.nxr;
T_opt           = run_options.T_opt;
mutmat          = run_options.mutmat;
nday            = run_options.nday;
n_tseries_loc   = run_options.n_tseries_loc;
t_occupied      = run_options.t_occupied;

ii=round(linspace(0,npopn,parp_size+1)).*nlineages; % divide populations among n proc
idiv = sort(diff(ii),'descend');       % get size of blocks on each proc
ii=[0 cumsum(idiv)];                   % sort ii as for idiv
 
spmd(run_options.parp_size)
    % create index to assign 'x' to workers
    codistr = codistributor1d(2, idiv, [nxr nxc]);

    % assign 'x' to workers according to codistr
    xD=codistributed(x,codistr);
    XD=xD.*0; % also initialise X
    t_occD = XD; % and t_occupied

    T_optD=codistributed(T_opt,codistr); % distribute

    % replicate mutation array along diagonal according to size of block
    if run_options.mutation
        mutmat=kron(eye(idiv(labindex)./npopn),mutmat); % order is important (identity first)
    end
    
    % indexing needs to be done manually 
    % https://uk.mathworks.com/help/parallel-computing/codistributed.colon.html
    indx=(ii(labindex)+1):(ii(labindex+1));
    
    tseries_xD  = codistributed(zeros(n_tseries_loc,nxc,nday),codistr);
    
    % Initialise local arrays
    x      = getLocalPart(xD);
    X      = getLocalPart(XD);
    T_opt  = getLocalPart(T_optD);
    tser_x = getLocalPart(tseries_xD);
    t_occ  = getLocalPart(t_occD);
    nxc    = size(x,2);
    
end

who