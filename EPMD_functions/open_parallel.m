function open_parallel(parp_size)
% Initialise parallel pool...
% 
%  Syntax:
%    open_parallel(parp_size)
%
%  (parp_size is number of processors to use)

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

    if nargin==0
        myCluster = parcluster('local');
        parp_size=myCluster.NumWorkers;
    end


    p = gcp('nocreate');
    if isempty(p)
        parpool(parp_size);
    elseif p.NumWorkers~=parp_size
        delete(gcp);
        parpool(parp_size);
    end

end