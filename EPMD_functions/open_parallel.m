function open_parallel(parp_size)

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