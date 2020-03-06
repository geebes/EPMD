function D2 = bray_curtis(ZI,ZJ)

D2=zeros(size(ZI'));

C  = sum(min(ZI,ZJ),2); % sum of minimum count of each species between sites
SS = sum(ZI+ZJ,2);      % sum of each species across both sites

D2 = 1-2.*C./SS;