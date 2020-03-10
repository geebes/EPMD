function D2 = bray_curtis(ZI,ZJ)
% Calculate Bray-Curtis Dissimilarity
%
% bray_curtis takes ZI and ZJ as inputs and returns the Bray-Curtis dissimilarity matrix.
% ZI is a 1-by-n vector containing a single row from the metacommunity
% matrix (representing 1 local community).
% 
% ZJ is the full metacommunity matrix (representing all sampled communities)
%
% Bray-Curtis returns D2, which is an m-by-1 vector of dissimilarities. 
% The jth element of D2 is the distance between the community ZI and community ZJ(j,:).

D2=zeros(size(ZI'));

C  = sum(min(ZI,ZJ),2); % sum of minimum count of each species between sites
SS = sum(ZI+ZJ,2);      % sum of each species across both sites

D2 = 1-2.*C./SS;