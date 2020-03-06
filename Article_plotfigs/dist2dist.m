function D2 = dist2dist(ZI,ZJ)

i=size(ZI,2)-size(ZJ,1)+1;
D2=ZI(i:end)';

