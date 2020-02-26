function D2 = jaccard(ZI,ZJ)


D2=zeros(size(ZI'));

n_intersect = sum(ZI>0&ZJ>0,2);
n_union     = sum(ZI>0|ZJ>0,2);

D2 = 1 - n_intersect ./ n_union;




