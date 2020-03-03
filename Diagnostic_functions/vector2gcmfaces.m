function fld=vector2gcmfaces(data,boxfacnum,ix,iy,iz)

fld=gcmfaces(5);

tmp{1} = zeros(90,270).*NaN;
tmp{2} = zeros(90,270).*NaN;
tmp{3} = zeros(90,90).*NaN;
tmp{4} = zeros(270,90).*NaN;
tmp{5} = zeros(270,90).*NaN;

for f=1:5
    ii=find(boxfacnum==f & iz==1);
    
    ind = sub2ind(size(tmp{f}),ix(ii),iy(ii));
    
    tmp{f}(ind) = data(ii);
    
    eval(['fld.f' num2str(f) '=tmp{f};']);
    
end