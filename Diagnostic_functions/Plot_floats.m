clf
ax=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
colormap(flipud(gray))

geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]); % Very SLOW!!!!!

A=ocean.B;

B=speye(size(TM));

B=B(:,ocean.sample_points);

n=10;
B=repmat(B,n,1);
B=reshape(B,60646,94*n);


tmax=1e5;

indx=zeros(tmax,size(B,2));

[~,I]=max(B,[],1);
indx(1,:)=I;
lat=ocean.lat(indx(1,:));
lon=ocean.lon(indx(1,:));
h=[];
scatterm(lat,lon,10,'r','filled');


for t=2:4*10*365
    B2=(A*B).*rand(size(B));
    [~,I]=max(B2,[],1);
    
    indx(t,:)=I;
    B=sparse(I,1:numel(I),1,size(B,1),size(B,2));

end
    



lat=ocean.lat(indx(1:t,:));
lon=ocean.lon(indx(1:t,:));

delete(h)
h=plotm(lat,lon,'k','LineW',1);


title(t)

cd = [uint8(flipud(gray(t))*255)].';
cd(4,:)=fliplr(cd(3,:)); % opacity
for i=1:numel(h)
    set(h(i).Edge, 'ColorBinding','interpolated', 'ColorData',cd)
end
