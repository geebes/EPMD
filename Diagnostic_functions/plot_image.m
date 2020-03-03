function [ax] = plot_image(x,mygrid,ocean)


% ax=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
ax=axesm ('eqaazim','frame','on','FlineWidth',0.5,'Origin',[-46.283333,-86.666665]);

ax.Position=ax.Position+[-0.1 -0.2 2 4].*0.01;

for i=1:3
    fld=vector2gcmfaces(x(:,i),ocean.iface,ocean.ix,ocean.iy,ocean.iz);
    [lon lat Xmap(:,:,i)]=convert2pcol(mygrid.XC,mygrid.YC,fld);
end


Xmap(Xmap==0)=NaN;

sh=surfacem(lat, lon, Xmap);

axis off

end