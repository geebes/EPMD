function [ax] = plot_vector(x,cscale,mygrid,ocean)


ax=axesm ('mollweid','frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);

ax.Position=ax.Position+[-0.1 -0.2 2 4].*0.01;

fld=vector2gcmfaces(x,ocean.iface,ocean.ix,ocean.iy,ocean.iz);

[lon lat Xmap]=convert2pcol(mygrid.XC,mygrid.YC,fld);

Xmap(Xmap==0)=NaN;

if cscale=='log'
    Xmap=log10(Xmap);
end

sh=surfacem(lat, lon, Xmap);

axis off

end