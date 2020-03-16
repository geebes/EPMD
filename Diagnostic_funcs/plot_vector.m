function [ax] = plot_vector(x,cscale,mygrid,ocean,proj)
% plot ocean surface vector ('x') as a 2D map with projection 'proj'
% 
% plot_vector(x,cscale,mygrid,ocean,proj)
%
% cscale = 'lin' or 'log'
% mygrid and ocean are predefined structural arrays

ax=axesm (proj,'frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
% ax=axesm ('eqaazim','frame','on','FlineWidth',0.5,'Origin',[-46.283333,-86.666665]);

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