function [ax,surfaces] = plot_vector(x,cscale,mygrid,ocean,proj)
% plot ocean surface vector ('x') as a 2D map with projection 'proj'
% 
% plot_vector(x,cscale,mygrid,ocean,proj)
%
% cscale = 'lin' or 'log'
% mygrid and ocean are predefined structural arrays

persistent cachedKey cachedMesh

ax=axesm (proj,'frame','on','FlineWidth',0.5,'MapLonLimit',[-245 115],'MapLatLimit',[-90 90]);
% ax=axesm ('eqaazim','frame','on','FlineWidth',0.5,'Origin',[-46.283333,-86.666665]);

ax.Position=ax.Position+[-0.1 -0.2 2 4].*0.01;

fld = vector2gcmfaces(x,ocean.iface,ocean.ix,ocean.iy,ocean.iz);
surfaces = gobjects(mygrid.nFaces,1);
projection = getm(ax);
isMollweide = strcmp(projection.mapprojection,'mollweid');
meshKey = {mygrid.dirGrid,mygrid.nFaces,projection};
reuseMesh = isequaln(meshKey,cachedKey);
if ~reuseMesh
    mesh = cell(mygrid.nFaces,1);
end
meshProperties = {'XData','YData','ZData','UserData','ButtonDownFcn'};

for face = 1:mygrid.nFaces
    values = fld{face};
    [ni,nj] = size(values);
    if strcmp(cscale,'log')
        values(values==0) = NaN;
        values = log10(values);
    end
    colours = nan(ni+1,nj+1);
    colours(1:ni,1:nj) = values;
    if reuseMesh
        geometry = cachedMesh{face};
        colours(geometry.seamCells) = NaN;
        meshArgs = [meshProperties; reshape(geometry.properties,1,[])];
        surfaces(face) = surface(meshArgs{:},'Parent',ax, ...
            'CData',colours,'FaceColor','flat','EdgeColor','none');
    else
        [lon,lat] = readCorners(mygrid.dirGrid,face,ni+1,nj+1);
        geometry.seamCells = [];
        if isMollweide
            % Retain all projected vertices and hide only cells that wrap.
            centre = projection.origin(2);
            wrappedLon = mod(lon-centre+180,360)-180;
            cornerLon = cat(3,wrappedLon(1:ni,1:nj), ...
                wrappedLon(2:end,1:nj),wrappedLon(2:end,2:end), ...
                wrappedLon(1:ni,2:end));
            [rows,cols] = find(max(cornerLon,[],3)-min(cornerLon,[],3)>180);
            geometry.seamCells = sub2ind(size(colours),rows,cols);
            colours(geometry.seamCells) = NaN;
            % Keep exact seam vertices on the same side as the mask.
            wrappedLon(wrappedLon==-180) = -180+1e-8;
            [mapX,mapY] = projfwd(projection,lat,wrappedLon+centre);
            surfaces(face) = surface('Parent',ax,'XData',mapX,'YData',mapY, ...
                'ZData',zeros(size(mapX)),'CData',colours, ...
                'FaceColor','flat','EdgeColor','none');
        else
            surfaces(face) = surfacem(lat,lon,colours,'FaceColor','flat');
        end
        geometry.properties = get(surfaces(face),meshProperties);
        mesh{face} = geometry;
    end
end
if ~reuseMesh
    cachedKey = meshKey;
    cachedMesh = mesh;
end

% Set the view from the projected mesh, independent of the map frame and overlays.
xCoordinates = get(surfaces,'XData');
yCoordinates = get(surfaces,'YData');
xMin = min(cellfun(@(v) min(v(:),[],'omitnan'),xCoordinates),[],'omitnan');
xMax = max(cellfun(@(v) max(v(:),[],'omitnan'),xCoordinates),[],'omitnan');
yMin = min(cellfun(@(v) min(v(:),[],'omitnan'),yCoordinates),[],'omitnan');
yMax = max(cellfun(@(v) max(v(:),[],'omitnan'),yCoordinates),[],'omitnan');
xlim(ax,[xMin xMax]);
ylim(ax,[yMin yMax]);
axis(ax,'off');
end

function [lon,lat] = readCorners(gridDir,face,ni,nj)
% MITgrid records 6 and 7 contain the full XG and YG corner arrays.
filename = fullfile(gridDir,sprintf('tile%03d.mitgrid',face));
[file,message] = fopen(filename,'r','ieee-be');
assert(file>=0,'Cannot open %s: %s',filename,message);
cleanup = onCleanup(@() fclose(file));
nVertices = ni*nj;
assert(fseek(file,5*nVertices*8,'bof')==0,'Cannot seek to corners in %s.',filename);
[corners,count] = fread(file,[nVertices 2],'double=>double');
assert(count==2*nVertices,'Incomplete corner data in %s.',filename);
lon = reshape(corners(:,1),ni,nj);
lat = reshape(corners(:,2),ni,nj);
end
