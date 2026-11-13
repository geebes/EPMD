function area = native_ecco_cell_area(ocean)
% Read ECCO tracer-cell areas (RAC) in native population-row order.

gridDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'nctiles_grid');
faceSizes = [90 270; 90 270; 90 90; 270 90; 270 90];
area = zeros(numel(ocean.iface),1);
for face = 1:size(faceSizes,1)
    inFace = ocean.iface(:)==face;
    if any(inFace)
        tileArea = readTileArea(gridDir,face,faceSizes(face,:)+1);
        indices = sub2ind(size(tileArea),ocean.ix(inFace),ocean.iy(inFace));
        area(inFace) = tileArea(indices);
    end
end
assert(all(isfinite(area) & area>0), ...
    'Expected a positive ECCO surface area for every native cell.');
end

function tileArea = readTileArea(gridDir,face,tileSize)
% RAC is record five of the big-endian MITgrid tile.
filename = fullfile(gridDir,sprintf('tile%03d.mitgrid',face));
[file,message] = fopen(filename,'r','ieee-be');
assert(file>=0,'Cannot open %s: %s',filename,message);
cleanup = onCleanup(@() fclose(file));
assert(fseek(file,4*prod(tileSize)*8,'bof')==0, ...
    'Cannot seek to cell areas in %s.',filename);
[tileArea,count] = fread(file,tileSize,'double=>double');
assert(count==prod(tileSize),'Incomplete cell areas in %s.',filename);
end
