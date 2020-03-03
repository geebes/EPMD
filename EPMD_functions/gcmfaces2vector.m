function vctr = gcmfaces2vector(data,TM_boxes)
% Convert transport matrix grid to vector
%
% Syntax:
%    [vector] = gcmfaces2vector(data,TM_boxes)
%
% data and TM_boxes are output from gcmfaces
% (see https://github.com/MITgcm/gcmfaces)

%  Copyright (C) 2020 Ben Ward <b.a.ward@soton.ac.uk>

facenum = TM_boxes.boxfacenum(TM_boxes.izBox==1);

vctr=[];
for f=1:5
    % index face cells at surface
    indx = find(TM_boxes.izBoxFace{f}==1);
    % extract non NaNs
    vctr = [vctr; data{f}(~isnan(data{f}))];
    
end


