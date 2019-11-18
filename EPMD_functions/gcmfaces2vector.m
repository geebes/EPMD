function vctr = gcmfaces2vector(data,TM_boxes)

facenum = TM_boxes.boxfacenum(TM_boxes.izBox==1);

vctr=[];
for f=1:5
    % index face cells at surface
    indx = find(TM_boxes.izBoxFace{f}==1);
    % extract non NaNs
    vctr = [vctr; data{f}(~isnan(data{f}))];
    
end


