clear

load('~/Transport_Matrices/MITgcm_ECCO_v4/Matrix13/Data/boxes.mat');
load('~/Transport_Matrices/MITgcm_ECCO_v4/grid.mat');

Ib = find(izBox==1); % surface boundary

QCARBON =  [2.566093668483939E-012,  7.789193721295913E-012,  2.364354020783690E-011,  7.176827455596964E-011,...
   2.178474622439914E-010,  6.612603841985462E-010,  2.007208581666503E-009,  6.092738029662402E-009,  1.849407034084822E-008,...
   2.178474622439914E-010,  6.612603841985462E-010,  2.007208581666503E-009,  6.092738029662402E-009,  1.849407034084822E-008,...
   2.178474622439914E-010,  6.612603841985462E-010,  2.007208581666503E-009,  6.092738029662402E-009,  1.849407034084822E-008,...
   5.613742723010097E-008,  1.704011436062454E-007,  5.172404788573358E-007,  1.570046463929721E-006,  4.765763701139335E-006,...
   1.446613471441439E-005,  2.007208581666504E-009,  6.092738029662406E-009,  1.849407034084824E-008,  5.613742723010101E-008,...
   1.704011436062455E-007,  5.172404788573361E-007,  1.570046463929721E-006,  4.765763701139337E-006,  1.446613471441440E-005,...
   4.391091684330804E-005,  2.007208581666504E-009,  6.092738029662406E-009,  1.849407034084824E-008,  5.613742723010101E-008,...
   1.704011436062455E-007,  5.172404788573361E-007,  1.570046463929721E-006,  4.765763701139337E-006,  1.446613471441440E-005,...
   4.391091684330804E-005,  1.332884461596117E-004,  4.045875412494648E-004,  1.228096532375123E-003,  3.727799151140559E-003,...
   1.131546759143489E-002,  3.434729222834825E-002]'; % mmol C / cell
 
dname='~/Transport_Matrices/MITgcm_ECCO_v4_data/interp/PlanktonBiomass/';
lst=dir([dname 'c*.0001.nc']);

%%

% set GUD output coordinates
xx=0.25:0.5:359.75;
yy=-89.75:0.5:89.75;
zz=unique(Zboxnom);

% index to TM vector
i=dsearchn(xx',Xboxnom);
j=dsearchn(yy',Yboxnom);
k=dsearchn(zz ,Zboxnom);  
ii=sub2ind([720,360,50],i,j,k);

% extract and save data for each species
for specID=1:length(lst)
    fname=lst(specID).name;
    disp(['Processing ' fname ' (' num2str(specID) ' of ' num2str(length(lst)) ')'])
    
    data = ncread([dname fname],fname(1:3));

    % Process data month by month
    for mn=1:12
        % extract data matrix and convert from C biomass to abundance
        X=data(:,:,:,mn)./QCARBON(specID);
        % remove negatives
        X(X<0)=0;
        % reshape
        X=cat(1,X(xx>180,:,:),X(xx<=180,:,:));
        % extract to vector format
        abund=X(ii);
        
        % place in output structure
        abundance(:,mn) = double(abund);
    end
    
    sname=['../GUD_forcing/GUD_X' num2str(specID,'%02.0f') '_abundance.mat'];
    save(sname,'abundance')
end


%%
fname='~/Transport_Matrices/MITgcm_ECCO_v4_data/interp/PhysicalOceanography/THETA.0001.nc';
data = ncread(fname,'THETA');

for mn=1:12
    % extract data matrix
    X=data(:,:,:,mn);
    % reshape
    X=cat(1,X(xx>180,:,:),X(xx<=180,:,:));
    % extract to vector format
    theta(:,mn) = double(X(ii));
end

sname=['../GUD_forcing/Theta.mat'];
save(sname,'theta')
%%































