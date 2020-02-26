clear
clc


Biodiversity_type = 'Metagenomic'; % 'Metagenomic or OTU
Size_class = 3%1:6; % 1 to 6

fnames = dir('~/Data/Literature/Richter_2020/*dissimilarity*');

if startsWith(Biodiversity_type,'Metagenomic')
    fnames = dir('~/Data/Literature/Richter_2020/*metagenomic_dissimilarity*');
elseif startsWith(Biodiversity_type,'OTU')
    fnames = dir('~/Data/Literature/Richter_2020/*OTU_dissimilarity*');
else
    error('Biodiversity_type must be ''Metagenomic'' or ''OTU''.')
end

Community_filename  = fnames(Size_class).name;
Env_filename        = 'supplementary_table_02.environmental_parameters.mat';
Geographic_filename = 'supplementary_table_16.geographic_distance.mat';

pathname='~/Data/Literature/Richter_2020/';
addpath(pathname)
%% Environmental data
load(Env_filename);
Env_data      = Data;
Env_stations  = char(Env_data.Station);
Env_stationID = str2num(Env_stations(:,6:end));
Env_layer     = Env_data.Depth;
Env_Temp      = normalize(Env_data.TemperatureC);
Env_Nitr      = normalize(Env_data.NO2NO3molL);
Env_Phos      = normalize(Env_data.PhosphatemolL);
Env_Iron      = normalize(Env_data.Ironmmolm3);
Env_Sili      = normalize(Env_data.SilicateWOA13molL);
Env_Chlo      = normalize(Env_data.Chlorophyllmgm3);
OceanBasin    = Env_data.Oceanandsearegion;
Biome         = Env_data.Marinebiome;
Latitude      = Env_data.Latitude;
Longitude     = Env_data.Longitude;

EnvV=[Env_Temp Env_Nitr Env_Iron Env_Sili Env_Chlo];
EnvV(isnan(EnvV))=0;

% isSO=find(startsWith(cellstr(OceanBasin),'Southern Ocean'));
% EnvV(isSO,:)=0;

%% Community data
load(Community_filename);
Community_data      = Data;
Community_distance  = Community_data{:,2:end};
Community_stations  = char(Community_data.VarName1);
Community_stationID = str2num(Community_stations(:,1:3));
Community_layer     = cellstr(Community_stations(:,5:end));

%% Geographic data
load(Geographic_filename);
Geographic_data      = Data;
Geographic_distance  = Geographic_data{:,2:end};
Geographic_distance(isnan(Geographic_distance))=0;
Geographic_stations  = char(Geographic_data.VarName1);
Geographic_stationID = str2num(Geographic_stations(:,1:3));
insurface            = find(startsWith(Community_layer,'SRF'));

%% Collate distance matrices
% find intersect of environmental and community data
[C,ia,ib] = intersect(Env_stationID,Community_stationID(insurface));

OBind = grp2idx(OceanBasin(ia));
EnvVars=EnvV(ia,:);

Edist = squareform(pdist(EnvVars));
Mdist = Community_distance(ib,ib);
Gdist = Geographic_distance(C,C);

species=(OceanBasin(ia));

opts = statset(...
    'OutputFcn',@(optimValues,state) KLLogging(optimValues,state,species),...
    'TolFun',1e-10,...
    'MaxIter',10000);

Y = tsne(Mdist,'Algorithm','exact',...
    'NumDimensions',2,...
    'Distance',@dist2dist,...
    'Perplexity',20,...
    'Options',opts)
