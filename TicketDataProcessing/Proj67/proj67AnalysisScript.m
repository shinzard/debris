close all;
clear all;

%load DEBRIS_ANALSIS_AL_COMPLETE % has dups...
load DEBRIS_AL_COMPLETE

%idx = find((project == 6 | project == 7) & (file ~= 50) );
idx = find(project == 6 | project == 7);
baseidx = idx;

duration = [40676:40726];

lat = lat(baseidx);
lon = lon(baseidx);
loadTime = loadTime(baseidx);
towerTime = towerTime(baseidx);
QC = QC(baseidx);
truckId = truckId(baseidx);
haulMi = haulMi(baseidx);
loadVolume = loadVolume(baseidx);

eqData = equality(lat,lon,duration,loadTime,QC,truckId,haulMi, ...
                   capacity,loadVolume,trucks);


% ASSUMPTION: trucks only active in restricted region for
% duration....dayData contains all data....
[effectiveness, efficiency, equality,sizeF,familiarityF,haul,teamHist,teamDay] = ...
fluidity3(dayData,trucks,loadTime,truckId,subcont, ...
          loadVolume,eqData,QC,capacity, dayService, ...
          0,['Proj67Data_',num2str(now), '.csv'],duration,0, ...
          towerTime,lat,lon,loadPercent, haulMi);

% Remove confusion....13 March 2013
fluidF = 1 - familiarityF;