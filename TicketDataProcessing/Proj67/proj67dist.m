close all;
clear all;

load DEBRIS_AL_COMPLETE

baseidx = find(project == 6 | project == 7);

duration = [40676:40726];
durationIdx= [4, 8, 10, 14, 19, 24, 31, 33, 34, 44]; % from
                                                     % systemModelValidationData
                                                     % script
day = duration(durationIdx);
numDays = length(unique(day));
% extract data
baseidx = baseidx(find(~isnan(match(day, floor(loadTime(baseidx))))));
lat = lat(baseidx);
lon = lon(baseidx);
loadTime = loadTime(baseidx);
towerTime = towerTime(baseidx);
QC = QC(baseidx);
truckId = truckId(baseidx);
haulMi = haulMi(baseidx);
loadVolume = loadVolume(baseidx);
subcont = subcont(baseidx);
tdsr = tdsr(baseidx);

tdsrs = unique(tdsr);
numTdsrs = length(tdsrs);

%eqData = equality(lat,lon,duration,loadTime,QC,truckId,haulMi, ...
%                   capacity,loadVolume,trucks);
 
drt = QC + 1e6*floor(loadTime);
drts = unique(drt);%eqData.QCDay; 
subs = subcont;

for i = 1:length(day)
    disp(sprintf('Day: %d', day(i)));
    
    teams = unique(QC(find(floor(loadTime)==day(i))));
    disp('Mean Distances');
    disp(sprintf('Team\t 601\t 602'));
    for j = 1:length(teams)
        d601 = mean(haulMi(find(floor(loadTime)==day(i) & tdsr==601 ...
                                & QC == teams(j))));
        d602 = mean(haulMi(find(floor(loadTime)==day(i) & tdsr==602 ...
                                & QC == teams(j))));

        disp(sprintf('%d\t %2.2f\t %2.2f',j, d601, d602));
    end        
end
  


