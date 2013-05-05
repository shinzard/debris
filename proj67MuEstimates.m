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

% initialize variables
loopTime = NaN.*zeros(1,length(baseidx));
travelTime = (towerTime'-loadTime')*24;
lambda = NaN.*zeros(numDays,numTdsrs);
mu = NaN.*zeros(numDays,numTdsrs);
wEst = NaN.*zeros(numDays,numTdsrs);
W = NaN.*zeros(length(baseidx),numTdsrs);
cnt = 1;

travelTime(find(travelTime>12)) = NaN;   % ignore overnight tickets

% Look at travel/wait time estimates....
for i = 1:numDays
    dayIdx = find(floor(towerTime)==day(i));
    dayTrucks = truckId(dayIdx);
    todaysTrucks = unique(dayTrucks);
    
    for j = 1:length(todaysTrucks)
        thisTruckData = dayIdx(find(truckId(dayIdx) == ...
                                    todaysTrucks(j)));
        if length(thisTruckData) > 1    % else, cannot estimate
                                        % anything...
            [sortedTData,ticketOrderIdx] = sort(loadTime(thisTruckData));

            for k = 1:(length(thisTruckData)-1)
                thisIdx = ticketOrderIdx(k);
                nextIdx = ticketOrderIdx(k+1);

                if drt(thisTruckData(nextIdx)) == ...
                        drt(thisTruckData(thisIdx)) % if returned to same team...
                    loopTime(thisTruckData(thisIdx)) = ...
                        (loadTime(thisTruckData(nextIdx))- ...
                        loadTime(thisTruckData(thisIdx)))*24;
                    
                    if loopTime(thisTruckData(thisIdx)) < 0.05 % trailer...
                        loopTime(thisTruckData(thisIdx)) == NaN; % will be picked up on the next iteration
                    end
                end
            end
        end
    end
    
    for k = 1:numTdsrs
        lambda(i,k) = length(find(tdsr(dayIdx) == tdsrs(k)))/12; % open for 12hrs...
    end
     
end

colors = ['br'];
figure;
for k = 1:numTdsrs
    idx = find(tdsr == tdsrs(k));
    W(idx,k) = loopTime(idx) - travelTime(idx)*2;
    plot(idx, loopTime(idx) - travelTime(idx)*2, [colors(k),'.']) 
    hold on;
end

W(find(W<0)) = NaN;                     % other ideas what to do here??

for i = 1:numDays
    dayIdx = find(floor(towerTime)==day(i));
    wEst(i,:) = nanmean(W(dayIdx,:));
    mu(i,:)=1./wEst(i,:) + lambda(i,:);
end


disp('Mu Estimates')
muEst = nanmean(mu);
disp(sprintf('TDSR\t mu'));
for k = 1:numTdsrs
    disp(sprintf('%d\t %2.2f', tdsrs(k), muEst(k)));
end


