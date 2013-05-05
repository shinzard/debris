close all;
clear all;

%load DEBRIS_ANALSIS_AL_COMPLETE % has dups...
load DEBRIS_AL_COMPLETE
dayDataORIG = dayData;                  % this is overwritten below...

idx = find(project == 24);
baseidx = idx;

duration = [40683:40784];

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

eqData = equality(lat,lon,duration,loadTime,QC,truckId,haulMi, ...
                   capacity,loadVolume,trucks);

 
drt = QC + 1e6*floor(loadTime);
drts = eqData.QCDay; %unique(drt);
subs = subcont;

% initialize variables
drtHaul = NaN.*zeros(1,length(drts));
drtDelay = NaN.*zeros(1,length(drts));
drtHaulRel = NaN.*zeros(1,length(drts));
drtDelayRel = NaN.*zeros(1,length(drts));
drtPercAlloc = NaN.*zeros(1,length(drts));
drtDay = NaN.*zeros(1,length(drts));

numSubs = zeros(1,length(duration));
numTrucks= zeros(1,length(duration));
numDrts = zeros(1,length(duration));
numDrtsSib = zeros(1,length(duration));
numDrtsOther = zeros(1,length(duration));

% -----------------------------------------------------------------
% Characterize percent allocation (not looking at switches directly
% but rather a process analogous to optimal method...
% -----------------------------------------------------------------
teamIdx = 1;
% Loop through all days of mission
for i = 1:length(duration)
    dayData = find(floor(loadTime) == duration(i));
    dayDrts = unique(drt(dayData));
    dayDrts = dayDrts(find(~isnan(match(drts,dayDrts))));% make sure consistent with
                                                         % equality DRT list...

    numSubs(i) = length(unique(subs(dayData)));
    numDrts(i) = length(dayDrts);
    numTrucks(i) = length(unique(truckId(dayData)));
    
    sibCount = 0;
    otherCount = 0;
    totalTrucks = 0;
    
    % loop through all teams on this day
    for j = 1:numDrts(i)
        tmpIdx = find(drt(dayData) == dayDrts(j)); 
        tmpSubs = unique(subs(dayData(tmpIdx)));
        SIB_TEAM = false;
        
        if length(tmpSubs) == 1
          if strcmp(tmpSubs, 'HRB')     % HRB is the main sub (75%)
              sibCount = sibCount + 1;  % only analyze these
                                        % teams??
              SIB_TEAM = true;
          else
              otherCount = otherCount + 1;
          end
        end
        
        teamtrucks = unique(truckId(dayData(tmpIdx))); 

        nTrucks = 0;
        % Find initial allocation:
        for k = 1:length(teamtrucks);
            % find all loads for truck k on day i
           idx = find(truckId(dayData) == ...
                      teamtrucks(k));
           [t, idx2] = sort(loadTime(dayData(idx)));
           if QC(dayData(idx2(1))) == QC(dayData(tmpIdx(1)))
               nTrucks = nTrucks + 1;
           end
        end
            
        totalTrucks = totalTrucks + length(teamtrucks);
        
        if SIB_TEAM
            teamMean = NaN.*zeros(1,length(teamtrucks));
            truckCnt = 1;

            % Loop through all truck data to get mean delay
            for k = teamtrucks'
                % Ticket data for truck k on day i for team j
                idx2 = find(truckId(dayData(tmpIdx)) == k);
            
                if length(idx2) > 1
                    teamMean(truckCnt) = ...
                        median(diff(sort(loadTime(dayData(tmpIdx(idx2))))))*24; 
                    % NOTE: in hrs
                    % switched to median 1/23/2013
                    if teamMean(truckCnt) == 0
                        tmpIdx
                        idx2
                        dayData(tmpIdx(idx2))
                        loadTime(dayData(tmpIdx(idx2))) ...
                            - duration(i)
                        truckId(dayData(tmpIdx(idx2)))
                        haulMi(dayData(tmpIdx(idx2)))
                        project(dayData(tmpIdx(idx2)))
                    end
                end
                truckCnt = truckCnt + 1;
            end
            drtHaul(teamIdx) = nanmean(haulMi(dayData(tmpIdx)));
            drtPercAlloc(teamIdx) = nTrucks/numTrucks(i);
            drtDelay(teamIdx) = nanmean(teamMean);
            drtDay(teamIdx) = duration(i);
        end
        teamIdx = teamIdx + 1;
    end
    
    % Calculate relative team characteristics
    todaysDrts = find(drtDay == duration(i));
    meanHaul = mean(drtHaul(todaysDrts));
    meanDelay = mean(drtDelay(todaysDrts));
    for j=1:length(todaysDrts)
        drtHaulRel(todaysDrts(j)) = drtHaul(todaysDrts(j))./meanHaul;
        drtDelayRel(todaysDrts(j)) = drtDelay(todaysDrts(j))./meanDelay;
    end
    
    numDrtsSib(i) = sibCount;
    numDrtsOther(i) = otherCount;
    
    %    if totalTrucks ~= numTrucks(i)
    %        totalTrucks
    %        numTrucks(i)
    %        warning('somethings wrong...');
    %    end
end

figure, plot([1:length(duration)], numSubs, 'b.');
hold on;
plot([1:length(duration)], numTrucks, 'r.');
figure, plot([1:length(duration)], numDrts, 'g.');
hold on;
plot([1:length(duration)], numDrtsSib, 'bs');
plot([1:length(duration)], numDrtsOther, 'ks');

figure, subplot(211);
plot(drtHaul, drtPercAlloc, 'b.');
title('Allocation by Haul Distance');
subplot(212), plot(drtDelay, drtPercAlloc, 'b.');
title('Allocation by Round-trip Delay');

% Look only at 6-team (siboney only) days
dS = decisions(QC, towerTime, loadTime, lat, lon,truckId,haulMi, subcont, duration(find(numDrtsSib==6)), ...
               0, 100, eqData.QCDay, tdsr);

%dS.include = ~isnan(match(drts,dS.teamID));% make sure consistent with
                                           % equality DRT list...


% indices (out of 314 teams) of 60 teams of interest (10 days of
% 6-team data in isolated region)
teamIdx = find((~isnan(dS.teamSize) & (~isnan(drtDay))));

% ASSUMPTION: trucks only active in restricted region for
% duration....dayData contains all data....
fData = fluidity3(dayDataORIG,trucks,loadTime,truckId,subcont, ...
                  loadVolume,eqData,QC,capacity, dayService, ...
                  0,[],duration,0, towerTime,lat,lon,loadPercent, haulMi,0);

% available performance data (only 51/60...)
perfDataIdx = find(~isnan(match(teamIdx,fData.teamIdx)));
%teamIdx = fData.teamIdx(perfDataIdx); % uncomment to have same
% dimensions...but will likely want all of these for the model (all
% numeric travelTime...)

% Extract all relevant measures 
throughput = 10.^(fData.effectiveness(perfDataIdx));
size = fData.sizeF(perfDataIdx);
day = duration(fData.dayF(perfDataIdx));
fluidity = fData.fluidF(perfDataIdx);
cap = fData.capacity(perfDataIdx);
travelTime = dS.teamTravel(:,teamIdx);
loopTime = dS.teamDelay(:,teamIdx);
probs = dS.teamTdsrProb(:,teamIdx);

disp('Validation Data by Team');
disp(sprintf(['Day\t Team\t Size\t AvgCap\t Fluidy\t T.put\t 601\t 602\t T(601)\t ' ...
              'T(602)']));

s=zeros(1,length(unique(day))*6);
f=zeros(1,length(unique(day))*6);

cnt = 1;

for i = unique(day)
    idx = find(day==i);
    numTeams = length(idx);
    for j = 1:numTeams
        disp(sprintf('%d\t %d\t %d\t %2.2f\t %2.2f\t %2.1f\t %d\t %d\t %2.2f\t %2.2f',i,j, ...
                     size(idx(j)),cap(idx(j))/size(idx(j))/45, ...
                     fluidity(idx(j)),throughput(idx(j)),probs(1,idx(j)), ...
                     probs(2,idx(j)), travelTime(1,idx(j)), ...
                     travelTime(2,idx(j))));
        s(cnt) = size(idx(j));
        f(cnt) = fluidity(idx(j));
        cnt = cnt + 1;
    end
    
    for j = (numTeams+1):6
        disp(sprintf('%d\t %d',i,j));   % blank lines
    end
end

disp('Summary of Validation Data');
disp(sprintf('Day\t Size\t\t AvgCap\t\t Fluidy\t\t Trips'));

for i = unique(day)
    idx = find(day==i);
    disp(sprintf('%d\t %2.2f(%d)\t %2.2f(%2.2f)\t %2.2f(%2.2f)\t %2.2f(%d)\t',i, ...
                 nanmean(size(idx)),range(size(idx)),nanmean(cap(idx)./size(idx)),range(cap(idx)./size(idx)), ...
                 nanmean(fluidity(idx)),range(fluidity(idx)), ...
                 nanmean(probs(1,idx)+probs(2,idx)),range(probs(1,idx)+probs(2,idx))));
end


% Look at travel/wait time estimates....
for i = 1:length(unique(tdsr))
    for j = 1:length(unique(day))
        tmp = [1:6]+6*(j-1);
        figure(120), plot(j,nanmean(travelTime(i,tmp)), 'b.'), hold on, 
        plot(j,nanmean(loopTime(i,tmp)), 'rs'), plot(j,nanmean(travelTime(i,tmp)*2), 'g.');
        legend({'One-way Travel Time', 'Loop Time', ['One-way ' ...
                            'Travel*2']});
        
        figure(121), plot(j,nanmean(loopTime(i,tmp) - travelTime(i,tmp)*2), ...
                          'b.');
        hold on;
        title(['Difference between travel*2 and loop -- Estimated wait time ' ...
               'at TDSR']);

        lambda = sum(probs(:,tmp),2)/12;
        W = nanmean(loopTime(i,tmp) - travelTime(i,tmp)*2,2);
        disp(sprintf('Mu estimates: %2.2f, %2.2f',1./W + lambda))
    end
end
