close all;
clear all;

PLOT = 1;
LOOP = false;                               % 1 - use loop times
                                        % 0 - use travel times
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
trucks = unique(truckId);

% initialize variables
alphas = 0.94; %linspace(0.001, 0.999); %[0.1, 0.5, 0.9];
betterCounter = zeros(1,length(alphas));% metric for choosing
                                        % alpha; higher is better!
worseCounter = zeros(1,length(alphas));         
betterSizeCounter = zeros(1,length(alphas));% metric for choosing
                                        % alpha; higher is better!
worseSizeCounter = zeros(1,length(alphas));         
xorCounter = zeros(1,length(alphas));
orCounter = zeros(1,length(alphas));
andCounter = zeros(1,length(alphas));

cnt2 = 1;                               % alpha counter
loopTime = NaN.*zeros(1,length(baseidx));

% Get travel times (always good)
travelTime = (towerTime'-loadTime')*24;
travelTime(find(travelTime>12)) = NaN;   % ignore overnight tickets

% Get loop times (may not be able to calc...)
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
end

for alpha = alphas                      % EWMA forgetting parameter

% clear all vars
d.truckAssignment = NaN.*zeros(length(trucks),length(baseidx)); % maxSize
d.ewma = NaN.*zeros(6,length(baseidx)); % maxSize
d.size = NaN.*zeros(6,length(baseidx)); % maxSize
d.haul = NaN.*zeros(6,length(baseidx)); % maxSize
d.qcday = NaN.*zeros(6,length(baseidx));% maxSize
d.time = NaN.*zeros(1,length(baseidx)); % maxSize
d.from = NaN.*zeros(1,length(baseidx)); % maxSize
d.to = NaN.*zeros(1,length(baseidx));   % maxSize
d.price = NaN.*zeros(1,length(baseidx));% maxSize
ewmaFilter = NaN.*zeros(1,6);
timeDiff = [];
cnt = 1;                                % event counter for d

teamSize = [];
teamHaul = [];
teamWait = [];

logitTime = [];
logitSize = [];
logitHaul = [];
logitDec = [];

% Look at EWMA estimates at each switch
for i = 1:numDays
    dayIdx = find(floor(loadTime)==day(i));
    daysTeams = unique(drt(dayIdx));
    keep = zeros(1,length(daysTeams));
    % remove non-siboney teams
    for j = 1:length(daysTeams)
        keep(j) = all(strcmp(subcont(dayIdx(find(drt(dayIdx)== ...
                                                 daysTeams(j)))),'SIB'));
    end
    daysTeams = daysTeams(find(keep));
    
    if length(daysTeams) ~= 6
        %        warning('not a 6-team day!!');
        continue;
    end
      
    % Look at today's tickets in time-order
    [~,ticketOrder] = sort(loadTime(dayIdx));

     % Clear filters daily
     ewmaFilter = NaN.*zeros(1,6);

     if PLOT
         figure;
     end
     
     start = cnt;
     
     for j = ticketOrder'
         thisTruck = find(trucks == truckId(dayIdx(j)));
         thisTeam = drt(dayIdx(j));
         teamIdx = find(daysTeams == thisTeam);
         
         if isempty(teamIdx)            % from removed team
             continue;
         end

         % Update appropriate EWMA filter...
         if LOOP
             tmp = loopTime(dayIdx(j));
         else
             tmp = travelTime(dayIdx(j));
         end
         
         if isnan(ewmaFilter(teamIdx))   % first ticket
             ewmaFilter(teamIdx) = tmp;
         elseif ~isnan(tmp)
             ewmaFilter(teamIdx) = alpha*tmp + ...
                 (1-alpha)*ewmaFilter(teamIdx);
         end

         if thisTeam ~= d.truckAssignment(thisTruck,cnt) % truck
                                                         % switched
             cnt = cnt + 1;
             
             % Copy vectors
             d.truckAssignment(:,cnt) = d.truckAssignment(:,cnt-1);
             d.qcday(:,cnt) = d.qcday(:,cnt-1);
             d.haul(:,cnt) = d.haul(:,cnt-1);

             lastteamIdx = find(daysTeams == ...
                                d.truckAssignment(thisTruck,cnt-1));
             if ~isempty(lastteamIdx)
                 d.from(cnt) = lastteamIdx;
             end
             d.to(cnt) = teamIdx;       % this will always be in
                                        % the set

             if cnt == start+1          % new day!
                 d.size(:,cnt-1) = zeros(6,1);
                 d.size(:,cnt) = zeros(6,1);
                 d.haul(:,cnt) = NaN.*zeros(6,1);
             else
                 d.size(:,cnt) = d.size(:,cnt-1);
             end
             
             % Update
             d.truckAssignment(thisTruck,cnt) = drt(dayIdx(j));
             d.time(cnt) = loadTime(dayIdx(j));
             d.qcday(teamIdx,cnt) = drt(dayIdx(j));
             d.size(teamIdx,cnt) = d.size(teamIdx,cnt) + 1;
             d.size(lastteamIdx,cnt) = d.size(lastteamIdx,cnt) - 1;
             d.ewma(:, cnt) = ewmaFilter;% grab current filter
                                         % values
             d.qcday(teamIdx,cnt) = thisTeam;
             tmp = haulMi(dayIdx(j));
             tmp2 = d.haul(lastteamIdx,cnt);
             d.haul(teamIdx,cnt) = tmp;

             if tmp >=15 & tmp2 < 15    % look at difference in
                                        % known price amounts
                 d.price(cnt) = -1;
             elseif tmp < 15 & tmp2 >= 15
                 d.price(cnt) = 1;
             else
                 d.price(cnt) = 0;
             end
         end
     end

     % loop time estimates
     tIdx = [start:cnt];
     ptime = d.time(tIdx) - day(i);

     
     changesIdx = find(~isnan(d.from(tIdx)));

     % collect inter-decision times
     timeDiff = [timeDiff, diff(sort(d.time(tIdx(changesIdx))))];

     tmp = d.ewma(:,tIdx(changesIdx));
     tmpS = d.size(:,tIdx(changesIdx));
     tmpH = d.haul(:,tIdx(changesIdx));

     % Data at switch
     fromTimes = zeros(1,length(changesIdx));
     toTimes = zeros(1,length(changesIdx));
     fromSize = zeros(1,length(changesIdx));
     toSize = zeros(1,length(changesIdx));
     fromHaul = zeros(1,length(changesIdx));
     toHaul = zeros(1,length(changesIdx));

     % Extract switch-time data
     for j = 1:length(changesIdx)
         fromTimes(j) = tmp(d.from(tIdx(changesIdx(j))),j);
         toTimes(j) = tmp(d.to(tIdx(changesIdx(j))),j);
     
         fromSize(j) = tmpS(d.from(tIdx(changesIdx(j))),j)+1; % need to correct 
         toSize(j) = tmpS(d.to(tIdx(changesIdx(j))),j)-1;

         fromHaul(j) = tmpH(d.from(tIdx(changesIdx(j))),j);
         toHaul(j) = tmpH(d.to(tIdx(changesIdx(j))),j);

         % Collect data for logistic regression...(positive is more
         % logical for all measures)
         Ltmp = zeros(6,6);
         LtmpS = zeros(6,6);
         LtmpH = zeros(6,6);

         for k = 1:6
             Ltmp(k,:) = tmp(k,j) - tmp(:,j);
             LtmpS(k,:) = tmpS(k,j) - tmpS(:,j);
             LtmpH(k,:) = tmpH(k,j) - tmpH(:,j);
         end
         
         % Hardcoded...
         noDiagIdx = [2, 3, 4, 5, 6, 7, 9,10,11,12,13,14,16,17,18, ...
                      19,20,21,23,24,25,26,27,28,30,31,32,33,34,35];

         logitTime = [logitTime, Ltmp(noDiagIdx)];
         logitSize = [logitSize, LtmpS(noDiagIdx)];
         logitHaul = [logitHaul, LtmpH(noDiagIdx)];
         decLogit = zeros(6,6);
         decLogit(d.from(tIdx(changesIdx(j))), d.to(tIdx(changesIdx(j)))) ...
             = 1;
         logitDec = [logitDec, decLogit(noDiagIdx)];
     end
     
     
     
     if PLOT
         subplot(311);
         % current 'state of knowledge'
         plot(ptime,d.ewma(:,tIdx));
         hold on;

         % old team
         plot(ptime(changesIdx), fromTimes, 'rd');
     
         % new team
         plot(ptime(changesIdx), toTimes, 'gd');
         
         axis([0,1,0,5]);
         
         subplot(312), 
         y = [];
         for j=1:length(daysTeams)
             [x,tmp] = stairs(ptime, d.size(j,tIdx));
             y = [y,tmp];
         end
         plot(x,y');
         axis([0,1,0,20]); 
         
         subplot(313),
         plot(ptime,d.haul(:,tIdx));
         axis([0,1,0,30]); 
     end

     % Accumulate 'successful' switches
     betterCounter(cnt2) = betterCounter(cnt2) + sum(fromTimes>=toTimes);
     worseCounter(cnt2) = worseCounter(cnt2) + sum(fromTimes<toTimes);
     betterSizeCounter(cnt2) = betterSizeCounter(cnt2) + sum(fromSize>toSize);
     worseSizeCounter(cnt2) = worseSizeCounter(cnt2) + sum(fromSize<=toSize);
     xorCounter(cnt2) = xorCounter(cnt2) + sum(xor(fromSize>toSize, ...
                                                   fromTimes>=toTimes));
     orCounter(cnt2) = xorCounter(cnt2) + sum(or(fromSize>toSize, ...
                                                 fromTimes>=toTimes));
     andCounter(cnt2) = andCounter(cnt2) + sum(and(fromSize>toSize, ...
                                                   fromTimes>=toTimes));
     
     %max(d.size(:,tIdx)')
     %ceil(median(d.size(:,tIdx)'))
     %ceil(mean(d.size(:,tIdx)'))
     
     teamSize = [teamSize, ceil(nanmean(d.size(:,tIdx)'))];
     teamHaul = [teamHaul, nanmean(d.haul(:,tIdx)')];
     teamWait = [teamWait, nanmean(d.ewma(:,tIdx)')];
end
cnt2 = cnt2 + 1;
end

figure; 
plot(alphas, betterCounter./(betterCounter + worseCounter), 'b.');
title('Proportion of Rational Decisions by EWMA Factor');
xlabel('Forgetting Factor');
ylabel('Proportion');

placefigures;

% Trim data
d.truckAssignment = d.truckAssignment(:,[1:cnt]);
d.ewma = d.ewma(:,[1:cnt]);
d.size = d.size(:,[1:cnt]);
d.haul = d.haul(:,[1:cnt]);
d.qcday = d.qcday(:,[1:cnt]);
d.time = d.time(:,[1:cnt]);
d.from = d.from(:,[1:cnt]);
d.to = d.to(:,[1:cnt]);
d.price = d.price(:,[1:cnt]);

% purge > 1 day diffs
timeDiff(timeDiff>1) = NaN;
figure, hist(timeDiff*24);
title('Inter-Decision Times')
xlabel('Hours');

[obs,y] = hist(timeDiff*24)
edges=cumsum(repmat(y(2)-y(1),10,1));
expected = sum(y)*(1-expcdf(edges,0.44));
chisq = sum((obs'-expected).^2./expected);

figure, stem([-20:20], xcorr(timeDiff(find(~isnan(timeDiff))),20,'coeff'))
title('Autocorrelation of time differences');

figure, scatter(teamHaul, teamWait, teamSize*100, 'b.')