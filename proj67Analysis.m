close all;
clear all;

%load DEBRIS_ANALSIS_AL_COMPLETE % has dups...
load DEBRIS_AL_COMPLETE

%idx = find((project == 6 | project == 7) & (file ~= 50) );
idx = find(project == 6 | project == 7);
baseidx = idx;

duration = [40676:40726];
%[min(floor(towerTime(idx))):max(floor(towerTime(idx)))];

 
tdsr602 = find(tdsr(idx)==602);
tdsr601 = find(tdsr(idx)==601);
t602 = sort(towerTime(idx(tdsr602)));
t601 = sort(towerTime(idx(tdsr601)));

% Remove truck/trailor (18 second cutoff)
idx = find([1, diff(t601)'*24]>0.005); 
t601 = t601(idx);

if 0

plotArrivals(t601, floor(t601))
%plotArrivals(t602, floor(t602))
%%
days = unique(floor(t601)');
for i = days
    idx = find(floor(t601)==i);
    if length(idx)>150
        figure, hist(diff(t601(idx)*24),100);
        title(sprintf('Day: %d',i));
    end
end
%% 
end
drt = QC(baseidx) + 1e6*floor(loadTime(baseidx));
%QC(baseidx)*1e6+floor(towerTime(baseidx));
drts = unique(drt);
subs = subcont(baseidx);


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
    dayData = find(floor(loadTime(baseidx)) == duration(i));
    dayDrts = unique(drt(dayData));

    numSubs(i) = length(unique(subs(dayData)));
    numDrts(i) = length(dayDrts);
    numTrucks(i) = length(unique(truckId(baseidx(dayData))));
    
    sibCount = 0;
    otherCount = 0;
    totalTrucks = 0;
    
    % loop through all teams on this day
    for j = 1:numDrts(i)
        tmpIdx = find(drt(dayData) == dayDrts(j)); 
        tmpSubs = unique(subs(dayData(tmpIdx)));
        SIB_TEAM = false;
        
        if length(tmpSubs) == 1
          if strcmp(tmpSubs, 'SIB')     % Siboney is the main sub
              sibCount = sibCount + 1;  % only analyze these
                                        % teams??
              SIB_TEAM = true;
          else
              otherCount = otherCount + 1;
          end
        end
        
        teamtrucks = unique(truckId(baseidx(dayData(tmpIdx)))); 

        nTrucks = 0;
        % Fina initial allocation:
        for k = 1:length(teamtrucks);
            % find all loads for truck k on day i
           idx = find(truckId(baseidx(dayData)) == ...
                      teamtrucks(k));
           [t, idx2] = sort(loadTime(baseidx(dayData(idx))));
           if QC(baseidx(dayData(idx2(1)))) == QC(baseidx(dayData(tmpIdx(1))))
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
                idx2 = find(truckId(baseidx(dayData(tmpIdx))) == k);
            
                if length(idx2) > 1
                    teamMean(truckCnt) = ...
                        median(diff(sort(loadTime(baseidx(dayData(tmpIdx(idx2)))))))*24; 
                    % NOTE: in hrs
                    % switched to median 1/23/2013
                    if teamMean(truckCnt) == 0
                        tmpIdx
                        idx2
                        baseidx(dayData(tmpIdx(idx2)))
                        loadTime(baseidx(dayData(tmpIdx(idx2)))) ...
                            - duration(i)
                        truckId(baseidx(dayData(tmpIdx(idx2))))
                        haulMi(baseidx(dayData(tmpIdx(idx2))))
                        project(baseidx(dayData(tmpIdx(idx2))))
                    end
                end
                truckCnt = truckCnt + 1;
            end
            drtHaul(teamIdx) = nanmean(haulMi(baseidx(dayData(tmpIdx))));
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

decisionsStruc = decisions(QC(baseidx), towerTime(baseidx), ...
                           loadTime(baseidx), lat(baseidx), ...
                           lon(baseidx),truckId(baseidx), ...
                           haulMi(baseidx), subcont(baseidx), duration, ...
                           0, 100);

% -----------------------------------------------------------------
% Look at net flows between teams
% -----------------------------------------------------------------
numDec = length(decisionsStruc.first)
flows.in = zeros(1,length(drts));
flows.out = zeros(1,length(drts));
for i = 1:numDec
   flows.out(decisionsStruc.first(i)) = ...
       flows.out(decisionsStruc.first(i)) + 1;
   flows.in(decisionsStruc.second(i)) = ...
       flows.in(decisionsStruc.second(i)) + 1;
end

flows.net = flows.in-flows.out
figure, bar(flows.net);
% Add alternate shading to show days
last = 0;
for i = 1:length(duration)
    idx = find(drtDay == duration(i));
    if ~last
        beg = idx(1)-0.5;
        en = idx(end)+0.5;
        patch('XData', [beg en en beg], 'YData', [-10 -10 10 10], ...
              'FaceColor', 'k', 'FaceAlpha', 0.05);
        last = 1;
    else
        last = 0;
    end
end

title('Net Vehicle Flow for each DRT');
ylabel('Net Vehicle Flows'), xlabel('DRT Id');

figure, subplot(211);
plot(drtHaulRel,flows.net, 'b.');
title('Relative Haul vs. Net Flow');
subplot(212);
plot(drtDelayRel,flows.net, 'b.');
title('Relative Delay vs. Net Flow');

% Look for significant trends
good = find(~isnan(drtDelayRel) & ~isnan(flows.net));
[model, p] = fit(drtDelayRel(good)', flows.net(good)', 'poly1')

good = find(~isnan(drtHaulRel) & ~isnan(flows.net));
[model, p] = fit(drtHaulRel(good)', flows.net(good)', 'poly1')

% look at fit within each day...
%for i = 1:length(duration)
%    idx = find(drtDay == duration(i));
%    good = find(~isnan(drtHaulRel(idx)) & ~isnan(flows.net(idx)));
%    [model, p] = fit(drtHaulRel(idx(good))', flows.net(idx(good))', 'poly1')
%    figure(200), subplot(211);
%    plot(i,model.p1, 'b.'), hold on;
%    subplot(212);
%    plot(i,model.p2, 'b.'), hold on;
%end


% 1 April 2013 work
drtNum = NaN*zeros(1,length(drts));

j = 1;
for i = 1:length(numDrts)
    drtNum(j:j+numDrts(i)-1)=numDrts(i);
    j = j + numDrts(i);
end

1./numDrts                              % uniform allocation percentages