% Load data, extract deicison model parameters
%
% 15 January 2013
% J.Brooks

clear all;
close all;

km2mi = 0.621371;

%load DEBRIS_ANALSIS_AL_COMPLETE
%load DEBRIS_ANALSIS_AL_PEAK_19feb %for verification
load DEBRIS_ANALYSIS_AL % for verification

% Remove unneeded data to save memory
clear CLEAN PLOT QCTmp ROE ROETmp SCRUB_ROE SCRUB_WEIGHT
clear SCRUB_WHITES SCRUB_ZERO_HAUL SYS_PLOTS TDSR_PLOTS
clear TRUCK_PLOTS ahaul ans arrivalRate attd avgLoadPercent bin
clear busyTime change color contentPercent contests 
clear corrected cyL cyT d data dataType day1 dayQCs 
clear directory edged effLoads endTime file filename good i idx
clear idxAM idxC idxL idxL interTimes interToday invoice j k l
clear lines linesI loadIdx loadPercentTmp loadsPerDay maxCapacity
clear meanService n noLine nonZero notcombo notsorted numQCs
clear numTdsrs numTicketsL numTicketsT numTrucksL numTrucksT
clear r ratio sitePercent startTime sub subInvoice subcont2
clear t tickets tmp tmp2 tmpI today total totalTickets totalTime 
clear totals truckIdTmp txt u utilization

clear equality % needed only for peak data (used for verification)

PLOT = 0;

% Overwrite to focus on main event (many zero days, very low activity..)
%duration = [40672:40817]; % Complete data set
duration = [40695:40720]; % original June data (exclude 40721 which
                          % has 10 trucks only?)

lonMin = -88.4;
lonMax = -85.4;
latMin = 32.22;
latMax = 35;

dayCount = length(duration);

% team identifier
QCday = QC + 1e6*floor(loadTime);
teams = unique(QCday);
teamDelay = NaN*zeros(1, length(teams));
teamHaul = NaN*zeros(1,length(teams));
teamLat = NaN*zeros(1,length(teams));
teamLon = NaN*zeros(1,length(teams));

disp(sprintf('Num unique QCdays: %d', length(teams)));

% ----------------------------------------
% Calculate mean team round-trip delays
% ----------------------------------------
for i = 1:length(duration)
    todayData = find(floor(loadTime) == duration(i) & ...
                     floor(towerTime) == duration(i) & ...
                     lat~=0  & lon~=0 & lat>latMin & lat<latMax ...
                     & lon>lonMin & lon<lonMax);

    todaysTeams = unique(QCday(todayData));
    %    todaysTrucks = unique(trucksId(todayData));

    for j = todaysTeams'
        teamData = find(QCday(todayData) == j);
        idx = todayData(teamData);
        teamtrucks = unique(truckId(idx));
        teamMean = NaN.*zeros(1,length(teamtrucks));
        truckCnt = 1;
        
        for k = teamtrucks'
            % Ticket data for truck k on day i for team j
            idx2 = idx(find(truckId(idx) == k));
            
            if length(idx2) > 1
                teamMean(truckCnt) = median(diff(sort(loadTime(idx2))))*24; 
                % NOTE: in hrs
                % switched to median 1/23/2013
            end
            truckCnt = truckCnt + 1;
        end

        % FOR DEBUGGING...
        %        figure, plot(teamMean, 'b.');
        
        teamIdx = find(teams == j);
        %        title(teamIdx);

        teamDelay(teamIdx) = nanmedian(teamMean); % switched from
                                                  % mean to median 1/23/2013
        teamHaul(teamIdx) = mean(haulMi(idx));
        teamSize(teamIdx) = length(teamtrucks);
        teamLat(teamIdx) = mean(lat(idx));
        teamLon(teamIdx) = mean(lon(idx));
        
        %figure, plot(lat(teamData),lon(teamData), 'b.');
        %hold on;
        %plot(teamLat(teamIdx), teamLon(teamIdx), 'r*');
        %title(teamIdx);
        
        if teamLat(teamIdx) == 0 || teamLon(teamIdx) == 0
            lat(teamData)
            lon(teamData)
            QC(teamData)
            floor(loadTime(teamData))
        end
    end
    
    if PLOT
        figure(1), scatter(i + 0.5*rand(1,length(todaysTeams)), ...
                           teamDelay(match(teams, todaysTeams)), 20, ...
                           teamHaul(match(teams, todaysTeams)), 'filled');
        hold on;
        figure(2), plot(i, length(todaysTeams), 'b.');
        hold on;
    end
end


% ----------------------------------------
% Look at switching behavior
% ----------------------------------------
numSwitchCounts = zeros(1,11);        % [0, ..., 10] switches (truck-days)

if PLOT
    figure(3);
    plot([0 7], [0 7], 'r');
    title('Switching Behavior -- Round-Trip');
    xlabel('First Team'), ylabel('Second Team');
    hold on;

    figure(4);
    plot([0 30], [0 30], 'r');
    title('Switching Behavior -- Haul');
    xlabel('First Team'), ylabel('Second Team');
    hold on;

    figure(5);
    plot([0 20], [0 20], 'r');
    title('Switching Behavior -- Size');
    xlabel('First Team'), ylabel('Second Team');
    hold on;
end

% One-way
sizeSmaller = 0;
roundTripBetter = 0;
haulBetter = 0;
anyBetter = 0;

% all interactions
better = zeros(1,8);
usubs = unique(subcont);
numSubs = length(usubs);
betterCont = zeros(numSubs,8);          % contingency table to
                                        % check for differences
                                        % between subcontractor
                                        % strategies 
decisions = 0;

% distances
distanceMoved = [];

for i = 1:length(duration)
    todayData = find(floor(loadTime) == duration(i) & ...
                     floor(towerTime) == duration(i) & ...
                     lat~=0  & lon~=0 & lat>latMin & lat<latMax ...
                     & lon>lonMin & lon<lonMax);

    %    todaysTeams = unique(QCday(todayData));
    todaysTrucks = unique(truckId(todayData));
    
    for j = todaysTrucks'
        %teamsWorked = unique(QC(todayData(find(truckId(todayData) = ...
        %                                      = j))));
        thisData = todayData(find(truckId(todayData) == j));
        [time,order] = sort(loadTime(thisData));
        qcs = QC(thisData);
        qcs = qcs(order);               % make sure properly sorted
        
        switches = find(diff(qcs) ~= 0);
        
        numSwitchCounts(length(switches)+1) = ...
            numSwitchCounts(length(switches)+1) + 1;
        
        % Uncomment to look at only single-switch cases
        % if length(switches)>1 continue; end

        for k = switches'               % look at all switches
            anyTmp = 0;
            timeTmp = 0;
            haulTmp = 0;
            sizeTmp = 0;
            
            qctmp1 = qcs(k) + 1e6*duration(i);
            qctmp2 = qcs(k+1) + 1e6*duration(i);

            first = find(teams == qctmp1);
            second = find(teams == qctmp2);
            
            if isnan(teamHaul(first)) || isnan(teamHaul(second))
                continue;
            end
            
            decisions = decisions + 1;
            
            decisionsStruc.first(decisions) = first;
            decisionsStruc.second(decisions) = second;
            decisionsStruc.truck(decisions) = j;
            
            if PLOT
                figure(3), plot(teamDelay(first), ...
                                teamDelay(second), 'b.');
                figure(4), plot(teamHaul(first), ...
                                teamHaul(second), 'b.');
                figure(5), plot(teamSize(first)+0.6*rand-0.3, ...
                                teamSize(second)+0.6*rand-0.3, ...
                                'b.');
            end

            distanceMoved = [distanceMoved, posdist(teamLat(first), ...
                                                    teamLon(first), ...
                                                    teamLat(second), ...
                                                    teamLon(second))*km2mi];

            if teamDelay(second) < teamDelay(first)
                roundTripBetter = roundTripBetter + 1;
                timeTmp = 1;
                anyTmp = 1;
            end
            
            if teamHaul(second) < teamHaul(first)
                haulBetter = haulBetter + 1;
                haulTmp = 1;
                anyTmp = 1;
            end
            
            if teamSize(second) < teamSize(first)
                sizeSmaller = sizeSmaller + 1;
                sizeTmp = 1;
                anyTmp = 1;
            end

            anyBetter = anyBetter + anyTmp;
            interactionTmp = haulTmp + 2*timeTmp + 4 *sizeTmp + 1;
            better(interactionTmp) = better(interactionTmp) + 1;
            
            subI = find(usubs == subcont(thisData(1)));% subcontractor index

            betterCont(subI, interactionTmp) = ... 
                betterCont(subI, interactionTmp) + 1;
            

        end
    end
end

if PLOT
    figure, bar([0:10], numSwitchCounts);
    title('Histogram of number of switches per truck-day');
    xlabel('Number of Switches');
    ylabel('Count');
    
    figure, hist(distanceMoved);
    title('Histogram of distance moved');
end

disp(sprintf('Number of decisions analyzed: %d', decisions));
disp(sprintf('Haul Smaller: %2.2f', haulBetter/decisions*100));
disp(sprintf('Round-Trip Time Smaller: %2.2f', roundTripBetter/ ...
             decisions*100));
disp(sprintf('Team Smaller: %2.2f', sizeSmaller/decisions*100));
disp(sprintf('Any of above: %2.2f', anyBetter/decisions*100));

% Analysis of contingency table (looking for differences in
% contractor strategy)

%idx = find(sum(betterCont,2)>50);      % average >5
 
 idx = find(all(betterCont' > 2));      % strict interpretation
 
 percCont = zeros(length(idx),8);
 cnt = 1;
 for i = idx
     percCont(cnt,:) = betterCont(i,:)./sum(betterCont(i,:));
     cnt = cnt + 1;
 end
 
 figure, bar(percCont', 1)
 
 