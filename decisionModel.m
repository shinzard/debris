% Load data, extract deicison model parameters
%
% 15 January 2013
% J.Brooks

clear all;

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

% pare down struc (many zero entries)

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
teamDelay = zeros(1, length(teams));
teamHaul = zeros(1,length(teams));

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
            % Ticket data for truck k on day i
            idx2 = todayData(find(truckId(todayData) == k));

            if length(unique(QC(idx2))) == 1 && length(idx2) > 1
                % truck only went to this
                % team
                teamMean(truckCnt) = mean(diff(sort(loadTime(idx2))))*24; % hrs
            end
            truckCnt = truckCnt + 1;
        end
        teamIdx = find(teams == j);
        teamDelay(teamIdx) = nanmean(teamMean);
        teamHaul(teamIdx) = mean(haulMi(idx));
    end
    
    figure(1), scatter(i + 0.5*rand(1,length(todaysTeams)), ...
                    teamDelay(match(teams, todaysTeams)), 20, ...
                       teamHaul(match(teams, todaysTeams)), 'filled');
    hold on;
    figure(2), plot(i, length(todaysTeams), 'b.');
    hold on;
end


% ----------------------------------------
% Look at switching behavior
% ----------------------------------------
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
        [time,order] = sort(loadTime(todayData(find(truckId(todayData) ...
                                                    == j))));
        qcs = QC(todayData(find(truckId(todayData) == j)));
        qcs = qcs(order);
        
        switches = find(diff(qcs) ~= 0);
        
        %if length(teamsWorked) == 2     % look only at first switch (easier 
                                        % to characterize)
            QCdayTmp = teamsWorked + 1e6*duration(i);
            teamIdx = match(teams, QCdayTmp)
           
            first = teamIdx(1);
            second = teamIdx(2);
            % NEED TO FIND FIRST, SECOND (AND MAKE SURE NOT
            % SWITCHING BACK AND FORTH
            figure(1), plot(teamDelay(first), ...
                            teamDelay(second), 'b.');
            hold on;

        end
    end
    
end