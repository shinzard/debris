function [effectiveness, efficiency, equality, sizeF, fluidF, haul, teamHist, dayF] ...
        = fluidity3(dayData,trucks,qcSpan,loadTime,truckId,subcont, ...
                   loadVolume,outData,QC, capacity, dayService, ...
                   VERB,filename,duration,PLOT,towerTime, lat, lon, ...
                    loadPercent, haulMi)
% Run debrisanalysis, generateQcData
disp('Running fluidity code');
%PLOT = 1;
close all;

t = collapseV(dayData, 'QCpercent');
d = collapseV(dayData, 'day');

% Look system-level fluidity: percentage of trucks switching on a given day
%idx = find(~isnan(d));
%days = unique(d(idx)); %[min(d):max(d)]

% 2 June 2012....restrict to same region
days = duration;

%%%%%for i = 1:length(days)
%%%%%    idx = find(d == days(i));
%%%%%    disp(sprintf('Day: %d, NumTrucks: %d', days(i), length(idx)));
%%%%%    fluid(i) = length(find(t(idx)<1))/length(idx);
%%%%%end

%figure, plot(fluid(1:end-1))
%xlabel('Days')
%ylabel('Percentage')
%title('Fluidity: Percentage of Trucks Switching DRT')

% Look at number of changes
%%%%%change = zeros(1, length(trucks));
%%%%%for i = 1:length(dayData)
%%%%%    change(i) = length(find(dayData(i).QCpercent < 1))/ ...
%%%%%        length(dayData(i).QCpercent);
%%%%%    
%%%%%end
%%%%%figure, hist(change, 20);

% look at number of teams for each truck
%%%%%numTeams = zeros(1, length(trucks));
%%%%%for i = 1:length(trucks)
%%%%%   idx =  find(truckId == trucks(i));
%%%%%   numTeams(i) = length(unique(QCday(idx)));
%%%%%end

%figure, hist(numTeams./dayService)
%xlabel('Number of Crews');
%ylabel('Count of Trucks');
%title('Histogram of Number of Crews Per Day of Service');

% Team history
QCday = QC + 1e6*floor(loadTime);

% 30 July 2012: Use data in outData from equality.m
% instead...already filtered
teams = outData.QCDay;
teamHist = zeros(1,length(teams))*NaN; 
teamDay = zeros(1,length(teams)); 
teamSize = zeros(1,length(teams)); 
teamMaxTruck = zeros(1,length(teams)); 
teamPerfEcyNEW = zeros(1,length(teams)); 
teamPerfEcy = zeros(1,length(teams)); 
teamPerfEff = zeros(1,length(teams)); 
teamPerfEq = zeros(1,length(teams)); 
teamHaul = zeros(1,length(teams)); 
teamCapacity = zeros(1,length(teams)); 
teamContractors = zeros(1,length(teams));
teamNestedness = zeros(1,length(teams));
anyInvalidTimes = zeros(1,length(teams));
anyInvalidCap = zeros(1,length(teams));
invalidTime = 0;
invalidCap = 0;
gt2invTickets = 0;
gt2 = 0;
perflCorrected = 0;

for i = 1:length(teams)
    % Find tickets for current team
    idx = find(QCday == teams(i));

    % Find trucks in current team
    [tr,uTrIdx, tmp] = unique(truckId(idx));
    trIdx = match(trucks, tr);
    
    teamDay(i) = floor(loadTime(idx(1)));

    %%%%%    dayIdx = find(days == teamDay(i));

    % Calculate nestedness -- currently unused...
    subs = unique(subcont(idx));
    teamContractors(i) = length(subs);
    
    try
        if length(subs) > 1
            [subCount, tmp] = hist(subcont(idx(uTrIdx)), subs);
            teamNestedness(i) = max(subCount)/length(tr);
            if teamNestedness(i) == 1
                subcont(uTrIdx)
                subs
                subCount
            end
        else
            teamNestedness(i) = 1;
        end
    catch
    end
    
    % Team size/max capacity (only want <100cy)
    n = length(tr);
    teamSize(i) = n;

    teamMaxTruck(i) = max(capacity(trIdx));
    
    eqDataTeamIdx = find(outData.QCDay == teams(i));
    if (n ~= outData.numTrucks(eqDataTeamIdx))
        disp(sprintf('WARNING: Inconsistent team size (%d)', i));
    end
    
    % Calculate performance (effectiveness)
    teamPerfEff(i) = sum(loadVolume(idx));
    
    % Calculate performance (efficiency)
    perfLs = [];
    numTrips = [];
    ahaul = [];

    curQC = teams(i) - 1e6*teamDay(i);
    %    idx = find(floor(loadTime)==teamDay(i) && QC == curQC);
    if n >= 10
        %        getMap = input(sprintf('Plot Map? (%d):', i),
        %        's');
        getMap = 'n';
        if strcmp(getMap, 'y')
            figure;
            color = 'rgbkymc';
            ctr = 1;
        end
    end
    
    for j = trIdx'
        dIdx = find(dayData(j).day == teamDay(i));
        if length(dIdx) > 1
            disp('WARNING...');
        end
        %        perfLs = [perfLs, dayData(j).perfL(dIdx)];
        % 28 Aug 2012: moved to later on after we try to fix zero entries.

        % Investigate large teams
        if n >= 10 & strcmp(getMap, 'y')
            idx = find(floor(loadTime)==teamDay(i)&truckId== ...
                       trucks(j));
            idx2 = find(floor(loadTime)==teamDay(i)&truckId== ...
                       trucks(j)&QC==curQC);
            [tmp,sIdx] = sort(loadTime(idx));
            [tmp,sIdx2] = sort(loadTime(idx2));            
            %            subplot(211);
            plot(lon(idx2(sIdx2)), lat(idx2(sIdx2)), color(ctr));
            hold on;
            plot(lon(idx2(sIdx2)), lat(idx2(sIdx2)), [color(ctr),'.']);
            title(sprintf('QC: %d, Day: %d', curQC, teamDay(i)));

            %            subplot(212);
            %            plot(lat(idx(sIdx)), lon(idx(sIdx)), color(ctr));
            %            hold on;
            %            plot(lat(idx(sIdx)), lon(idx(sIdx)), [color(ctr),'.']);
            %            title('All QC data for same day');
            ctr = ctr + 1;
            if ctr > 7
                ctr = 1;
            end
            hold on;
        end
        
        % Investigate/correct invalid time sequences
        if (dayData(j).perfL(dIdx) == 0)
            % get truck's tickets
            idxL = find(floor(loadTime)==teamDay(i)&truckId== ...
                       trucks(j));
            idxT = find(floor(towerTime)==teamDay(i)&truckId== ...
                       trucks(j));            
        
            % Loads picked up and delivered same day
            idxC = intersect(idxL, idxT);

            REMOVE_INDIV = 0;           % for some reason, this
                                        % makes less available data...
            if REMOVE_INDIV
                method = 0;
                % remove individually invalid tickets
                idxC = removerows(idxC, find(loadTime(idxC)> ...
                                             towerTime(idxC)));
                
                
                tmp = [];
                % build remaining time sequence
                for l = 1:length(idxC)
                    tmp = [tmp, loadTime(idxC(l)), ...
                           towerTime(idxC(l))];
                end
                
                % check if this fixes the sequence...
                valid = issorted(tmp);
            else
                valid = false;
            end
            
            if length(idxC) > 1 & ~valid % if not, try removing one
                                         % additional ticket at random
                method = 1;
                
                % look at random subset of various sizes
                removeOrder = randperm(length(idxC));
                valid = false;
                
                % try to remove one ticket...
                for tmpCtr = 1:length(idxC)
                    tmp = [];
                    removeTicketIdx = removeOrder(tmpCtr);
                    % build time sequence
                    for l = 1:length(idxC)
                        if l ~= removeTicketIdx
                            tmp = [tmp, loadTime(idxC(l)), ...
                                   towerTime(idxC(l))];
                        end
                    end

                    valid = issorted(tmp);
                    if valid
                        break;
                    end
                end

                % try two...
                if length(idxC) > 2 & ~valid
                    method = 2;
                    % look at random subset of various sizes
                    removeOrder1 = randperm(length(idxC));
                    removeOrder2 = randperm(length(idxC));
                    valid = false;      % redundant
                
                    % try to remove two tickets...
                    for tmpCtr1 = 1:length(idxC)
                        for tmpCtr2 = 1:length(idxC)
                            tmp = [];
                            removeTicketIdx1 = removeOrder1(tmpCtr1);
                            removeTicketIdx2 = removeOrder2(tmpCtr2);
                            % build time sequence
                            for l = 1:length(idxC)
                                if l ~= removeTicketIdx1 && ...
                                   l ~= removeTicketIdx2
                                    tmp = [tmp, loadTime(idxC(l)), ...
                                           towerTime(idxC(l))];
                                end
                            end
                        
                            valid = issorted(tmp);

                            % end when found
                            if valid
                                break;
                            end
                        end

                        % end when found
                        if valid
                            break;
                        end
                    end
                end                     % end try two
            else                        % only one ticket and
                                        % towerTime < loadTime!
                valid = false;
            end % length > 1 & ~valid
            
            if valid                    % correction worked; recalc perfL
                if method == 0          % do nothing...only had to remove 
                                        % individually invalid tickets
                elseif method == 1          % only one ticket removed
                    idxC = removerows(idxC, removeTicketIdx);
                else                    % two tickets removed
                    idxC = removerows(idxC, [removeTicketIdx1, ...
                                        removeTicketIdx2] );
                end

                dayData(j).perfL(dIdx) = ...
                    sum(loadPercent(idxC).*haulMi(idxC)./ ...
                        ((towerTime(idxC) - ...
                          loadTime(idxC))*24));
                perflCorrected = perflCorrected + 1;
                if dayData(j).perfL(dIdx) < 0 || any(loadTime(idxC)> ...
                                                     towerTime(idxC))
                    loadPercent(idxC)
                    haulMi(idxC)
                    towerTime(idxC) - loadTime(idxC)
                    error('WRONG TIME SEQUENCE');
                end
                    
            else
                gt2 = gt2 + 1;          % need to remove > 2
                if length(find(loadTime(idxC) > towerTime(idxC))) > 2
                    %disp(sprintf('%d INVALID TICKET(S)', ... 
                    %             length(find(loadTime(idxC) > ...
                    %                         towerTime(idxC)))));
                    gt2invTickets = gt2invTickets + 1;
                end
            end
%            figure, plot(loadTime(idx), 'b.');
        %            hold on, plot(towerTime(idx), 'r.')
        %            title(sprintf('Truck: %d, Day: %d', j, teamDay(i)));
        end                             % If perfL == 0

        perfLs = [perfLs, dayData(j).perfL(dIdx)];
        numTrips = [numTrips, dayData(j).numhauls(dIdx)];
        ahaul = [ahaul, dayData(j).ahaul(dIdx)];
    end
    
    if n >= 10 & strcmp(getMap, 'y');
        plot_google_map;                % can only do ~400/day!
    end
    
    if any(perfLs == 0)
        disp(sprintf(['WARNING: some invalid time sequences for ', ...
                         'team %d; THROWING AWAY'],i)); 
        % 6 August 2012: Throw away this data...
        anyInvalidTimes(i) = 1;

        if teamSize(i) > 1
            invalidTime = invalidTime + 1;
            %            perfLs
            %            teams(i) - 1e6*teamDay(i)
        end
    end
    
    if any(numTrips == 0)
        disp(sprintf(['WARNING: some numTrips = 0 for ', ...
                         'team %d'],i)); 
    end
    
    if any(ahaul == 0)
        disp(sprintf(['WARNING: some ahaul = 0 for ', ...
                         'team %d'],i)); 
    end
    
    if PLOT & VERB
        figure(10), plot(perfLs, 'b.'), hold on, plot(numTrips, 'r.')
        axis([-1, n + 1, -10, max(perfLs)+10]);
        title(sprintf('QCDay: %d, TeamSize, : %d',teams(i),n));
        pause(2);
        hold off;
    end
    
    % HOW MANY NAN'S ARE WE THROWING AWAY??
    if sum(isnan(perfLs) | isnan(numTrips))> 1
        disp(sprintf('NaNs: %d', sum(isnan(perfLs)|isnan(numTrips))> ...
                     1));
        disp(sprintf('Team Size: %d', teamSize(i)));
    end

    % Take into account total capacity
    % 28 July 2012: This is also in outData -- should compare
    teamCapacity(i) = sum(capacity(trIdx));
    if any(capacity(trIdx)==0)
        anyInvalidCap(i) = 1;
        disp(sprintf('Team %d contains truck with zero capacity.', i));
        if teamSize(i) > 0
            invalidCap = invalidCap + 1;
        end
    end

    eqDataTeamIdx = find(outData.QCDay == teams(i));
    if (teamCapacity(i) ~= outData.capability(eqDataTeamIdx))
        disp(sprintf('WARNING: Inconsistent team capacity (%d)', i));
    end
    
    teamPerfEcy(i) = mean(perfLs./numTrips);
    
    if isnan(teamPerfEcy(i))
        disp(sprintf('NaN efficiency for team %d', i));
    end
    %teamPerfEcy(i) = teamPerfEff(i)/sum(numTrips);
    %teamPerfEcy(i) = teamPerfEff(i)/teamCapacity(i);
    teamPerfEcyNEW(i) = sum(ahaul.*numTrips)/teamSize(i);
    
    % Equality extraction
    try
        teamPerfEq(i) = outData.remaining(find(outData.QCDay== ...
                                               teams(i)));
        teamHaul(i) =  outData.ahaul(find(outData.QCDay== ...
                                               teams(i)));
    catch % this shouldn't happen anymore (30 July 2012)
        teamPerfEq(i) = NaN;
    end
    

end

%%new section for familiarity calculation computer cumulative sum
%%on the fly
%days = unique(floor(loadTime)); %use days vector from above

% Look at total work history
cumhistory = zeros(length(trucks), length(trucks));

for i = 1:length(days)
    history = zeros(length(trucks), length(trucks));
    idx = find(floor(loadTime) == days(i));
    todayQC = unique(QC(idx));
    
    % Each iteration through this loop is another team-day
    for j = 1:length(todayQC)
        % NEW 30 July 2012
        thisQCday = todayQC(j) + 1e6*days(i);
        if isempty(find(outData.QCDay==thisQCday))
            disp(sprintf('Team already filtered out: %d', ...
                         thisQCday));
            %            break; %28 Aug 2012, should skip???
        end

        teamTrucks = unique(truckId(idx(find(QC(idx)== ...
                                             todayQC(j)))));
        % UPDATE n!! 31 July 2012
        n = length(teamTrucks);
        tmp = match(trucks,teamTrucks);

        teamCumHistory = 0;
        
        for k = 1:length(teamTrucks)-1
            history(tmp(k), tmp(k+1:end)) = 1;
            teamCumHistory = teamCumHistory + ...
                sum(cumhistory(tmp(k),tmp(k+1:end)));
        end
        
        % Calculate team familiarity (average dyad working history)
        teamIdx = find(teams == todayQC(j)+1e6*days(i));

        % Added 30 July 2012 (if + max truck only)
        if ~isempty(teamIdx)
            if ( (length(tmp) == 1) | ...
                 (teamCapacity(teamIdx)/teamSize(teamIdx) >= 80) | ...  % teamMaxTruck(teamIdx) >= 80 |...
                 (teamHaul(teamIdx) >= 30) | ...
                 (anyInvalidTimes(teamIdx))| ... % Added 6 aug 2012
                 (anyInvalidCap(teamIdx)) ) % added 13 aug 2012
                teamHist(teamIdx) = NaN;
                teamDay(teamIdx) = NaN; % to prevent single data points
                                        % from messing up the analysis
            else
                teamHist(teamIdx) = teamCumHistory;

                % Scale by number of pairs between members
                teamHist(teamIdx) = teamHist(teamIdx)/((n*(n-1))/2);
            end
        end
    end

    % Accumulate today's history
    cumhistory = cumhistory + history;
end

% Final history total days
history = cumhistory;

if PLOT
    figure, image(history)
end

% Look at final distribution of work history pairs
if PLOT
    id = find(history > 0);
    figure, hist(history(id), max(history(id)))
    title('Histogram of Truck Pair History')
    xlabel('Number of Days Working History')
    ylabel('Count');
end
    
% Average relationship length
numPartners = zeros(1, length(trucks));
for i = 1:length(trucks)
    numPartners(i) = length(find(history(i,:)));
    if isempty(numPartners(i))
        numPartners(i) = NaN;
    end
end
avgLen = sum(history')./numPartners./dayService;

if PLOT
    figure, hist(avgLen)
    title('Histogram of Average Relationship Length');
    xlabel('Average Length (normalized to active days)');
    ylabel('Count of Trucks');
    
    figure, plot(dayService, avgLen, 'b.')
    title('Avg. Relationship Length by Active Days')
    ylabel('Average relationship length (norm)');
    xlabel('Days of Active Service')
    
    medLen = ones(1,max(dayService));
    for i = 2:max(dayService)
        idx = find(dayService == i);
        medLen(i) = nanmedian(avgLen(idx));
    end
    hold on, plot(medLen, 'r')
end

    
%% END OF NEW SECTION
    
% May also want to do this mean by day...or median

if PLOT
    teamPerfEffRatio = teamPerfEff/mean(teamPerfEff);
    teamPerfEffDiff = teamPerfEff - mean(teamPerfEff);
    teamPerfEcyRatio = teamPerfEcy/mean(teamPerfEcy);
    
    figure, boxplot(teamHist, teamDay); 
    title('Average working history within teams');
    
    figure, boxplot(teamHist./(teamDay-min(teamDay)-10), teamDay); 
    title('Average working history within teams; Normalized');
    
    figure, boxplot(teamSize, teamDay);
    title('Team size')
    
    figure, boxplot(teamHist, teamSize);
    title('Team History by Team Size')
    xlabel('Team Size')
    ylabel('Average Days Pairwise Working History')
    
    figure, boxplot(teamPerfEffDiff, teamDay);
    title('Performance: Effectiveness Difference');
    
    figure, boxplot(teamPerfEffRatio, teamDay);
    title('Performance: Effectiveness Ratio');
    
    figure, subplot(121), boxplot(teamPerfEffRatio, teamSize)
    title('Effectiveness')
    xlabel('Team Size')
    subplot(122), boxplot(teamPerfEcyRatio, teamSize)
    title('Efficiency')
    xlabel('Team Size')
    
end

%% NEW 30 July 2012
disp('Recoding Section');

% Recoded data;
%idx = find(teamSize>1); % cannot compare fluidity of a single
                        % truck...

% 30 July 2012: better, since we now have additional conditions...
idx = find(~isnan(teamHist));

%%%%%codeLimits = quantile(teamSize(idx), [0.25, 0.75]);
%%%%%
%%%%%teamSizeCode = (teamSize(idx) <= codeLimits(1)) ...
%%%%%    + 2*(teamSize(idx) > codeLimits(1)) ...
%%%%%    + (teamSize(idx) > codeLimits(2) );

% Scale familiarity (working history) to number of days in mission
% 28 July 2012: Should we scale this??
%teamHist2 = teamHist./(teamDay-min(teamDay)+1);
% moved to below calc of dayF

% this (may) cause large differences between peak analysis and complete mission

%%%%%codeLimits = quantile(teamHist2(idx), [0.25, 0.75]);
%%%%%
%%%%%fluidCode = (teamHist2(idx) <= codeLimits(1)) ...
%%%%%    + 2*(teamHist2(idx) > codeLimits(1)) ...
%%%%%    + (teamHist2(idx) > codeLimits(2) );
%%%%%
%%%%%idx2 = find(fluidCode ~= 2 & teamSizeCode~=2);

% Data is log-normal!

%disp('Efficiency ANOVA')
%anovan(log10(teamPerfEcy(idx(idx2))), ...
%       {fluidCode(idx2)' teamSizeCode(idx2)'}, 'model', 'interaction', ...
%       'varnames', {'Fluidity', 'Team Size'})

%disp('Effectiveness ANOVA')
%anovan(log10(teamPerfEff(idx(idx2))), {fluidCode(idx2)' teamSizeCode(idx2)'}, ...
%       'model', 'interaction', 'varnames', {'Fluidity', ['Team ' ...
%                    'Size']})


% 28 July 2012: Why are we throwing away so much data??? Is it the
% log10 transform...what do we gain by having that guy?
good = find(~isnan(log10(teamPerfEcy(idx))) & ... 
            ~isinf(log10(teamPerfEcy(idx))) & ...
            ~isnan(log10(teamPerfEq(idx)*100000)) & ...
            ~isinf(log10(teamPerfEq(idx)*100000)) & ...
            ~isnan(log10(teamPerfEff(idx))) & ...
            ~isinf(log10(teamPerfEff(idx))));
bad = find(isnan(log10(teamPerfEcy(idx))) | ... 
           isinf(log10(teamPerfEcy(idx))) | ...
           isnan(log10(teamPerfEq(idx)*100000)) | ...
           isinf(log10(teamPerfEq(idx)*100000)) | ...
           isnan(log10(teamPerfEff(idx))) | ...
           isinf(log10(teamPerfEff(idx))));

figure, plot(teamPerfEff(idx), 'b.'), hold on;
plot(teamPerfEcy(idx), 'r.');
plot(teamPerfEq(idx)*100000, 'g.');
legend({'Effectiveness', 'Efficiency', 'Equality'});

teamPerfEcy(idx(bad))
teamPerfEff(idx(bad))
teamPerfEq(idx(bad))

disp(sprintf('Good Ratio: %2.2f', length(good)/length(idx)));


effectiveness = log10(teamPerfEff(idx(good)));
efficiency = log10(teamPerfEcy(idx(good)));
equality = log10(teamPerfEq(idx(good))*100000); % per 100k cubic yards

%fluidF = teamHist2(idx(good));
fluidF = teamHist(idx(good));
sizeF = teamSize(idx(good));
haul = teamHaul(idx(good));
dayF = teamDay(idx(good))-days(1)+1;
teamHist = teamHist(idx(good));

% Scale familiarity (working history) to number of days in mission
% 28 July 2012: Should we scale this?? (USED TO BE CALLED TEAMHIST2)
fluidF = fluidF./(dayF-min(dayF)+1);

if min(dayF) ~= 1
    disp('Minimum teamDay is not 1');
end

if VERB
    figure, subplot(331), histfit(effectiveness);
    title('Effectiveness (log_{10} transformed)');
    subplot(332), plot(efficiency, effectiveness, 'b.');
    subplot(333), plot(equality, effectiveness, 'b.');
    
    subplot(334), plot(effectiveness, efficiency, 'b.');
    subplot(335), histfit(efficiency);
    title('Efficiency (log_{10} transformed)');
    subplot(336), plot(equality, efficiency, 'b.');
    
    subplot(337), plot(effectiveness, equality, 'b.');
    subplot(338), plot(efficiency, equality, 'b.');
    subplot(339), histfit(equality);
    title('Equality (log_{10} transformed)')
    
    [h,p] = lillietest(efficiency)
    [h,p] = lillietest(effectiveness)
    [h,p] = lillietest(equality)
    
    figure, qqplot(efficiency)
    title('Efficiency')
    figure, qqplot(effectiveness)
    title('Effectiveness')
    figure, qqplot(equality)
    title('Equality')
end

%%%%%codeLimit = quantile(sizeF, [0.5]);
%%%%%sizeCode = sizeF > codeLimit;
%%%%%codeLimit = quantile(fluidF, [0.5]);
%%%%%fluidCode = fluidF > codeLimit;

if PLOT
    figure, hist(sizeCode), title('Size');
    figure, hist(fluidCode), title('Fluid');
end

if VERB
    disp('Efficiency');
    [h,p] = lillietest(efficiency(find(sizeCode==0&fluidCode==0)))
    [h,p] = lillietest(efficiency(find(sizeCode==0&fluidCode==1)))
    [h,p] = lillietest(efficiency(find(sizeCode==1&fluidCode==0)))
    [h,p] = lillietest(efficiency(find(sizeCode==1&fluidCode==1)))
    
    disp('Effectiveness');
    [h,p] = lillietest(effectiveness(find(sizeCode==0&fluidCode==0)))
    [h,p] = lillietest(effectiveness(find(sizeCode==0&fluidCode==1)))
    [h,p] = lillietest(effectiveness(find(sizeCode==1&fluidCode==0)))
    [h,p] = lillietest(effectiveness(find(sizeCode==1&fluidCode==1)))
end

csvwrite(filename, [effectiveness', efficiency', equality', teamCapacity(idx(good))', sizeF', fluidF', haul',dayF']);

disp(sprintf(['Teams with at least one truck with invalid time ' ...
              'sequence: %d'], invalidTime));
disp(sprintf('Teams with at least one truck with zero capacity: %d', ...
             invalidCap));

disp(sprintf('More complicated cases: %d', gt2));
disp(sprintf('Of these, %d, have >2 invalid tickets', ...
             gt2invTickets))

disp(sprintf('Corrected perfLs: %d', perflCorrected));
%teamContractors
%teamMaxTruck