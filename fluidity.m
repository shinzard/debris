function [effectiveness, efficiency, equality, sizeF, fluidF, haul] ...
        = fluidity(dayData,trucks,qcSpan,loadTime,truckId,subcont, ...
                   loadVolume,outData,VERB,filename)
    
% Run debrisanalysis, generateQcData
disp('Running fluidity code');
PLOT = 0;
close all;

t = collapseV(dayData, 'QCpercent');
d = collapseV(dayData, 'day');

% Look system-level fluidity: percentage of trucks switching on a given day
idx = find(~isnan(d));
days = unique(d(idx)); %[min(d):max(d)]

for i = 1:length(days)
    idx = find(d == days(i));
    disp(sprintf('Day: %d, NumTrucks: %d', days(i), length(idx)));
    fluid(i) = length(find(t(idx)<1))/length(idx);
end

%figure, plot(fluid(1:end-1))
%xlabel('Days')
%ylabel('Percentage')
%title('Fluidity: Percentage of Trucks Switching DRT')

% Look at number of changes
change = zeros(1, length(trucks));
for i = 1:length(dayData)
    change(i) = length(find(dayData(i).QCpercent < 1))/ ...
        length(dayData(i).QCpercent);
    
end
figure, hist(change, 20);

% Look at total work history
history = zeros(length(trucks), length(trucks), length(days));
for i = 1:length(qcSpan)
    for j = 1:length(qcSpan(i).day)
        tmp = find(qcSpan(i).trucks(:,j));
        for k = 1:length(tmp)-1
            % is this double counting pairs? -- YES!
            %history(tmp(k), tmp(k+1:end)) = history(tmp(k), tmp(k+1:end)) + 1;
            currDay = find(days == qcSpan(i).day(j));
            history(tmp(k), tmp(k+1:end), currDay) = 1;
        end
    end
end
cumhistory = cumsum(history, 3);
history = sum(history, 3);

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
%avgLen = sum(history')./numPartners./dayService;

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

% Look at number of teams for each truck
numTeams = zeros(1, length(trucks));
QCday = QC + 1e6*floor(loadTime);
for i = 1:length(trucks)
   idx =  find(truckId == trucks(i));
   numTeams(i) = length(unique(QCday(idx)));
end

%figure, hist(numTeams./dayService)
%xlabel('Number of Crews');
%ylabel('Count of Trucks');
%title('Histogram of Number of Crews Per Day of Service');

% Team history
teams = unique(QCday);
teamHist = zeros(1,length(teams)); 
teamDay = zeros(1,length(teams)); 
teamSize = zeros(1,length(teams)); 
teamPerfEcyOLD = zeros(1,length(teams)); 
teamPerfEcy = zeros(1,length(teams)); 
teamPerfEff = zeros(1,length(teams)); 
teamPerfEq = zeros(1,length(teams)); 
teamHaul = zeros(1,length(teams)); 
teamCapacity = zeros(1,length(teams)); 
teamContractors = zeros(1,length(teams));
teamNestedness = zeros(1,length(teams));
for i = 1:length(teams)
    idx = find(QCday == teams(i));
    [tr,uTrIdx, tmp] = unique(truckId(idx));
    trIdx = match(trucks, tr);
    teamDay(i) = floor(loadTime(idx(1)));
    dayIdx = find(days == teamDay(i));

    % Calculate nestedness
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
    
    % Team size
    n = length(tr);
    teamSize(i) = n;
    
    % Calculate performance (effectiveness)
    teamPerfEff(i) = sum(loadVolume(idx));
    
    % Calculate performance (efficiency)
    perfLs = [];
    numTrips = [];
    ahaul = [];
    
    for j = trIdx'
        dIdx = find(dayData(j).day == teamDay(i));
        if length(dIdx) > 1
            disp('WARNING...');
        end
        perfLs = [perfLs, dayData(j).perfL(dIdx)];
        numTrips = [numTrips, dayData(j).numhauls(dIdx)];
        ahaul = [ahaul, dayData(j).ahaul(dIdx)];
    end
    
    % HOW MANY NAN'S ARE WE THROWING AWAY??
    if sum(isnan(perfLs) | isnan(numTrips))> 1
        disp(sprintf('NaNs: %d', sum(isnan(perfLs)||isnan(numTrips))> ...
                     1));
        disp(sprintf('Team Size: %d', teamSize(i)));
    end

    % Take into account total capacity
    teamCapacity(i) = sum(capacity(trIdx));

    teamPerfEcy(i) = nanmean(perfLs./numTrips);
    %teamPerfEcy(i) = teamPerfEff(i)/sum(numTrips);
    %teamPerfEcy(i) = teamPerfEff(i)/teamCapacity(i);
    teamPerfEcyNEW(i) = sum(ahaul.*numTrips)/teamSize(i);
    
    % Equality extraction
    try
        teamPerfEq(i) = outData.remaining(find(outData.QCDay== ...
                                               teams(i)));
        teamHaul(i) =  outData.ahaul(find(outData.QCDay== ...
                                               teams(i)));
    catch
        teamPerfEq(i) = NaN;
    end
    

    % Calculate team familiarity (average dyad working history)
    if length(tr) == 1
        teamHist(i) = NaN;
    elseif dayIdx > 1
        for j = 1:length(tr)
            for k = j:length(tr)
                teamHist(i) = teamHist(i) +  cumhistory(trIdx(j), ...
                                                     trIdx(k),dayIdx-1);
            end
        end
        % Scale by number of pairs between members
        teamHist(i) = teamHist(i)/((n*(n-1))/2);
    end
end

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

% Recoded data;
idx = find(teamSize>1); % cannot compare fluidity of a single truck...
codeLimits = quantile(teamSize(idx), [0.25, 0.75]);

teamSizeCode = (teamSize(idx) <= codeLimits(1)) ...
    + 2*(teamSize(idx) > codeLimits(1)) ...
    + (teamSize(idx) > codeLimits(2) );

% Scale familiarity (working history) to number of days in mission
teamHist2 = teamHist./(teamDay-min(teamDay)+1);
codeLimits = quantile(teamHist2(idx), [0.25, 0.75]);

fluidCode = (teamHist2(idx) <= codeLimits(1)) ...
    + 2*(teamHist2(idx) > codeLimits(1)) ...
    + (teamHist2(idx) > codeLimits(2) );

idx2 = find(fluidCode ~= 2 & teamSizeCode~=2);

% Data is log-normal!

%disp('Efficiency ANOVA')
%anovan(log10(teamPerfEcy(idx(idx2))), ...
%       {fluidCode(idx2)' teamSizeCode(idx2)'}, 'model', 'interaction', ...
%       'varnames', {'Fluidity', 'Team Size'})

%disp('Effectiveness ANOVA')
%anovan(log10(teamPerfEff(idx(idx2))), {fluidCode(idx2)' teamSizeCode(idx2)'}, ...
%       'model', 'interaction', 'varnames', {'Fluidity', ['Team ' ...
%                    'Size']})

good = find(~isnan(log10(teamPerfEcy(idx))) & ... 
            ~isinf(log10(teamPerfEcy(idx))) & ...
            ~isnan(log10(teamPerfEq(idx))) & ...
            ~isinf(log10(teamPerfEq(idx))) & ...
            ~isnan(log10(teamPerfEff(idx))) & ...
            ~isinf(log10(teamPerfEff(idx))));

disp(sprintf('Good Ratio: %2.2f', length(good)/length(idx)));

effectiveness = log10(teamPerfEff(idx(good)));
efficiency = log10(teamPerfEcy(idx(good)));
equality = log10(teamPerfEq(idx(good))*100000); % per 100k cubic yards
fluidF = teamHist2(idx(good));
sizeF = teamSize(idx(good));
haul = teamHaul(idx(good));

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

codeLimit = quantile(sizeF, [0.5]);
sizeCode = sizeF > codeLimit;
codeLimit = quantile(fluidF, [0.5]);
fluidCode = fluidF > codeLimit;

figure, hist(sizeCode), title('Size');
figure, hist(fluidCode), title('Fluid');

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

csvwrite(filename, [effectiveness', efficiency', equality', teamCapacity(good)', sizeF', fluidF', haul']);
