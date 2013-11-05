% This script calculates team familiarity metric based on prior
% working history of team members (i.e., trucks); extracted from
% fluidity3.m 
% 
%  J.Brooks
%  21 Oct 2013
%  Last Modified: 21 Oct 2013
%  
function [teamHist,teamHaul,teamCapacity] = familiarity(tickets,teams,duration)
disp('Calculating team familiariy measure...');

% Unpack ticket data
loadTime = tickets.loadTime;
QC = tickets.QC;
truckId = tickets.truckId;
trucks = tickets.trucks;
capacity = tickets.capacity;
haulMi = tickets.haulMi;
QCday = QC + 1e6*floor(loadTime);

% Allocate memory
teamCapacity = zeros(1,length(teams));
teamHaul = zeros(1,length(teams));
teamHist = zeros(1,length(teams));

%days = sort(unique(floor(loadTime)));   % double check sorted!
days = sort(duration);

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
        teamIdx = find(teams == thisQCday);

        if isempty(teamIdx)
            disp(sprintf('Team already filtered out: %d', ...
                         thisQCday));
            % continue to include history even though filtered...
        end

        teamTrucks = unique(truckId(idx(find(QC(idx)== ...
                                             todayQC(j)))));
        % UPDATE n!! 31 July 2012
        n = length(teamTrucks);
        trIdx = match(trucks,teamTrucks);

        teamCumHistory = 0;
        
        for k = 1:length(teamTrucks)-1
            history(trIdx(k), trIdx(k+1:end)) = 1;
            teamCumHistory = teamCumHistory + ...
                sum(cumhistory(trIdx(k),trIdx(k+1:end)));
        end
        
        % Other team-level vars
        teamCapacity(teamIdx) = max(capacity(trIdx));
        teamHaul(teamIdx) = nanmean(haulMi(find(QCday==thisQCday)));
        
        % Calculate team familiarity (average dyad working history)
        if ~isempty(teamIdx)
            if ( (length(trIdx) == 1) | ... % single truck 'teams'
                 (teamCapacity(teamIdx) >= 90) | ...% large trucks
                 (teamHaul(teamIdx) >= 30) ) % team with long hauls
                 % THESE ARE NO LONGER RELEVANT
                 %(anyInvalidTimes(teamIdx))| ... % Added 6 aug 2012
                 %(anyInvalidCap(teamIdx)) ) % added 13 aug 2012

                %teamHist(teamIdx) = teamCumHistory/((n*(n-1))/2);
                teamHist(teamIdx) = NaN;% to prevent these from
                                        % messing up the analysis
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
