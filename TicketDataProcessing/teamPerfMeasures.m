% This script calculates the new equality measure: a measure of
% workload inequity between trucks in a team
% 
%  J.Brooks
%  21 Oct 2013
%  Last Modified: 21 Oct 2013
%  
function [outData,teamSize] = teamPerfMeasures(tickets,teams,trucks)
disp('Calculating team performance measures...');

% Unpack ticket data
QC = tickets.QC;    
loadTime = tickets.loadTime;
loadPercent = tickets.loadPercent;
truckId = tickets.truckId;
trucks = tickets.trucks;
haulMi = tickets.haulMi;

% team identifier
QCday = QC + 1e6*floor(loadTime);

% Allocate memory for measures
teamPerfEq = zeros(1,length(teams)); 
teamPerfEff = zeros(1,length(teams));
teamPerfEcy = zeros(1,length(teams));
teamSize = zeros(1,length(teams)); 

for i = 1:length(teams)
    % Find tickets for current team
    idx = find(QCday == teams(i));

    % Find trucks in current team
    [tr,uTrIdx, tmp] = unique(truckId(idx));
    trIdx = match(trucks, tr);
    
    n = length(tr);
    
    loadsHauled = zeros(1,n);
    effLoadsHauled = zeros(1,n);
    haul = zeros(1,n);
    
    for j = 1:n
        loads = idx(find(truckId(idx)==tr(j)));
        loadsHauled(j) = length(loads);
        effLoadsHauled(j) = sum(loadPercent(loads));
        haul(j) = nanmean(haulMi(loads));
    end
    
    % New equality measure (Coulter for loads hauled)
    %    teamPerfEq(i) = 1/n*sqrt(sum((loadsHauled./mean(loadsHauled) - 1).^2));

    % New equality measure (Coulter for effective loads hauled)
    teamPerfEq(i) = 1/sqrt(n*n-1)*sqrt(sum((effLoadsHauled./mean(effLoadsHauled) ...
                                   - 1).^2));
    
    teamSize(i) = n;
    
    % New effectiveness measure
    teamPerfEff(i) = sum(effLoadsHauled);

    % New efficiency measure
    teamPerfEcy(i) = sum(effLoadsHauled.*haul)/n;
end

outData.eq = teamPerfEq;
outData.eff = teamPerfEff;
outData.ecy = teamPerfEcy;