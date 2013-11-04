% This script calculates the new equality measure: a measure of
% workload inequity between trucks in a team
% 
%  J.Brooks
%  11 Oct 2013
%  
function outData = equalityWorkload(QC,loadTime,loadPercent,truckId,trucks)
    
% team identifier
QCday = QC + 1e6*floor(loadTime);
%teams = outData.QCDay;                  % why restrict ourselves??
teams = unique(QCday);
disp(sprintf('Num unique QCdays: %d', length(teams)));

teamPerfEq = zeros(1,length(teams)); 
teamPerfEq2 = zeros(1,length(teams));
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
    
    for j = 1:n
        loadsHauled(j) = length(find(truckId(idx)==tr(j)));
        effLoadsHauled(j) = sum(loadPercent(find(truckId(idx)==tr(j))));
    end
    
    % New equality measure (Coulter for loads hauled)
    %    teamPerfEq(i) = 1/n*sqrt(sum((loadsHauled./mean(loadsHauled) - 1).^2));

    % New equality measure (Coulter for effective loads hauled)
    teamPerfEq(i) = 1/n*sqrt(sum((effLoadsHauled./mean(effLoadsHauled) ...
                                   - 1).^2));
    
    teamSize(i) = n;
    
    % New effectiveness measure
    teamPerfEff(i) = sum(effLoadsHauled);

    % New efficiency measure
    teamPerfEcy(i) = sum(effLoadsHauled)/n;
end

outData.eq = teamPerfEq;
outData.size = teamSize;
outData.eff = teamPerfEff;
outData.ecy = teamPerfEcy;