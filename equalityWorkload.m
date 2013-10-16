% This script calculates the new equality measure: a measure of
% workload inequity between trucks in a team
% 
%  J.Brooks
%  11 Oct 2013
%  
function outData = equalityWorkload(QC,loadTime,truckId,trucks,outData)
    
% team identifier
QCday = QC + 1e6*floor(loadTime);
teams = outData.QCDay;

teamPerfEq = zeros(1,length(teams)); 

for i = 1:length(teams)
    % Find tickets for current team
    idx = find(QCday == teams(i));

    % Find trucks in current team
    [tr,uTrIdx, tmp] = unique(truckId(idx));
    trIdx = match(trucks, tr);
    
    n = length(tr);
    
    loadsHauled = zeros(1,n);
    
    for j = 1:n
        loadsHauled(j) = length(find(truckId(idx)==tr(j)));
    end
    
    teamPerfEq(i) = 1/n*sqrt(sum((loadsHauled./mean(loadsHauled) - 1).^2));

end

outData.eq = teamPerfEq;
outData.loads = loadsHauled