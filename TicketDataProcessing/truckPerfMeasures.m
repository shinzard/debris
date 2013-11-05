% This script calculates effectiveness performance data at the
% truck level (used for hierarchical Null model for IMME proposal
% 
%  J.Brooks
%  30 Oct 2013
%  Last Modified: 30 Oct 2013
%  
function [outData] = truckPerfMeasures(tickets)
disp('Calculating truck performance measures...');

% Unpack ticket data
QC = tickets.QC;    
loadPercent = tickets.loadPercent;
truckId = tickets.truckId;
trucks = tickets.trucks;
day = floor(tickets.loadTime);
days = unique(floor(tickets.loadTime));
sub = match(unique(tickets.subcont), tickets.subcont);

% Allocate memory for measures
id = [];
qc = [];
d = [];
eff = [];
nL = [];
sc = [];

for i = 1:length(trucks)
    for j = 13:length(days)             % ignore first 12 days
        % Find tickets for current truck and day
        idx = find(truckId == trucks(i) & day == days(j));
        
        QCs = unique(QC(idx));

        % Split out loads by team
        
        %        for k = 1:length(QCs)
        k = 1;
        if ( length(QCs) == 1 )         % only keep trucks working
                                        % consistently for one team
            loads = idx(find(QC(idx)==QCs(k)));
            id = [id, trucks(i)];
            qc = [qc, QCs(k)];
            d = [d, days(j)];
            nL = [nL, length(loads)];
            eff = [eff, sum(loadPercent(loads))];
            sc = [sc, sub(idx(1))];     % NOTE assuming all same
                                        % subcontractor...
            if length(unique(sub(idx)))>1
                warning('more than one sub');
            end
        end
    end
end

outData.id = id;
outData.qc = qc;
outData.d = d;
outData.eff = eff;
outData.nL = nL;
outData.sc = sc;