function [outData] = ...
        equality(lat,lon,duration,loadTime,QC,truckId,...
                 haulMi,capacity,loadVolume,trucks)
% This function calculates equality measures of performance
%
% 19 January 2012
% Modified: 28 May 2012
% J.Brooks

% LOAD all data into workspace prior to running.
disp('Running equality code');

GRID_SIZE = 0.01;

% 24 May 2012
% overwrite to look only at contiguous data
%duration = [40672:40817];%[40390:40981];
% moved to analysisScript

%nonz = find(lat > 0);
%latMax = max(lat(nonz));
%latMin = min(lat(nonz));
latMin = 32.22;
latMax = 35;
latDim = ceil((latMax-latMin)/GRID_SIZE);

%nonz = find(lon < 0);
%lonMin = min(lon(nonz));
%lonMax = max(lon(nonz));
lonMin = -88.4;
lonMax = -85.4;
lonDim = ceil((lonMax-lonMin)/GRID_SIZE);

latGrid = [latMin: GRID_SIZE: latMax+GRID_SIZE]; 
lonGrid = [lonMin: GRID_SIZE: lonMax+GRID_SIZE];
debrisGrid = zeros(length(latGrid), length(lonGrid), length(duration));

dayCount = length(duration);

outData.day = [];
outData.remaining = [];
outData.QCDay = [];
outData.ahaul = [];
% Added 28 July 2012.
outData.capability = [];
outData.anumTrips = [];
outData.numTrucks = [];


% team identifier
QCday = QC + 1e6*floor(loadTime);
teams = unique(QCday);
disp(sprintf('Num unique QCdays: %d', length(teams)));

for j = 1:length(duration)
    todayData = find(floor(loadTime)==duration(end-dayCount+1) & ...
                     lat~=0  & lon~=0 & lat>latMin & lat<latMax ...
                     & lon>lonMin & lon<lonMax);
    todaysTeams = unique(QCday(todayData));
    teamGrids = cell(1,length(todaysTeams));
    
    %    disp(sprintf('Today: %d\t Num.Teams: %d\n', ... 
    %                 duration(end-dayCount+1), ...
    %                 length(todaysTeams)));
        
    if dayCount < length(duration)
        debrisGrid(:,:,dayCount) = debrisGrid(:,:,dayCount+1);
    end
    
    for i = 1:length(todayData)
        latIdx = find(latGrid < lat(todayData(i)), 1, 'last');
        lonIdx = find(lonGrid < lon(todayData(i)), 1, 'last');
        
        teamIdx = find(todaysTeams == QCday(todayData(i)));

        % Need to account for > 1 team per grid
        gridIdx = sub2ind(size(debrisGrid(:,:,1)), ...
                          latIdx,lonIdx); 

        % This handles >1 grid per team
        teamGrids{teamIdx} = [teamGrids{teamIdx}, gridIdx];
            
        debrisGrid(latIdx, lonIdx, dayCount) = ...
            debrisGrid(latIdx, lonIdx, dayCount) + ...
            loadVolume(i);
    end
    
    todayGrid = debrisGrid(:,:,dayCount);
    capabilities = zeros(1,length(todaysTeams));
    numTrips = zeros(1,length(todaysTeams));
    numTrucks = zeros(1,length(todaysTeams));
    ahaul = zeros(1,length(todaysTeams));
    
    for i = 1:length(todaysTeams)
        idx = find(QCday(todayData) == todaysTeams(i));
        
        % find number of trucks in team
        uniqueTrucks= unique(truckId(todayData(idx)));
        
        % add trucks' capacities
        capabilities(i) = sum(capacity(match(trucks, ...
                                             uniqueTrucks)));

        % Average number of trips made
        anumTrips(i) = length(idx)/length(uniqueTrucks);
        
        % Average haul distance
        tmp = find(isnan(haulMi(todayData(idx)))==0);
        ahaul(i) = mean(haulMi(todayData(idx(tmp))));
        %ahaul(i) = nanmean(haulMi(todayData(idx)));

        % Number of trucks
        numTrucks(i) = length(uniqueTrucks);
    end
    
    eqPerf = zeros(1,length(todaysTeams));
    grids = [];
    for i = 1:length(todaysTeams)
        %eqPerf(i) = sum(todayGrid(unique(teamGrids{i})))/...
        %    capabilities(i)/numTrips(i); 
        eqPerf(i) = numTrucks(i)/sum(todayGrid(unique(teamGrids{i})));

        % 29 Aug 2012: New equality metric idea:
        %eqPerf(i) = capabilities(i)/sum(todayGrid(unique(teamGrids{i})));
        % can make this modification in R script: eq*capacity/size

        grids = [grids, unique(teamGrids{i})];
    end
    
    if length(unique(grids)) < length(grids)
        disp(sprintf('Team overlap: %2.2f', 1-length(unique(grids))/ ...
                     length(grids))); 
    end
    
    outData.QCDay = [outData.QCDay, todaysTeams'];
    outData.remaining = [outData.remaining, eqPerf];
    outData.day = [outData.day, repmat(dayCount, ...
                                       1,length(eqPerf))]; 
    outData.ahaul = [outData.ahaul, ahaul];
    
    % Added 28 July 2012
    outData.capability = [outData.capability, capabilities];
    outData.anumTrips = [outData.anumTrips, anumTrips];
    outData.numTrucks = [outData.numTrucks, numTrucks];

    
    %    figure(10), plot(eqPerf, capabilities, 'b.');
    %    figure(20), plot(numTrips, 'b.'), hold on;
    %    plot(todayGrid(teamGrids), 'r.')
    %    title(sprintf('day: %d', j));
    %    pause(5);

    
    %    figure, plot(capabilities, todayGrid(teamGrids), 'b.')
    %    xlabel('total truck capacity')
    %    ylabel('debris remaining in grid');
        
    dayCount = dayCount - 1;
end

%figure, boxplot(outData.remaining, outData.day);