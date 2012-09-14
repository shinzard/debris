function debrisanalysis(TRUCK_PLOTS, TDSR_PLOTS, SYS_PLOTS, RANGE)
% This script analyzes debris ticket data
%
% J.Brooks
% 13 May 2011

CLEAN = 1;
if nargin == 0
    TRUCK_PLOTS = 0;
    TDSR_PLOTS = 0;
    SYS_PLOTS = 0;
    RANGE = [];
end

% Load and Extract Data
if CLEAN
    %    clear all;
    close all;
    
    % Parameters
    % Set to true to pull clean copy
    %    TRUCK_PLOTS = 1;
    %    SYS_PLOTS = 1;  
    SCRUB_WHITES = 0; % removes WHITES for  haulpass and ASH, MUL for Storm
    SCRUB_ROE = 1; % For storm data only
    SCRUB_WEIGHT = 1; % For storm data only
    SCRUB_ZERO_HAUL = 1; % removes tickets with zero haul miles and appropriate project
    
%    dataType = 'HAULPASS';
%    data = csvread('/home/james/Documents/Research/Debris Task Assignment/USACE/Ike Data/DR-1792-LA_Cameron_Debris_Load_Records_mod.csv',1,0);
    dataType = 'STORM';
%   directory = 'M:/Documents/Research/Debris Task Assignment/USACE/AL Tornados/CompleteData/Cleaned/';
    directory = '/home/james/Documents/Research/Debris Task Assignment/USACE/AL Tornados/InitialData/';
    filename = dir(directory);
    
     project = [];
     tdsr = [];
     towerTime = [];
     QA = [];
     loadTime = [];
     QC = [];
     truckId = [];
     maxCapacity = [];
     loadPercent = [];
     loadVolume = [];
     ROE = [];
     lat = [];
     lon = [];
     haulMi = [];
     contents = [];
     subcont = [];
     subcont2 = [];
     file = [];
     invoice = [];
     subInvoice = [];
     
     totalTickets = 0;
     roeScrub = 0;
     whiteScrub = 0;
     haulScrub = 0;
     qcScrub = 0;
     rangeScrub = 0;
     weightScrub = 0;
    
    for i = 3:length(filename)
        try
            [data,txt] = xlsread([directory, '/', filename(i).name]);
        catch
            data = csvread([directory, '/', filename(i).name]);
        end
        %data = xlsread(filename(i),1,0);
        
        if strcmp(dataType,'HAULPASS')
            disp('HAULPASS FORMAT');
            truckId = data(:,1);
            contents = data(:,12);
            disp(sprintf('Initial data: %d', length(contents)));
            
            if SCRUB_WHITES
                good = find(truckId~=0 & contents~=4);
                disp(sprintf('Whites good: %d', length(good)));
            else
                good = find(truckId~=0);
                disp('Whites NOT SCRUBBED');
            end
            
            haulMi = data(good,6);
            if SCRUB_ZERO_HAUL
                good = good(find(haulMi > 0));
                disp(sprintf('Haul Mi scrubbed: %d', length(good)));
            end
            
            loadTime = data(good, 7);
            if ~isempty(RANGE)
                good = good(find(floor(loadTime)>=RANGE(1) & ...
                    floor(loadTime)>=RANGE(1)));
                disp(sprintf('Range size: %d', length(good)));
            end
            
            truckId = data(good,1);
            maxCapacity = data(good,2);
            loadPercent = data(good,3)/100;d
            loadVolume = data(good,4);
            site = data(good,5);
            tdsr = site;
            haulMi = data(good,6);
            loadTime = data(good,7);S
            towerTime = data(good,8);
            timeD1 = data(good,9);
            serverTime = data(good,10);
            timeD2 = data(good,11);
            contents = data(good,12);
            
        elseif strcmp(dataType, 'STORM')
            truckIdTmp = data(:,7);
            QCTmp = data(:,6);
            
            disp(sprintf('Initial data: %d', length(QCTmp)));
            good = find(truckIdTmp~=0 & QCTmp<123456789 & QCTmp~=0);       
            disp(sprintf('NonZero Truck, QC: %d', length(good)));
            qcScrub = qcScrub + length(good);
            clear truckIdTmp QCTmp;
            
            
            if SCRUB_ROE
                ROETmp = data(good, 11); 
                good = good(find(ROETmp==0));
                disp(sprintf('ROE Good: %d', length(good)));
                roeScrub = roeScrub + length(good);
                clear ROETmp;
            end
            
            if SCRUB_WEIGHT
                loadPercentTmp = data(good,9);        
                good = good(find(loadPercentTmp ~=0));
                disp(sprintf('Weight Good: %d', length(good)));
                weightScrub = weightScrub + length(good);
                clear loadPercentTmp;
            end

            if SCRUB_WHITES
                contentsTmp = data(good,15);
                good = good(find(contentsTmp < 3));
                disp(sprintf('Whites Good: %d', length(good)));
                whiteScrub = whiteScrub + length(good);
                clear contentsTmp
            end     
            
            if SCRUB_ZERO_HAUL
                haulMiT = data(good,6);
                projectT = data(good, 1);
                try
                    fileT = data(good, 18);
                catch
                    disp('no file data');
                end
                contentsT = data(good, 15);
                % Haul Mi non-zero, non-MATOC ROW only contract, no wet
                % debris
                
                try
                good = good(find(haulMiT > 0 & ...
                    fileT > 28 & ... % non-MATOC ROW only
                    floor(projectT) ~= 37 & ... % MATOC PPDR
                    projectT ~= 34 & ... % loose stumps
                    projectT ~= 35 & ... % lake martin (wet debris)
                    projectT ~= 36 & ... % P&J PPDR
                    projectT ~= 38 & ... % leaners/hangers
                    projectT ~= 39 & ... % stump extraction
                    ~strcmp(contentsT, 'Material') & ...  % No idea what these are (very few)
                    ~strcmp(contentsT, 'MUL') & ...  % Mulch loads
                    ~strcmp(contentsT, 'ASH')));  % Ash loads
                catch
                good = good(find(haulMiT > 0 & contentsT < 3));
                end
                    
                disp(sprintf('Data Good: %d', length(good)));
                haulScrub = haulScrub + length(good);
                clear haulMiT projectT fileT contentsT;
            end
            
            if ~isempty(RANGE)
                loadTimeT = data(good, 5);
                good = good(find(floor(loadTimeT)>=RANGE(1) & ...
                    floor(loadTimeT)<=RANGE(2)));
                disp(sprintf('Range size: %d', length(good)));
                rangeScrub = rangeScrub + length(good);
                clear loadTimeT;
            end
            
            
            project = [project; data(good, 1)];
            tdsr = [tdsr; data(good, 2)];
            towerTime = [towerTime; data(good,3)];
            QA = [QA; data(good, 4)];
            loadTime = [loadTime; data(good,5)];
            QC = [QC; data(good, 6)];
            truckId = [truckId; data(good,7)];
            maxCapacity = [maxCapacity; data(good,8)];
            loadPercent = [loadPercent; data(good,9)];
            loadVolume = [loadVolume; data(good,10)];
            ROE = [ROE; data(good, 11)];
            lat = [lat; data(good, 12)];
            lon = [lon; data(good, 13)];
            haulMi = [haulMi; data(good,14)];
            try % These are in complete dataset only
                contents = [contents; txt(good,15)];
                subcont = [subcont; txt(good, 16)];
                subcont2 = [subcont2; txt(good, 17)];
                file = [file; data(good, 18)];
                invoice = [invoice; data(good, 19)];
                subInvoice = [subInvoice; data(good, 20)];
            catch
                contents = [contents; data(good, 15)];
                subcont = [subcont; data(good, 16)];
                subcont2 = [subcont2; data(good, 17)];
            end
        end
        disp(sprintf('File: %s -- Percent Data Kept: %2.2f%%', ...
                 filename(i).name, length(good)/size(data,1)*100));
        totalTickets = totalTickets + size(data,1);
        clear data;
    end
    
    disp(sprintf('\n Total Tickets: %d -- Total Percent Data Kept: %2.2f%%', ...
        totalTickets, length(truckId)/totalTickets*100));
    disp(sprintf('Good QC: %d\nGood ROE: %d\nGood Weight: %d\n Good White: %d\nGood Haul: %d\nGood Range: %d\n', ...
        qcScrub, roeScrub, weightScrub, whiteScrub, haulScrub, rangeScrub));
    
end
close all;

% Allocate truck data matrices
trucks = unique(truckId)';
capacity = zeros(1,length(trucks));
numLoads = zeros(1,length(trucks));
effLoads = zeros(1,length(trucks));
avgLoadPercent = zeros(1,length(trucks));
dayService = zeros(1,length(trucks));
contentPercent = zeros(length(trucks),4);
sitePercent = zeros(length(trucks),4);
productivity = zeros(1,length(trucks));
loadsPerDay = zeros(1,length(trucks));
total = zeros(1,length(trucks));
startTime = zeros(1,length(trucks));
endTime = zeros(1,length(trucks));
dayData(length(trucks)).ahaul = [];
dayData(length(trucks)).effSpeed = [];
dayData(length(trucks)).numhauls = [];
dayData(length(trucks)).sameday = [];
dayData(length(trucks)).QCpercent = [];
dayData(length(trucks)).roundTrip = [];

corrected = 0;
notsorted = 0;

for i = [1:length(trucks)]
    loadIdx = find(truckId == trucks(i));

    if strcmp(dataType, 'HAULPASS')
       cnd = find(contents(loadIdx) == 1);
        veg = find(contents(loadIdx) == 2);
        mixed = find(contents(loadIdx) == 3);
        whites = find(contents(loadIdx) == 4);
        site1 = find(site(loadIdx) == 1);
        site2 = find(site(loadIdx) == 2);
        site3 = find(site(loadIdx) == 3);
        site4 = find(site(loadIdx) == 4);
    end
    capacity(i) = maxCapacity(loadIdx(1)); 
    numLoads(i) = length(loadIdx);
    effLoads(i) = sum(loadVolume(loadIdx))/capacity(i);
    avgLoadPercent(i) = mean(loadPercent(loadIdx));
    
    try
        endTime(i) = max([serverTime(loadIdx); towerTime(loadIdx); ...
                          loadTime(loadIdx)]);
        startTime(i) = min([serverTime(loadIdx); towerTime(loadIdx)]);
    catch                        
        endTime(i) = max([loadTime(loadIdx); towerTime(loadIdx)]);
        startTime(i) = min([towerTime(loadIdx); ...
                            loadTime(loadIdx)]);
    end
    
    % Count active days
    dayService(i) = 0;
    
    for j=floor(startTime(i)):floor(endTime(i))
        if ~isempty(find(floor(loadTime(loadIdx))==j))
            dayService(i) = dayService(i) + 1;
        end
    end
    
    % Initialize daily data structures
    dayData(i).day = zeros(1, dayService(i)).*NaN;
    dayData(i).ahaul = zeros(1,dayService(i)).*NaN;
    dayData(i).effSpeed = zeros(1,dayService(i)).*NaN;
    dayData(i).numhauls = zeros(1,dayService(i)).*NaN;
    dayData(i).sameday = zeros(1,dayService(i)).*NaN;
    dayData(i).QCpercent = zeros(1,dayService(i)).*NaN;
    dayData(i).QCmaj = zeros(1,dayService(i)).*NaN;
    dayData(i).roundTrip = zeros(1,dayService(i)).*NaN;
    dayData(i).perfL = zeros(1,dayService(i)).*NaN;
    dayData(i).perfW = zeros(1,dayService(i)).*NaN;
    
    
    % Find daily data for each truck
    if exist('QC')
        QCs = unique(QC);
    end
    
    k = 1;
    for j=floor(startTime(i)):floor(endTime(i))
        idxL = find(floor(loadTime(loadIdx))==j);
        idxT = find(floor(towerTime(loadIdx))==j);
        
        % Loads picked up and delivered same day
        idxC = loadIdx(intersect(idxL, idxT));
        idxAM = find(loadTime(idxC)-j < 0.3);
        
        if ~isempty(idxAM)
            loadTime(idxC(idxAM)) = loadTime(idxC(idxAM)) + 0.5;
            corrected = corrected + length(idxAM);
        end
        
        [tmp, tmpI] = sort(loadTime(idxC));
        idxC = idxC(tmpI);
        
        if ~isempty(idxC)
            dayData(i).day(k) = j;
            dayData(i).ahaul(k) = mean(haulMi(idxC));
            dayData(i).effSpeed(k) = mean(haulMi(idxC)./(towerTime(idxC) - ...
                                                         loadTime(idxC))/24);
            dayData(i).numhauls(k) = length(haulMi(idxC)>0);
            dayData(i).sameday(k) = length(idxC)/length(idxL);
            
            if exist('QC')
                tmp = disagreement([truckId(idxC), QC(idxC)]);
                dayData(i).QCpercent(k) = tmp.MajorityPerc;
                dayQCs = unique(QC(idxC));
                dayData(i).QCmaj(k) = dayQCs(tmp.MajorityID);
            end
            dayData(i).roundTrip(k) = ...
                mean(diff(sort(loadTime(idxC))))*24;
            
            % Check proper time order
            tmp = [];
            for l = 1:length(idxC)
                tmp = [tmp, loadTime(idxC(l)), towerTime(idxC(l))];
            end
            
            if ~issorted(tmp)
                notsorted = notsorted + 1;
                dayData(i).perfL(k) = 0; 
                dayData(i).perfW(k) = 0;                
            else
                dayData(i).perfL(k) = sum(loadPercent(idxC).*haulMi(idxC)./ ...
                                          ((towerTime(idxC) - ...
                                            loadTime(idxC))*24));
                dayData(i).perfW(k) = sum(haulMi(idxC(1:end-1))./((loadTime(idxC(2:end))- ... 
                                                                  towerTime(idxC(1:end-1)))*24));
            end
            
            if dayData(i).perfL(k) < 0 | dayData(i).perfW(k) < 0
                disp('LESS THAN ZERO')
                trucks(i)
                j
                loadPercent(idxC)
                haulMi(idxC)
                loadTime(idxC)- j
                towerTime(idxC) - j
            end
            
            
            % Increment day
            k = k + 1;
        end
    end
    
    
    loadsPerDay(i) = numLoads(i)/dayService(i);
    total(i) = sum(loadVolume(loadIdx));
    productivity(i) = total(i)/capacity(i)/dayService(i);
    if strcmp(dataType, 'HAULPASS')
        contentPercent(i,1) = length(cnd)/numLoads(i)*100;
    contentPercent(i,2) = length(veg)/numLoads(i)*100;
    contentPercent(i,3) = length(mixed)/numLoads(i)*100;
    contentPercent(i,4) = length(whites)/numLoads(i)*100;
        sitePercent(i,1) = length(site1)/numLoads(i)*100;
        sitePercent(i,2) = length(site2)/numLoads(i)*100;
        sitePercent(i,3) = length(site3)/numLoads(i)*100;
        sitePercent(i,4) = length(site4)/numLoads(i)*100;
    end
end

disp(sprintf('Number of time stamps corrected: %d', corrected));
disp(sprintf('Number out of order: %d', notsorted));

if TRUCK_PLOTS
    figure, plot([1:length(trucks)], numLoads, 'b.');
    text(10, 10, sprintf('Average: %2.2f\nStd. Dev.: %2.2f',...
                         mean(numLoads), std(numLoads)));
    title('Number of Loads'), xlabel('Truck ID'), ylabel('Count');
    
    figure, plot([1:length(trucks)], effLoads, 'b.');
    text(10, 10, sprintf('Average: %2.2f\nStd. Dev.: %2.2f',...
                         mean(effLoads), std(effLoads)));
    title('Effective Loads'), xlabel('Truck ID'), ylabel('Count');
    
    figure, plot([1:length(trucks)], avgLoadPercent, 'b.');
    text(10, .5, sprintf('Average: %2.2f\nStd. Dev.: %2.2f',...
                         mean(avgLoadPercent), std(avgLoadPercent)));
    title('Average load percent'), xlabel('Truck ID'), ylabel('Avg. Percentage');
    
    figure, plot([1:length(trucks)], dayService, 'b.');
    text(10, 5, sprintf('Average: %2.2f\nStd. Dev.: %2.2f',...
                        mean(dayService), std(dayService)));
    title('Time (days) in Service'), xlabel('Truck ID'), ylabel('Number of Days');
    
    figure, bar([1:length(trucks)], contentPercent, 'stacked'), ...
        axis([0,length(trucks),0,100]);
    if  strcmp(dataType, 'HAULPASS')
        legend({'C & D', 'Veg', 'Mixed', 'Whites'},'Location', ...
               'SouthWest');
    else
        legend({'C & D', 'Veg', 'Ash', 'Mul'},'Location', ...
               'SouthWest');
    end
    title('Proportion of loads by contents')
    xlabel('Truck ID'), ylabel('Percentage');
    
    figure, bar([1:length(trucks)], sitePercent, 'stacked')
    axis([0,length(trucks),0,100]);
    legend({'McManus', 'Romero', 'Wilkerson', 'Empty'},'Location', 'SouthWest');
    title('Proportion of site'), xlabel('Truck ID'), ylabel('Percentage');
    
    figure, plot([1:length(trucks)], productivity,'b.');
    text(10, 5, sprintf('Average: %2.2f\nStd. Dev.: %2.2f',...
                        mean(productivity), std(productivity)));
    title('Truck Productivity (eff. loads/day of service)')
    xlabel('Truck ID'), ylabel('Full Loads/Day');
    
    figure, plot([1:length(trucks)], loadsPerDay, 'b.')
    text(10, 5, sprintf('Average: %2.2f\nStd. Dev.: %2.2f',...
                        mean(loadsPerDay), std(loadsPerDay)));
    title('Truck Productivity (loads/day of service)')
    xlabel('Truck ID'), ylabel('Actual Loads/Day');
end

% Allocate System Metrics
duration = unique(floor(loadTime)); % many zero entries due to one
                                    % early ticket
%[floor(min(startTime)):floor(max(endTime))];
numTrucksL = zeros(1,length(duration));
numTrucksT = zeros(1,length(duration));
numTicketsL = zeros(1,length(duration));
numTicketsT = zeros(1,length(duration));
numTdsrs = zeros(1,length(duration));
cyL = zeros(1,length(duration));
cyT = zeros(1,length(duration));
pfL = zeros(1,length(duration));
pfT = zeros(1,length(duration));
attd = zeros(1,length(duration));
ahaul = zeros(1, length(duration));
numQCs = zeros(1, length(duration));

for j = 1:length(duration)
    idxL = find(floor(loadTime)==duration(j));
    idxT = find(floor(towerTime)==duration(j));
    
    % Loads picked up and delivered same day
    idxC = intersect(idxL, idxT);
    
    % Number of QCs active (loading sites??)
    if exist('QC')
        numQCs(j) = length(unique(QC(idxL)));
    end
    
    % Number of tickets issued
    numTicketsL(j) = length(idxL);
    numTicketsT(j) = length(idxT);
    
    % Number of trucks active
    numTrucksL(j) = length(unique(truckId(idxL)));
    numTrucksT(j) = length(unique(truckId(idxT)));
    
    % Number of TDSRs active
    numTdsrs(j) = length(unique(tdsr(idxT)));
    
    % Total cubic yards
    cyL(j) = sum(loadVolume(idxL));
    cyT(j) = sum(loadVolume(idxT));
    
    % Percent full
    pfL(j) = mean(loadPercent(idxL));
    pfT(j) = mean(loadPercent(idxT));
    
    % Average time to deliver
    attd(j) = mean(towerTime(idxC) - loadTime(idxC));
    
    % Average haulage
    ahaul(j) = mean(haulMi(idxC));
end

if SYS_PLOTS

    figure, plot(duration - duration(1), numTrucksL, 'b')
    title(['Number of trucks active']), xlabel('Time(days)')  
    hold on; plot(duration - duration(1), numTrucksT, 'r')
    legend({'loading', 'tower'});
    
    figure, plot(duration - duration(1), numQCs, 'b')
    title(['Number of QCs active']), xlabel('Time(days)')  
    
    figure, plot(duration - duration(1), numTdsrs, 'b')
    title(['Number of TDSRs active']), xlabel('Time(days)')  
    
    figure, plot(duration - duration(1), numTicketsL, 'b')
    title(['Number of tickets']), xlabel('Time(days)')  
    hold on; plot(duration - duration(1), numTicketsT, 'r')
    legend({'loaded', 'tower'});
    
    figure, plot(duration - duration(1), numTicketsL./numTrucksL, ...
                 'b')
    hold on, plot(duration - duration(1), numTicketsT./numTrucksT, ...
                 'r')
    legend({'loading', 'tower'});
    title('Productivity: avg tickets/truck'), xlabel('Time(days)');
    
    figure, plot(duration-duration(1), cyL, 'b')
    hold on, plot(duration-duration(1), cyT, 'r')
    legend({'loaded', 'tower'});
    title('Cubic Yards of Debris'), xlabel('Time(days)')  
    
    figure, plot(duration-duration(1), cyL./numTrucksL, 'b')
    hold on, plot(duration-duration(1), cyT./numTrucksT, 'r')
    legend({'loaded', 'tower'});
    title('Cubic Yards of Debris Per Truck'), xlabel('Time(days)')  
    
    figure, plot(duration-duration(1), pfL, 'b')
    hold on, plot(duration-duration(1), pfT, 'r')
    legend({'loaded', 'tower'});
    title('Average Percent Full'), xlabel('Time(days)')  
    
    figure, plot(duration-duration(1), attd*24, 'b')
    title('Average Time to deliver'), xlabel('Time(days)')  
    ylabel('Time (hrs)');
    
    figure, plot(duration-duration(1), ahaul, 'b');
    title('Average Haul Distance'), xlabel('Time(days)')  
    ylabel('Haul Distance (mi)')
    
    figure, plot(duration-duration(1), ahaul./(attd*24), 'b.')
    title('Effective delivery speed')
    ylabel('Speed (mph)');
    xlabel('Time(days)') 

    % figure, bar(trucks, [startTime;endTime], 'stacked', 'wb');

end

% Get TDSR statistics
if strcmp(dataType, 'HAULPASS')
    % Site metrics
    sites = [1:3]
    edges = [0:1/6:3, 30];
    color = ['bgr'];
    day1 = min(towerTime)
    totals = zeros(1,length(sites));
    days = zeros(1,length(sites));
    for i = sites
        tickets = find(site == i);
        [tmp, idx] = sort(towerTime(tickets));
        [n, bin] = histc(diff(sort(towerTime(tickets)*24)),edges); 
        
        if SYS_PLOTS
            figure, bar([edges(1:end-1), 3+1/6],n,'histc');
            title('Interarrival ticket time'), xlabel('Time(hrs)');
            text(2, 3/4*max(n), sprintf('Total Tickets: %d', sum(n)));
            figure(20), plot(tmp - day1, cumsum(loadVolume(idx)), ...
                             color(i)), hold on;
        end
        
        days(i) = tmp(end) - tmp(1);
        totals(i) = sum(loadVolume(tickets));
        
        nonZero = find(haulMi(tickets) ~= 0);
        if SYS_PLOTS
            figure, hist(haulMi(tickets(nonZero))), title('HaulMiles');
            text(15, 50, sprintf('Mean: %2.2f', mean(haulMi(tickets(nonZero)))));
            text(15, 30, sprintf('Median: %2.2f', ...
                                 median(haulMi(tickets(nonZero)))));
        end
    end

    if SYS_PLOTS
        figure(20), title('Total Debris Collected at Each TDSR');
        xlabel('Time (day)'), ylabel('Cumulative cubic yards')
        legend({'McManus', 'Romero', 'Wilkerson'});
        text(10, 10000, sprintf('Total Collected: %2.0f cy', sum(loadVolume)));
        text(10, 8000, sprintf('Productivity:\n    McManus : %2.0f cy/day\n    Romero: %2.0f cy/day\n    Wilkerson: %2.0f cy/day',...
                               totals(1)/days(1), totals(2)/days(2), ...
                               totals(3)/days(3)));
    end
else 

    % Site metrics
    tdsrs = unique(tdsr);
    edges = [0:0.5/60:3];          % in hours
    color = ['bgrmkcy'];
    day1 = min(towerTime);
    totals = zeros(1,length(tdsrs));
    days = zeros(1,length(tdsrs));
    p = zeros(1,length(tdsrs));
    l = zeros(1,length(tdsrs));
    r = zeros(1,length(tdsrs));
    disp(sprintf('TDSR\tNum.Tickets\tNum.Trucks\tTotalDebris\tAvg.InterTime\tDays\tAvg.Haul\tPerc<5min\tRatio\tUtil.\tArrRate\tMu'))
    for i = 1:length(tdsrs)
        tickets = find(tdsr == tdsrs(i));
        interTimes = [];
        utilization = NaN*zeros(1,length(duration));
        arrivalRate = NaN*zeros(1,length(duration));
        meanService = NaN*zeros(1,length(duration));
        ratio = NaN*zeros(1,length(duration));

        % Calculate TDSR queue properties
        for j = 1:length(duration)
            today = find(floor(towerTime(tickets)) == duration(j));

            if length(today) > 10
                [tmp, idx] = sort(towerTime(tickets(today)));
                
                % Check for truck/trailer combo
                tmp2 = diff(tmp)*24;
                % combo = find(tmp2 < 1/60 & ...
                %          diff(loadTime(tickets(today(idx)))) < ...
                %           1/60 );
                notcombo = find(tmp2 > 1/60 | ...
                            diff(loadTime(tickets(today(idx)))) > ...
                            1/60 );
                                
                interToday = tmp2(notcombo); % hours
                
                interTimes = [interTimes; interToday];
                
                lines = interToday < 5/60; % hours
                linesI = find(interToday < 5/60);
                noLine = length(find(lines==0));
                
                ratio(j) = (length(lines)-noLine)/length(lines);

                meanService(j) = mean(interToday(linesI));
                
                busyTime = ( (1 + noLine)*meanService(j)) + ...
                    lines'*interToday;
                
                totalTime = (tmp(end) - tmp(1))*24;
                
                if totalTime > 2
                    arrivalRate(j) = length(notcombo)/totalTime;
                end

                utilization(j) = busyTime/totalTime;
                days(i) = days(i) + 1;
            end
        end

        if TDSR_PLOTS
            figure(200)
            plot(duration - duration(1), arrivalRate, ... 
                 color(mod(i,length(color)-1)+1) );
            hold on;
        end
        
        % Calculate the average utilization for this tdsr
        r(i) = nanmean(ratio);
        p(i) = nanmean(utilization);
        l(i) = nanmean(arrivalRate);
        u(i) = nanmean(meanService);
        
        if length(interTimes)>1
            [n, bin] = histc(interTimes,[-inf, edges, inf]); 
            if TDSR_PLOTS
                figure, bar([-1, edges, 4], n,'histc');
                title(sprintf('Interarrival ticket time: %d', tdsrs(i)))
                xlabel('Time(hrs)');
                text(2, 3/4*max(n), sprintf('Total Tickets: %d', ...
                                            sum(n)));
                % figure(20), plot(tmp - day1, cumsum(loadVolume(idx)), ... 
                %                  color(mod(i, length(color)-1)+1)), hold on;
            end
            
            totals(i) = sum(loadVolume(tickets));
            
            nonZero = find(haulMi(tickets) ~= 0 & haulMi(tickets) ...
                           < 56);
            
            if TDSR_PLOTS
                figure, hist(haulMi(tickets(nonZero)));
                title(sprintf('HaulMiles: %d', tdsrs(i)));
                text(15, 50, sprintf('Mean: %2.2f', mean(haulMi(tickets(nonZero)))));
                text(15, 30, sprintf('Median: %2.2f', ...
                                     median(haulMi(tickets(nonZero)))));
            end
            % Display table 
            disp(sprintf('%d\t%d\t%d\t%5.2f\t%2.2f\t%d\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f', ...
                         tdsrs(i), ...
                         length(tickets), length(unique(truckId(tickets))), ...
                         totals(i), mean(interTimes), days(i), ...
                         mean(haulMi(tickets(nonZero))), n(2)/sum(n), ...
                         r(i), p(i), l(i), u(i)));
        end
    end
    if TDSR_PLOTS
        figure(20), title('Total Debris Collected at Each TDSR');
        xlabel('Time (day)'), ylabel('Cumulative cubic yards')
        %legend(num2cell(tdsrs));
        text(10, 10000, sprintf('Total Collected: %2.0f cy', ...
                                sum(loadVolume)));
    end
end

%qcSpan = generateQcData(QC, lat, lon, loadTime, truckId, trucks);

% Save all data to base workspace...
tmp = whos;
for i = 1:length(tmp)
    assignin('base', tmp(i).name, eval(tmp(i).name))
end