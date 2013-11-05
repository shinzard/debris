% Load data, run all scripts needed to generate csv file for
% statistical analysis (R)
%
% 28 May 2012
% J.Brooks

clear all;

%load DEBRIS_AL_COMPLETE                 % LAST 
%load DEBRIS_ANALSIS_AL_COMPLETE         % be sure to exclude
                                        % file < 18
%%load DEBRIS_ANALSIS_AL_PEAK_19feb %for verification
load DEBRIS_ANALYSIS_AL % for verification (PILOT DATA)

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
clear RANGE i idxT p

clear equality % needed only for peak data (used for verification)


% make structure of remaining ticket data
tickets = [];
x = whos;
for i = 1:length(x)
    if strcmp(x(i).name, 'tickets')
        continue;
    end
    tickets = setfield(tickets, x(i).name, eval(x(i).name));
end

retain tickets;                         % clear everything else

[perf] = truckPerfMeasures(tickets);

% Write Data
filename = ['truckPerf_',num2str(now), '.csv']
fid = fopen(filename, 'w+');            % create and overwrite data
fprintf(fid, '%s\n', 'TruckID,QC,Day,Perf,NumLoads,SubCont');
fprintf(fid, '%d, %d, %d, %6.2f, %d, %d\n', [perf.id; perf.qc;, ...
                    perf.d; perf.eff; perf.nL; perf.sc]);
fclose(fid);

% --------------------
% MAKE PLOTS
% --------------------
d = floor(tickets.loadTime);
days = unique(d);

teams = zeros(1,length(days));
nTrucks = zeros(1,length(days));
fluid = zeros(1,length(days));

for i = 13:length(days)
    idx = find(d == days(i));
    teams(i) = length(unique(tickets.QC(idx)));
    thisTrucks = unique(tickets.truckId(idx));
    nTrucks(i) = length(thisTrucks);
    notsingleTeam = zeros(1,nTrucks(i));
    for j = 1:nTrucks(i)
        tmp = idx(find(tickets.truckId(idx)==thisTrucks(j)));
        notsingleTeam(j) = (length(unique(tickets.QC(tmp))) > 1);
    end
    fluid(i) = sum(notsingleTeam)/nTrucks(i);
end

t = days(13:end-1) - min(days(13:end-1)) + 1;
figure, plot(t, fluid(13:end-1))
xlabel('Days')
ylabel('Percentage')
title('Fluidity: Percentage of Trucks Switching DRT')

figure, plotyy(t, teams(13:end-1), t, nTrucks(13:end-1)./teams(13:end-1))
xlabel('Days')
ylabel('Number DRTs/Average DRT Size')
title('Number DRTs and Average Team Size');
legend({'DRTs', 'Size'})