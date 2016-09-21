% Load data, run all scripts needed to generate csv file for
% statistical analysis (R)
%
% 28 May 2012
% J.Brooks

clear all;

load Data/DEBRIS_AL_COMPLETE
%load DEBRIS_ANALSIS_AL_COMPLETE         % be sure to exclude
                                        % file < 18
%%load DEBRIS_ANALSIS_AL_PEAK_19feb %for verification
%load DEBRIS_ANALYSIS_AL % for verification

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

% Overwrite to focus on main event (many zero days, very low activity..)
duration = [40672:40817]; % Complete data set
%duration = [40695:40720]; % original June data (exclude 40721 which
                          % has 10 trucks only?)

[teams, sub] = selectData(tickets, duration);

[performance,teamSize] = teamPerfMeasures(tickets,teams);

% note here, duration is the extent over which to calculate
% familiarity (can be different than duration sent to selectData);
% will only record the data for the selected teams
[familiarity,teamHaul,teamCapacity] = familiarity(tickets,teams,duration);

% Write Data
filename = ['CompleteData_NewMeasures_',num2str(now), '.csv']
fid = fopen(filename, 'w+');            % create and overwrite data
fprintf(fid, '%s\n', ['Team,Effectiveness,Efficiency,Inequality,Size,' ...
              'Familiarity,Haul,MaxCap,Sub']);
fprintf(fid, '%d, %6.2f, %6.2f, %6.2f, %d, %6.2f, %6.2f, %6.2f, %d\n', [teams'; ...
        performance.eff; performance.ecy; performance.eq; teamSize; ...
        familiarity; teamHaul; teamCapacity; sub]);
fclose(fid);

% ------------------------------
% OLD SCRIPT
% ------------------------------
%eqData = equality(lat,lon,duration,loadTime,QC,truckId,haulMi, ...
%                   capacity,loadVolume,trucks);

% New equity measure: workload equity
%eqData2 = equalityWorkload(QC,loadTime,truckId,trucks,eqData);

%fData = fluidity3(dayData,trucks,loadTime,truckId,subcont, ...
%                  loadVolume,eqData,QC,capacity, dayService, ...
%                  0,['CompleteData_',num2str(now), '.csv'],duration,0, ...
%                  towerTime,lat,lon,loadPercent, haulMi);