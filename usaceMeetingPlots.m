% load DEBRIS_ANALYSIS_AL % contains duplicates...
load DEBRIS_AL_COMPLETE

figure, boxplot(loadPercent, floor(loadTime)-min(floor(loadTime)), 'plotstyle', 'compact', 'labelorientation', 'horizontal')
title('Boxplot of Load Percent Data by Day')
xlabel('Day from First Ticket')
ylabel('Percent Full')

meanpf = zeros(1, length(duration));

for i = 1:length(duration)
    idx = find(floor(loadTime) == duration(i));
    meanpf(i) = nanmean(loadPercent(idx));
end

hold on;

plot(duration-duration(1), meanpf(:,1))

medpf = zeros(1, length(duration));
for i = 1:length(duration)
    idx = find(floor(loadTime) == duration(i));
    medpf(i) = nanmedian(loadPercent(idx));
end
plot(duration-duration(1), medpf)

subs = unique(subcont);
subcont = match(subs,subcont);

struc = groupId(truckId, QC, loadTime, subcont, haulMi, ...
0.9, 1);

figure, hist(struc.crewSize, [1:max(struc.crewSize)])
xlabel('Crew Size')
ylabel('Count')
title('Histogram of Daily Crew Size')

%set(gca, 'TickDir', 'in')


subSize = zeros(1,29); 
for i = 1:29
    idx = find(subcont == subs(i));
    subSize(i) = length(unique(truckId(idx)));
end

figure, hist(subSize, [5:10:max(subSize)])
axis([0 130 0 6.5])
title('Histogram of Subcontractor Size')
xlabel('Size (number of trucks)')
ylabel('Count')

save GROUP_ID_OUT struc


QCDay = QC*1e6 + floor(loadTime);
for i = 1:length(trucks)
    idx = find(truckId == trucks(i));
    crewsPerDay(i) = length(unique(QCDay(idx)))/length(unique(floor(loadTime(idx))));
end

figure, hist(crewsPerDay, [1.05:0.1:2.2])
title('Histogram of Number of Crews Per Day of Service')
xlabel('Number of Crews (avg)');
ylabel('Count (trucks)')

for i = 1:length(tdsrs)
    idx = find(tdsr == tdsrs(i));
    tdsrTickets(i) = length(idx);
    tdsrTotals(i) = sum(loadVolume(idx));
    tdsrDays(i) = length(unique(floor(towerTime(idx))));
end

tdsrCounts = [];
tdsrId = [];

for i = duration(13:end)
    idx = find(floor(towerTime) == i);
    for j = 1:length(tdsrs)
        idx2 = find(tdsr(idx) == tdsrs(j));
        tdsrCounts = [tdsrCounts, length(idx2)];
        tdsrId = [tdsrId, tdsrs(j)];
    end
end

figure, boxplot(tdsrCounts, tdsrId, 'plotstyle', 'compact')
title('Boxplot of Tickets Per Day at Each TDSR')
xlabel('TDSR ID')
ylabel('Number of Tickets')

load OWNER_DATA
owners = unique(ownerData(:,2));

for i = 1:length(owners)-1 % -1 avoids 9999 blank indicator
    oSize(i) = length(find(ownerData(:,2)==owners(i)));
end

figure, hist(oSize, [1:max(oSize)])
axis([0 max(oSize) 0 450])
title('Histogram of Owner Size')
xlabel('Number of Trucks');
ylabel('Count');
