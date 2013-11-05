close all;
clear all;

load fireIslandSimData_20Feb2013

%boxplot([unifA(:,2);randA(:,2)+1], [unifA(:,1);randA(:,1)+1], 'colors', hsv(7));

boxplot([unifA(:,2)], [unifA(:,1)], 'colors', 'b', 'labels', {'','','','','','','',''}, 'sym', 'b+');
hold on;
axis manual;
boxplot([randA(:,2)], [randA(:,1)], 'colors', 'r');
axis([0 9 6 20]);
xlabel('Number of Trucks');
ylabel('Mean Throughput Rate');
title('System Throughput');

figure, boxplot([unifA(:,2)], [unifA(:,1)], 'colors', 'b', 'labels', {'','','','','','','',''});
hold on;
axis manual;
boxplot([unifB(:,2)], [unifB(:,1)], 'colors', 'r');
axis([0 9 6 20]);
xlabel('Number of Trucks');
ylabel('Mean Throughput Rate');
title('System Throughput');

figure, boxplot([unifAvaried(:,2)], [unifAvaried(:,1)], 'colors', 'b', 'labels', {'','','','','','','',''});
hold on;
axis manual;
boxplot([randAvaried(:,2)], [randAvaried(:,1)], 'colors', 'r');
axis([0 9 6 20]);
xlabel('Number of Trucks');
ylabel('Mean Throughput Rate');
title('System Throughput');

unifAmeans = grpstats(unifA(:,2), unifA(:,1));
randAmeans = grpstats(randA(:,2), randA(:,1));
figure,
plot([6:6:48], unifAmeans, 'b');
hold on;
plot([6:6:48], randAmeans, 'r');

%gscatter(unifA(:,2), unifA(:,1), [repmat('A',1,40), repmat('B',1,40)]);