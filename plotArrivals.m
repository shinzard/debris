function plotArrivals(times, day)
%PLOTARRIVALS(TIMES, GROUP) - this function plots the count of arrivals 
% along time intervals given by the second input
%   TIMES - vector of arrival times
%   GROUP - vector of base time (difference will be plotted)

% TODOs:
% - handle truck/trailer combos...

close all;

times = sort(times);

days = unique(day);
rate = zeros(1,length(days));
open = zeros(1,length(days));

for i = 1:length(days)
    idx = find(day == days(i));
    figure(1);
    plot(times(idx), [1:length(idx)]);
    hold on;
    
    open(i) = times(idx(end))-times(idx(1));
    rate(i) = length(idx)/(open(i)*24);
    
    ideal = (times(idx)-times(idx(1)))*rate(i)*24;
    plot(times(idx), ideal, 'r--')
    
    figure(2);
    plot(times(idx), [1:length(idx)]-ideal', 'b.');
    hold on;
end

figure(1);
title('Arrival Process')
xlabel('Time')
ylabel('Arrival Count')

figure(2);
title('Arrival Process Error')
xlabel('Time')
ylabel('Arrival Count Difference from Ideal')

figure, plot(days, rate);
title('Arrival Rates');
xlabel('Time');
ylabel('Arrivals per hour');

figure, plot(days, open);
title('Duration Open');
xlabel('Time');
ylabel('Time open');


% Do chi-sq test of statistical fit to Poisson...

end

