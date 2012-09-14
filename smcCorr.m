% Script for SMC Correspondence paper

% be sure to comment out AL and uncomment Ike
debrisanalysis(1,1,1)

for i = 1:length(dayData)
    meanPerf(i) = mean(dayData(i).perfL./dayData(i).numhauls);
end

figure, plot(meanPerf)

figure, plot(dayService, meanPerf, 'b.')

% use data fit for linear regression:
%y = 0.764*z + 14.2
%where z = (x - 5.63)/4.53


% there are NaNs that need to be looked into

% look at temporal effect per individual...
% could also plot one point for each day/individual