clear all;
close all;

%load 85Prices;                          % only T=0.73, need to
                                        % hardcode best below
                                        % (uncomment)

%load 85Prices_all_2;                    % contains both T = 0.73,
                                        % 0.2, optimal prices, and
                                        % optimal partition (_all
                                        % didn't have optimal
                                        % partition, prices)...both
                                        % have wrong flow, length
                                        % inequity measures....

%load 85Prices_final;                    % fixed inequity (missing
%T=0.2 with optimal prices)
%measures...
%load 85Prices_final_2;                    % fixed inequity
%measures (but length inequity measure was messed up for optimal prices...)

%load 85Prices_final_3;                  % has wrong flow inequity measure...

% --------------------------------------------------
%load 85Prices_final_4;                  % THIS WAS USED IN
                                        % DISSERTaTION FOR STUDY ONE....
% --------------------------------------------------

%load dynamic_2;                         % 20 days...

%load dynamic_3_nofamScaling;  % tried not scaling familiarity (20 days)


%load dynamic_60d_w;                       % 60 days, week price
                                          % updates. (still
                                          % constant travel times

%load dynamic_60d_w_dijUpdate;             % same as above but with
                                        % randomly varying
                                        % distances (may not have
                                        % full duration for optimal
                                        % prices; also missing optimal...)

%load dynamic_60d_w_dijUpdate_full;  % above corrected
load dynamic_60d_w_dijUpdate_01min      % 0.01 minimum dij

% correct optimal price strategy:
% priceStrategy(end-49:end) = numPrices + 2;

% Correct temperature...
%T(end-49:end) = 0.2;

h = zeros(2,numPrices);
for i = 1:numPrices
    h(1,i) = mean(throughput(priceStrategy == i & T == 0.73));
    h(2,i) = 1.96*std(throughput(priceStrategy == i & T == 0.73))/sqrt(reps);
    h_lt(1,i) = mean(throughput(priceStrategy == i & T == 0.2));
    h_lt(2,i) = 1.96*std(throughput(priceStrategy == i & T == 0.2))/sqrt(reps);

    r(1,i) = mean(revenue(priceStrategy == i & T == 0.73));
    r(2,i) = 1.96*std(revenue(priceStrategy == i & T == 0.73))/sqrt(reps);
    r_lt(1,i) = mean(revenue(priceStrategy == i & T == 0.2));
    r_lt(2,i) = 1.96*std(revenue(priceStrategy == i & T == 0.2))/sqrt(reps);

    f(1,i) = mean(meanFam(priceStrategy == i & T == 0.73));
    f(2,i) = 1.96*std(meanFam(priceStrategy == i & T == 0.73))/sqrt(reps);
    f_lt(1,i) = mean(meanFam(priceStrategy == i & T == 0.2));
    f_lt(2,i) = 1.96*std(meanFam(priceStrategy == i & T == 0.2))/sqrt(reps);

    s(1,i) = mean(maxSize(priceStrategy == i & T == 0.73));
    s(2,i) = 1.96*std(maxSize(priceStrategy == i & T == 0.73))/sqrt(reps);
    s_lt(1,i) = mean(maxSize(priceStrategy == i & T == 0.2));
    s_lt(2,i) = 1.96*std(maxSize(priceStrategy == i & T == 0.2))/sqrt(reps);

    m(1,i) = mean(meanMu(priceStrategy == i & T == 0.73));
    m(2,i) = 1.96*std(meanMu(priceStrategy == i & T == 0.73))/sqrt(reps);
    m_lt(1,i) = mean(meanMu(priceStrategy == i & T == 0.2));
    m_lt(2,i) = 1.96*std(meanMu(priceStrategy == i & T == 0.2))/sqrt(reps);

    I(1,i) = mean(inequity(priceStrategy == i & T == 0.73));
    I(2,i) = 1.96*std(inequity(priceStrategy == i & T == 0.73))/sqrt(reps);
    I_lt(1,i) = mean(inequity(priceStrategy == i & T == 0.2));
    I_lt(2,i) = 1.96*std(inequity(priceStrategy == i & T == 0.2))/sqrt(reps);

    I2(1,i) = mean(inequity2(priceStrategy == i & T == 0.73));
    I2(2,i) = 1.96*std(inequity2(priceStrategy == i & T == 0.73))/sqrt(reps);
    I2_lt(1,i) = mean(inequity2(priceStrategy == i & T == 0.2));
    I2_lt(2,i) = 1.96*std(inequity2(priceStrategy == i & T == 0.2))/sqrt(reps);

    I3(1,i) = mean(inequity3(priceStrategy == i & T == 0.73));
    I3(2,i) = 1.96*std(inequity3(priceStrategy == i & T == 0.73))/sqrt(reps);
    I3_lt(1,i) = mean(inequity3(priceStrategy == i & T == 0.2));
    I3_lt(2,i) = 1.96*std(inequity3(priceStrategy == i & T == 0.2))/sqrt(reps);

    I4(1,i) = mean(inequity4(priceStrategy == i & T == 0.73));
    I4(2,i) = 1.96*std(inequity4(priceStrategy == i & T == 0.73))/sqrt(reps);
    I4_lt(1,i) = mean(inequity4(priceStrategy == i & T == 0.2));
    I4_lt(2,i) = 1.96*std(inequity4(priceStrategy == i & T == 0.2))/sqrt(reps);
end

diffFromOpt = zeros(1,numPrices);
optPrice = [1, 0.9345, 0.9515, 1.0875, 0.9438, 1.1404]; %
                                                        % pre-calculated...
optPrice = optPrice./sum(optPrice);

for i = 1:numPrices
    diffFromOpt(i) = norm(priceOpt(i,:) - optPrice);
end

% Partition optimal
opt = mean(throughput(priceStrategy == numPrices+2)); % T
                                                      % irrelevant
                                                      % for this
                                                      % case
opt_ci = 1.96*std(throughput(priceStrategy == numPrices+2))/sqrt(reps);
opt_f = mean(meanFam(priceStrategy == numPrices+2)); % T
opt_s = mean(maxSize(priceStrategy == numPrices+2)); % T
opt_m = mean(meanMu(priceStrategy == numPrices+2)); % T

%[best, bestIdx] = max(h(1,:));
%best = 16.5407;                         % from separate simulation
                                        % runs using the optimal
                                        % prices...
                                        %best_ci = 0.1462;

% Optimal prices
best = mean(throughput(priceStrategy == numPrices+1 & T == 0.73));
best_ci = 1.96*std(throughput(priceStrategy == numPrices+1 & T == 0.73))/sqrt(reps);
best_f = mean(meanFam(priceStrategy == numPrices+2 & T == 0.73)); % T
best_s = mean(maxSize(priceStrategy == numPrices+2 & T == 0.73)); % T
best_m = mean(meanMu(priceStrategy == numPrices+2 & T == 0.73)); % T

best_lt = mean(throughput(priceStrategy == numPrices+1 & T == 0.2));
best_f_lt = mean(meanFam(priceStrategy == numPrices+2 & T == 0.2)); % T
best_s_lt = mean(maxSize(priceStrategy == numPrices+2 & T == 0.2)); % T
best_m_lt = mean(meanMu(priceStrategy == numPrices+2 & T == 0.2)); % T

%[sortedDiff, idx] = sort(diffFromOpt);
idx = [1 2 3];
% figure, errorbar([0, sortedDiff], [opt/best, h(1,idx)/best], [opt_ci/best, h(2,idx)/best], 'b', 'LineWidth', ...
%                  0.75);

figure, errorbar(h(1,idx)/best, h(2,idx)/best, 'b.');
hold on;
plot([0 numPrices+2], [opt/best, opt/best], 'b--');
plot([0 numPrices+2], [1, 1], 'b-.');
plot([0 numPrices+2], [1+best_ci/best, 1+best_ci/best], 'b:');
plot([0 numPrices+2], [1-best_ci/best, 1-best_ci/best], 'b:');
title(['Relative Throughput, Sorted by Distance from Optimal ' ...
       'Price']);
%xlabel('Price Strategy (sorted by distance from optimal)')
xlabel('Price Strategy');
ylabel('Relative Throughput');
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Uniform', 'Time', 'Distance'});             % remove numeric labels
axis([0 87 0.65 1.1]);

% --------------------------------------------------
% Revenue
% --------------------------------------------------
optR = mean(revenue(priceStrategy == numPrices+2)); % T irrelevant here
optR_ci = 1.96*std(revenue(priceStrategy == numPrices+2))/sqrt(reps);

bestR = mean(revenue(priceStrategy == numPrices+1 & T == 0.73));
bestR_lt = mean(revenue(priceStrategy == numPrices+1 & T == 0.2));

bestR_ci = 1.96*std(revenue(priceStrategy == numPrices+1 & T == 0.73))/sqrt(reps);

figure, errorbar(r(1,idx)/bestR, r(2,idx)/bestR, ...
                 'b.');
hold on;
plot([0 numPrices+2], [optR/bestR, optR/bestR], 'b--');
plot([0 numPrices+2], [1, 1], 'b-.');
plot([0 numPrices+2], [1+bestR_ci/bestR, 1+bestR_ci/bestR], 'b:');
plot([0 numPrices+2], [1-bestR_ci/bestR, 1-bestR_ci/bestR], 'b:');
title(['Relative Cost, Sorted by Distance from Optimal ' ...
       'Price']);
%xlabel('Price Strategy (sorted by distance from optimal)')
xlabel('Price Strategy');
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Uniform', 'Time', 'Distance'});             % remove numeric labels
ylabel('Relative Cost');
%set(gca, 'XTickLabel', []);             % remove numeric labels

% figure, plot(h(1,idx)/best, r(1,idx)/bestR, 'b.', 'MarkerSize', 15)
% hold on;
% plot(best/best, bestR/bestR, 'rs', 'MarkerSize', 10);
% plot(opt/best, optR/bestR, 'k*', 'MarkerSize', 10);
% title('Price-Throughput Frontier');
% xlabel('Mean Relative Throughput');
% ylabel('Mean Relative Revenue');
% legend({'Other Prices', 'Optimal Price', 'Optimal Partition'})

figure, plot([h(1,:);h_lt(1,:)]/best, [r(1,:);r_lt(1,:)]/bestR, 'b:')
hold on;
hh(1) = plot([h(1,:)]/best, [r(1,:)]/bestR, 'b.', 'MarkerSize', 15)
hh(2) = plot([h_lt(1,:)]/best, [r_lt(1,:)]/bestR, 'bo', 'MarkerSize', 5)

plot([best,best_lt]/best, [bestR,bestR_lt]/bestR, 'r:', 'MarkerSize', 10);

hh(3) = plot(best/best, bestR/bestR, 'rs', 'MarkerSize', 10, ...
             'MarkerFaceColor', 'r')
hh(4) = plot(best_lt/best, bestR_lt/bestR, 'rs', 'MarkerSize', 10)

hh(5) = plot(opt/best, optR/bestR, 'k*', 'MarkerSize', 10);
title('Price-Throughput Frontier');
xlabel('Mean Relative Throughput');
ylabel('Mean Relative Revenue');
legend(hh, {'Other Prices - High T', 'Other Prices - Low T', ['Optimal ' ...
                    'Prices - High T'], 'Optimal Prices - Low T', 'Optimal Partition'});
set(gcf, 'Position', [60 200 700 600])

% --------------------------------------------------
% Flow Equity
% --------------------------------------------------
optR = mean(inequity(priceStrategy == numPrices+2));
optR_ci = 1.96*std(inequity(priceStrategy == numPrices+2))/sqrt(reps);

bestR = mean(inequity(priceStrategy == numPrices+1 & T == 0.73));
bestR_lt = mean(inequity(priceStrategy == numPrices+1 & T == 0.2));

bestR_ci = 1.96*std(inequity(priceStrategy == numPrices+1 & T == 0.73))/sqrt(reps);

figure, plot([h(1,:);h_lt(1,:)]/best, [I(1,:);I_lt(1,:)]/bestR, 'b:')
hold on;
hh(1) = plot([h(1,:)]/best, [I(1,:)]/bestR, 'b.', 'MarkerSize', 15)
hh(2) = plot([h_lt(1,:)]/best, [I_lt(1,:)]/bestR, 'bo', 'MarkerSize', 5)

plot([best,best_lt]/best, [bestR,bestR_lt]/bestR, 'r:', 'MarkerSize', 10);

hh(3) = plot(best/best, bestR/bestR, 'rs', 'MarkerSize', 10, ...
             'MarkerFaceColor', 'r')
hh(4) = plot(best_lt/best, bestR_lt/bestR, 'rs', 'MarkerSize', 10)

hh(5) = plot(opt/best, optR/bestR, 'k*', 'MarkerSize', 10);
title('Flow Inequity-Throughput Frontier');
xlabel('Mean Relative Throughput');
ylabel('Mean Relative Flow Inequity');
legend(hh, {'Other Prices - High T', 'Other Prices - Low T', ['Optimal ' ...
                    'Prices - High T'], 'Optimal Prices - Low T', 'Optimal Partition'});
set(gcf, 'Position', [60 200 700 600])


% --------------------------------------------------
% Wait Equity
% --------------------------------------------------
optR = mean(inequity2(priceStrategy == numPrices+2));
optR_ci = 1.96*std(inequity2(priceStrategy == numPrices+2))/sqrt(reps);

bestR = mean(inequity2(priceStrategy == numPrices+1 & T == 0.73));
bestR_lt = mean(inequity2(priceStrategy == numPrices+1 & T == 0.2));

bestR_ci = 1.96*std(inequity(priceStrategy == numPrices+1 & T == 0.73))/sqrt(reps);

figure, plot([h(1,:);h_lt(1,:)]/best, [I2(1,:);I2_lt(1,:)]/bestR, 'b:')
hold on;
hh(1) = plot([h(1,:)]/best, [I2(1,:)]/bestR, 'b.', 'MarkerSize', 15)
hh(2) = plot([h_lt(1,:)]/best, [I2_lt(1,:)]/bestR, 'bo', 'MarkerSize', 5)

plot([best,best_lt]/best, [bestR,bestR_lt]/bestR, 'r:', 'MarkerSize', 10);

hh(3) = plot(best/best, bestR/bestR, 'rs', 'MarkerSize', 10, ...
             'MarkerFaceColor', 'r')
hh(4) = plot(best_lt/best, bestR_lt/bestR, 'rs', 'MarkerSize', 10)

hh(5) = plot(opt/best, optR/bestR, 'k*', 'MarkerSize', 10);
title('Wait Inequity-Throughput Frontier');
xlabel('Mean Relative Throughput');
ylabel('Mean Relative Wait Inequity');
legend(hh, {'Other Prices - High T', 'Other Prices - Low T', ['Optimal ' ...
                    'Prices - High T'], 'Optimal Prices - Low T', 'Optimal Partition'});
set(gcf, 'Position', [60 200 700 600])

% --------------------------------------------------
% Length Equity
% --------------------------------------------------
optR = mean(inequity3(priceStrategy == numPrices+2));
optR_ci = 1.96*std(inequity3(priceStrategy == numPrices+2))/sqrt(reps);

bestR = mean(inequity3(priceStrategy == numPrices+1 & T == 0.73));
bestR_lt = mean(inequity3(priceStrategy == numPrices+1 & T == 0.2));

bestR_ci = 1.96*std(inequity(priceStrategy == numPrices+1 & T == 0.73))/sqrt(reps);

figure, plot([h(1,:);h_lt(1,:)]/best, [I3(1,:);I3_lt(1,:)]/bestR, 'b:')
hold on;
hh(1) = plot([h(1,:)]/best, [I3(1,:)]/bestR, 'b.', 'MarkerSize', 15)
hh(2) = plot([h_lt(1,:)]/best, [I3_lt(1,:)]/bestR, 'bo', 'MarkerSize', 5)

plot([best,best_lt]/best, [bestR,bestR_lt]/bestR, 'r:', 'MarkerSize', 10);

hh(3) = plot(best/best, bestR/bestR, 'rs', 'MarkerSize', 10, ...
             'MarkerFaceColor', 'r')
hh(4) = plot(best_lt/best, bestR_lt/bestR, 'rs', 'MarkerSize', 10)

hh(5) = plot(opt/best, optR/bestR, 'k*', 'MarkerSize', 10);
title('Length Inequity-Throughput Frontier');
xlabel('Mean Relative Throughput');
ylabel('Mean Relative Length Inequity');
legend(hh, {'Other Prices - High T', 'Other Prices - Low T', ['Optimal ' ...
                    'Prices - High T'], 'Optimal Prices - Low T', 'Optimal Partition'});
set(gcf, 'Position', [60 200 700 600])

% --------------------------------------------------
% Utilization Equity
% --------------------------------------------------
optR = mean(inequity4(priceStrategy == numPrices+2));
optR_ci = 1.96*std(inequity4(priceStrategy == numPrices+2))/sqrt(reps);

bestR = mean(inequity4(priceStrategy == numPrices+1 & T == 0.73));
bestR_lt = mean(inequity4(priceStrategy == numPrices+1 & T == 0.2));

bestR_ci = 1.96*std(inequity(priceStrategy == numPrices+1 & T == 0.73))/sqrt(reps);

figure, plot([h(1,:);h_lt(1,:)]/best, [I4(1,:);I4_lt(1,:)]/bestR, 'b:')
hold on;
hh(1) = plot([h(1,:)]/best, [I4(1,:)]/bestR, 'b.', 'MarkerSize', 15)
hh(2) = plot([h_lt(1,:)]/best, [I4_lt(1,:)]/bestR, 'bo', 'MarkerSize', 5)

plot([best,best_lt]/best, [bestR,bestR_lt]/bestR, 'r:', 'MarkerSize', 10);

hh(3) = plot(best/best, bestR/bestR, 'rs', 'MarkerSize', 10, ...
             'MarkerFaceColor', 'r')
hh(4) = plot(best_lt/best, bestR_lt/bestR, 'rs', 'MarkerSize', 10)

hh(5) = plot(opt/best, optR/bestR, 'k*', 'MarkerSize', 10);
title('Utilization Inequity-Throughput Frontier');
xlabel('Mean Relative Throughput');
ylabel('Mean Relative Utilization Inequity');
legend(hh, {'Other Prices - High T', 'Other Prices - Low T', ['Optimal ' ...
                    'Prices - High T'], 'Optimal Prices - Low T', 'Optimal Partition'});
set(gcf, 'Position', [60 200 700 600])


% --------------------------------------------------
% Mean Percent Difference
% --------------------------------------------------

mean((h(1,:)-h_lt(1,:))./h(1,:))*100
mean((I(1,:)-I_lt(1,:))./I(1,:))*100
mean((I2(1,:)-I2_lt(1,:))./I2(1,:))*100
mean((I3(1,:)-I3_lt(1,:))./I3(1,:))*100
mean((I4(1,:)-I4_lt(1,:))./I4(1,:))*100
