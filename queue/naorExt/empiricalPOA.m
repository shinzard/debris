clear all;

numTests = 1000;

Crange = [1, 5];
Rrange = [5, 15];
murange = [0.1 10];

indivRate = zeros(1,length(numTests));
thru = zeros(1,length(numTests));
socialRate = zeros(1,length(numTests));
socialThru = zeros(1,length(numTests));
ns = zeros(1,length(numTests));

Rvec = zeros(1,length(numTests));
Cvec = zeros(1,length(numTests));
mus = zeros(length(numTests),2);

for i = 1:numTests
    C = rand*diff(Crange) + min(Crange);
    mu1 = rand*diff(murange) + min(murange);
    mu2 = rand*diff(murange) + min(murange);
    minR = C*(1/mu1+1/mu2);
    R = rand*diff(Rrange) + min([Rrange, minR]);

    mu = [mu1, mu2];

    [tmp, ns(i), indivRate(i), thru(i), socialRate(i), socialThru(i)] ...
        = explore(C,R,0,mu,0);
    Rvec(i) = R;
    Cvec(i) = C;
    mus(i,:) = mu;
    i = i + 1
end

close all;

hist(socialThru./thru);
set(gca, 'fontsize', 20);
title('System Throughput Ratio');
xlabel('social/indiv');
ylabel('Count');

figure;
hist(socialRate./indivRate);
title('Price of Anarchy');
xlabel('social/indiv');
ylabel('Count');
