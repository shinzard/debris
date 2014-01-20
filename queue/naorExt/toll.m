Tvec = linspace(0,20,200);
Pvec = linspace(0.01,0.99,20);
maxToll = zeros(length(Pvec),length(Tvec));
maxSToll = zeros(length(Pvec),length(Tvec));
S = zeros(length(Pvec),length(Tvec));
I = zeros(length(Pvec),length(Tvec));
thruRatio = zeros(length(Pvec),length(Tvec));
i = 1;
for p = Pvec
    %    disp(sprintf('p = %f',p));
    j = 1;
    for T = Tvec
        %        disp(sprintf('T = %f',T));
        [maxToll(i,j),maxSToll(i,j), S(i,j), I(i,j), thruRatio(i,j)] = explore(3,5,T,p);
        j = j + 1;
    end
    i = i + 1;
end

close all;

mesh(Tvec, Pvec, maxToll)
title('Total Gain Rate');

figure;
mesh(Tvec, Pvec, maxSToll)
title('Maximum Social Toll');

figure;
mesh(Tvec, Pvec, S)
title('Social Optimal Threshold');

figure;
mesh(Tvec, Pvec, I)
title('Individual Threshold');

figure;
mesh(Tvec, Pvec, thruRatio)
title('Throughput Ratio');

figure;
[x,y]=max(maxToll');
plot(Pvec, Tvec(y));
title('Total Gain Rate');
