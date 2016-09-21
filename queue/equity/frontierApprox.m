clear all;
close all;

load('simplexResults-5_01-Oct-2013 13:21:42')
simplexPlot

hold on;

w1 = [1/3, 1/3, 1/3]

lambda = equityLinesearch(mus, w1, NUM_ENTITIES);% find zero
                                                 % inequity
                                                 % solution

weights = 4/(4-lambda)^2 + mus(2:4)./(mus(2:4)-lambda/3).^2;% from
                                                            % optimality
                                                            % conditions!

w2 = weights./sum(weights);             % normalize (put on simplex)

figure(7), hold on;
figure(3), hold on;

for alpha = linspace(0,1)
    wtmp = w1*alpha+w2*(1-alpha);
    [u,w,q,x,e]=optimalAssignmentWeights(mus, type, 0, wtmp, NUM_ENTITIES, 1);
    figure(7), plot(sum(x), sqrt(sum(((x./mean(x)-1).^2)))/sqrt(6), 'bs')
    figure(3), plot3(wtmp(1), wtmp(2), wtmp(3), 'bs');
end


% --------------------------------------------------
% Iterative method with small \beta
% --------------------------------------------------

for i = 1:200
    [u,w,q,x,e]=optimalAssignmentWeights(mus, type, 0, w1, NUM_ENTITIES, 1);
    i = sqrt((x./mean(x) - 1).^2);
    n = (x<mean(x))*2-1;
    w1 = w1 + 0.01.*i.*n;
    w1 = w1/sum(w1);
    figure(7), plot(sum(x), sqrt(sum(((x./mean(x)-1).^2)))/sqrt(6), 'ro');
    figure(3), plot3(w1(1), w1(2), w1(3), 'ro');
end

