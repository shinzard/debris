R = 5;
Cvec = linspace(0.01, 4, 100);
POA = zeros(1,length(Cvec));

i = 1;

for C = Cvec
    [NE(i), POA(i)] = exploreNE([C C],R, min(1000,round(100*1/C)), 0);
    i = i + 1
end

figure;
plot(Cvec, NE, 'b');
title('Number of NEs');
xlabel('Common cost');

figure;
plot(Cvec, POA, 'b');
title('Price of Anarchy');
xlabel('Common cost');
