% SCRIPT to compare partition methodology with empirically optimal
% solutions (originially for IEEE HST conference
% 
%  9 June 2012
%  J.Brooks

load /home/james/Documents/MATLAB/queueing/27Oct2012_InitialStudy.mat


type = [1 2 2];
partition = zeros(10,2);

for i = 1:NUM_TESTS
    mus = [mu3r(i), mu1, mu2r(i)];
    [u,w,q,x,e]=optimalAssignment(mus,type,0,N,1);
    [Nvec, p] = routing(type, mus, 0, x, N, 1);
    partition(i,:) = Nvec;
end


percentCorrect = zeros(1,11);
for i = 1:11
    percentCorrect(i) = length(find(abs(partition(:,1) - n1)<=i-1))/NUM_TESTS;
end

figure, plot([0:10], percentCorrect, 'b');
hold on, plot([0:10], percentCorrect, 'b.');
title('Percent of Solutions by Distance from Optimal');
xlabel('Partition Difference');
ylabel('Percent');