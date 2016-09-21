% Script to generate test data and optimal partitions to be
% compared to provably optimal solutions generated by PDQ

close all;
clear all;

type = [1, 1, 2];
numChains = length(find(type==1));
numCycles = length(find(type==2));

filenameroot = '2-1networks'

configs = dlmread([filenameroot, 'config.txt']);

NUM_TESTS = configs(1,1);
configs = configs(2:end,:);             % thow out first row

for i = 1:NUM_TESTS
    
    mu = configs(i,1:3);
    d = configs(i,4:end-1);
    d = zeros(size(d));                 % no travel times...
    N = configs(i,end);
    
    for METHOD = [1, 2, 4]
        [u,w,q,x,flag] = optimalAssignment(mu,type,d,N,METHOD);
        if flag == 1
            [Nvec, p] = routing(type, mu, d, x, N, METHOD);
        else    
            Nvec = ones(1,length(d));
        end
    
        dlmwrite([filenameroot, 'solution_', num2str(METHOD), '.txt'], ...
                 [Nvec, flag], '-append', 'delimiter', '\t');
    end
    
end


type = [1, 1, 2, 2];
numChains = length(find(type==1));
numCycles = length(find(type==2));

filenameroot = '2-2networks'

configs = dlmread([filenameroot, 'config.txt']);

NUM_TESTS = configs(1,1);
configs = configs(2:end,:);             % thow out first row

for i = 1:NUM_TESTS
    
    mu = configs(i,1:4);
    d = configs(i,5:end-1);
    d = zeros(size(d));                 % no travel times...
    N = configs(i,end);
    
    for METHOD = [1, 2, 4]
        [u,w,q,x,flag] = optimalAssignment(mu,type,d,N,METHOD);
        if flag == 1
            [Nvec, p] = routing(type, mu, d, x, N, METHOD);
        else    
            Nvec = ones(1,length(d));
        end
        dlmwrite([filenameroot, 'solution_', num2str(METHOD), '.txt'], ...
                 [Nvec, flag], '-append', 'delimiter', '\t');
    end
    
end


% type = [1, 1, 2, 2, 2];
% numChains = length(find(type==1));
% numCycles = length(find(type==2));

% filenameroot = '2-3networks'

% configs = dlmread([filenameroot, 'config.txt']);

% NUM_TESTS = configs(1,1);
% configs = configs(2:end,:);             % thow out first row

% for i = 1:NUM_TESTS
    
%     mu = configs(i,1:5);
%     d = configs(i,6:end-1);
%     d = zeros(size(d));                 % no travel times...
%     N = configs(i,end);
    
%     for METHOD = [1, 2, 4]
%         [u,w,q,x,flag] = optimalAssignment(mu,type,d,N,METHOD);
%         if flag == 1
%             [Nvec, p] = routing(type, mu, d, x, N, METHOD);
%         else    
%             Nvec = ones(1,length(d));
%         end
    
%         dlmwrite([filenameroot, 'solution_', num2str(METHOD), '_zeroD.txt'], ...
%                  [Nvec, flag], '-append', 'delimiter', '\t');
%     end
    
% end


