% This script calculates team familiarity metric based on prior
% working history of team members (i.e., trucks); extracted from
% fluidity3.m 
% 
%  J.Brooks
%  21 Oct 2013
%  Last Modified: 21 Oct 2013
%  
function [fam,n,cumHist,today] = familiarity(t,lastT,oldCumhist,oldToday,class,oldFam)

c = [1:length(oldFam)];                 % teams
n = zeros(1,length(c));                 % team size
teamHist = zeros(1,length(c));          % team cumulative history

today = oldToday;
cumHist = oldCumhist;
fam = oldFam;

% Team sizes
for i = 1:length(c)
    n(i) = length(find(class == c(i)));
end

if mod(lastT,12) < mod(t, 12)           % within same day
    for i = 1:length(c)
        idx = find(class == c(i));
        for j = 1:length(idx)
           for k = j+1:length(idx)
               today(idx(j), idx(k))= 1;
               today(idx(k), idx(j))= 1; % necessary only for
                                         % assignment problem (in
                                         % dispatcher.m) 
           end
        end
    end

% Update familiarity 
else                                    % new day
    % Accumulate today's history and reset
    cumHist = cumHist + today;
    today = zeros(size(today));
    
    % Sum and scale by number of pairs between members and duration
    for i = 1:length(c)
        idx = find(class == c(i));
        for j = 1:length(idx)
           for k = j+1:length(idx)
               teamHist(i) = teamHist(i) + cumHist(idx(j), idx(k));
           end
        end
        if (n(i)>1) && (t>12)
            fam(i) = teamHist(i)/((n(i)*(n(i)-1))/2)/floor(t/12);
        else
            fam(i) = 1;
        end
    end

end

