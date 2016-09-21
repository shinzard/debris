function lambda = equityLinesearch2(type, mus, dij, pr, prop, NUM_ENTITIES, verb)
% This script finds the minimum wait-time inequity flow allocation through
% a simple linesearch (assuming no travel delay...).
%
% 24 Sept 2013
% J.Brooks
%

% load simplexResults first and pass in....
if nargin < 7
    verb = 0;
end

stepSize = 0.0001;
lambda = 0;
alpha = min(mus);
num = 0;

central = find(type == 1);
parallel = find(type == 2);
%m = length(central);
%n = length(parallel);

if length(prop) > 1                     % fixed proportion (length
                                        % or wait time equity)
    while (num < NUM_ENTITIES)
        lambda = lambda + stepSize;
        centralFlows = (lambda.*prop)*pr;

        num = sum(centralFlows./(mus(central)-centralFlows)) + ... % central
              sum((lambda.*prop)./(mus(parallel)-lambda.*prop)) + ... % parallel
              sum(lambda.*prop*(pr.*dij));    % travel
        
    end
else                                    % equal wait times
   while (num < NUM_ENTITIES && all(mus(parallel)-alpha >= 0) && ...
          all((mus(parallel)-alpha)*pr < mus(central))) 

       alpha = alpha - stepSize;

       centralFlows = (mus(parallel)-alpha)*pr;

       num = sum(centralFlows./(mus(central)-centralFlows)) + ... % central
             sum((mus(parallel)-alpha)./alpha) + ... % parallel
             sum((mus(parallel)-alpha)*(pr.*dij));    % travel
   end
   
   % if no equitable solution can be found
   if (num > NUM_ENTITIES + 0.1 || any(mus(parallel)-alpha < 0) || ...
       any((mus(parallel)-alpha)*pr >= mus(central)))
       lambda = NaN;
       return;
   else 
       lambda = sum(mus(parallel)-alpha);
   end
end

if verb
    if length(prop)>1
        disp(sprintf('Proportions: %s', num2str(prop)));
        disp(sprintf('Indiv. Flows: %s', num2str(lambda.*prop)));
        disp(sprintf('Server load: %s', num2str((lambda.*prop)./mus(parallel))));
        disp(sprintf('Wait times: %s', num2str(1./(mus(parallel) - lambda.*prop))));
        disp(sprintf('Parallel server Lengths: %s', ...
                     num2str(lambda.*prop./(mus(parallel) - lambda.*prop))));
        disp(sprintf('Central server Lengths: %s', ...
                     num2str((lambda.*prop*pr)./(mus(central) - (lambda.*prop)*pr))));
        disp(sprintf('Total traveling: %s', num2str(lambda.*prop*dij)));
    else
        disp(sprintf('alpha: %2.2f', alpha));
        disp(sprintf('Indiv. Flows: %s', num2str(mus(parallel)-alpha)));
        disp(sprintf('Server load: %s', num2str((mus(parallel)-alpha)./mus(parallel))));
        disp(sprintf('Wait times: %s', num2str(1./(mus(parallel) - (mus(parallel)-alpha)))));
        disp(sprintf('Server Lengths: %s', num2str((mus(parallel)-alpha)./alpha)));
    end
    disp(sprintf('Total Flow: %s', num2str(lambda)));
    
    
    disp(sprintf('Total Num: %s', num2str(num)));
end
