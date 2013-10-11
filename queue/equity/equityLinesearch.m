function lambda = equityLinesearch(mus, prop, NUM_ENTITIES, verb)
% This script finds the minimum wait-time inequity flow allocation through
% a simple linesearch (assuming no travel delay...).
%
% 24 Sept 2013
% J.Brooks
%

% load simplexResults first and pass in....
if nargin < 4
    verb = 0;
end

stepSize = 0.00001;
lambda = 0;
alpha = min(mus);
num = 0;

if length(prop) > 1                     % fixed proportion (length
                                        % or wait time equity)
    while (num < NUM_ENTITIES)
        lambda = lambda + stepSize;
        num = lambda/(mus(1)-lambda) + sum(lambda.*prop./(mus(2:end)-lambda.*prop));
    end
else                                    % equal wait times
   while (num < NUM_ENTITIES && all(mus(2:end)-alpha >= 0) && ...
          sum(mus(2:end)-alpha) < mus(1)) 
       mus(2:end)-alpha;
       alpha = alpha - stepSize;
       num = sum(mus(2:end)-alpha)/(mus(1)-sum(mus(2:end)-alpha)) + sum((mus(2:end)-alpha)./alpha);
   end
   
   % if no equitable solution can be found
   if (num > NUM_ENTITIES + 0.1 || any(mus(2:end)-alpha < 0) || ...
       sum(mus(2:end)-alpha) >= mus(1))
       lambda = NaN;
       return;
   else 
       lambda = sum(mus(2:end)-alpha);
   end
end

if verb
    if length(prop)>1
        disp(sprintf('Proportions: %s', num2str(prop)));
        disp(sprintf('Indiv. Flows: %s', num2str(lambda.*prop)));
        disp(sprintf('Server load: %s', num2str((lambda.*prop)./mus(2:end))));
        disp(sprintf('Wait times: %s', num2str(1./(mus(2:end) - lambda.*prop))));
        disp(sprintf('Server Lengths: %s', num2str(lambda.*prop./(mus(2:end) - lambda.*prop))));
    else
        disp(sprintf('alpha: %2.2f', alpha));
        disp(sprintf('Indiv. Flows: %s', num2str(mus(2:end)-alpha)));
        disp(sprintf('Server load: %s', num2str((mus(2:end)-alpha)./mus(2:end))));
        disp(sprintf('Wait times: %s', num2str(1./(mus(2:end) - (mus(2:end)-alpha)))));
        disp(sprintf('Server Lengths: %s', num2str((mus(2:end)-alpha)./alpha)));
    end
    disp(sprintf('Total Flow: %s', num2str(lambda)));
    
    
    disp(sprintf('Total Num: %s', num2str(num)));
end
