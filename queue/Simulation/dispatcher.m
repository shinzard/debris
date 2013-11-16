% This function defines the dispatcher decision logic. It is called
% from simulateQ at random intervals as specified based on current
% system state and the parameters.
% 
% 4 Nov 2013
% J.Brooks
%  
%  Last Modified: 12 Nov 2013
%  
function newAssign = dispatcher(oldAssign, param, state)
n = length(oldAssign);    

% Fixed
switch(param.method)
  case 'Fixed'                          % FIXED
    newAssign = oldAssign;
  case 'Bounded'                        % Bounded Rationality
    ak = zeros(1,n*(n-1));              % all possible switches
    
    p = exp(ak/param.T)./sum(exp(ak/param.T));
    newAssign = oldAssign;

  case 'Random'                         % Random assignment
    
end
        