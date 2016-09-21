% This function maps team composition to resulting mean service
% rate
% 
% 4 Nov 2013
% J.Brooks
%  
%  Last Modified: 4 Nov 2013
%  
function mu = compositionEffect(familiarity, size, baseMu, param)

size = max(1,size);                     % for teams with nobody
                                        % currently...

mu = 10^( param.deltaA*(0.736 - familiarity) + ... % it's fluidity!
         param.deltaN*(log10(3.31/size)) )*baseMu;

mu = max(1,min(mu, 15));                % set reasonable limits
                                        % based on extent of
                                        % validation data
%mu = baseMu;    