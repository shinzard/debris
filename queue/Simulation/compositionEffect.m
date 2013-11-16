% This function maps team composition to resulting mean service
% rate
% 
% 4 Nov 2013
% J.Brooks
%  
%  Last Modified: 4 Nov 2013
%  
function mu = compositionEffect(familiarity, size, baseMu, param)
    mu = 1/10^( param.deltaA*(0.736 - familiarity) + ...
                param.deltaN*(log10(3.31/size)) )*baseMu;
    if isnan(mu)
        familiarity
        baseMu
        size
    end