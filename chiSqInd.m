function p = chiSqInd(count)
% CHISQIND - This function performs a chi-square independence test
% using a contingency table as input
%   
% 31 August 2011
% J.Brooks
%   
% H0: row and column methods of classification are independent
% H1: methods are dependent
% 
% Assumptions: independence of samples

    n = sum(sum(count));

    % Category 1
    u_hat = 1/n*sum(count');
    
    % Category 2
    v_hat = 1/n*sum(count);
    
    % Expected counts
    expected =  n*u_hat'*v_hat;
    
    % Test Statistic
    x = sum(sum(((count- expected).^2)./expected));
    
    % P-value
    p = 1 - chi2cdf(x, (size(count,1)-1)*(size(count,2)-1));
    
    if x > chi2inv(0.95, (size(count,1)-1)*(size(count,2)-1))
        disp(sprintf(['Reject H0 (not independent); P-value: %2.6f; ' ...
                      'Stat: %2.6f'], p, x));
    else
        disp(sprintf(['Fail to reject H0 (may be dependent); P-value: ' ...
        '%2.6f; Statistic: %2.6f'], p, x));       
    end
