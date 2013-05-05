
close all;
%proj67ewma;                             % get d structure

decIdx = find(~isnan(d.from));

temps = linspace(0.01,10);               % cf. 0.25 in Gonzalez2003

likelihood = zeros(1,length(temps));

for i = 1:length(temps)
    for j = 1:length(decIdx)
        tmp = d.ewma(:,decIdx(j));
        tmpH = d.haul(:,decIdx(j));
        tmpS = d.size(:,decIdx(j));
        Ltmp = zeros(6,6);
        price = (tmpH>=15)*1.2+(tmpH<15); % current policy...
                                          %price = ones(6,1);
        
        for k = 1:6
            x = price./tmp - price(k)/tmp(k);
            Ltmp(k,:) = sign(x).* x.^2; % notice opposite
                                        % from other analysis
                                        % because looking at inverses
        end
         
        Ltmp = Ltmp + diag(NaN.*ones(1,6));% remove diagonal...
 
        probs = exp(Ltmp/temps(i))/sum(nansum(exp(Ltmp/temps(i))));
        
        if ~isnan(probs(d.from(decIdx(j)), d.to(decIdx(j))))
            likelihood(i) = likelihood(i) + probs(d.from(decIdx(j)), ...
                                                  d.to(decIdx(j)));
        end
    end
end

figure, plot(temps, likelihood, 'b.')
title('Total likelihood')
xlabel('Temperature');
ylabel('Sum of likelihoods');

[mL,i] = max(likelihood);
disp(sprintf(['Best T: %2.2f; Total Likelihood: %2.2f, Relative to ' ...
              'Random: %2.2f'], temps(i),mL,mL/4.2));