
close all;
FORCE_RELOAD = 0;

if ~exist('d') || FORCE_RELOAD
    proj67ewma;                             % get d structure
end

decIdx = find(~isnan(d.from));

temps = linspace(0.01,3);               % cf. 0.25 in Gonzalez2003
prices = 1; %linspace(0.01,100,200);


bestT = zeros(1,length(prices));
mL = zeros(1,length(prices));
cnt = 1;

for p = prices
    likelihood = zeros(1,length(temps)); % need to reset!!
    
    for i = 1:length(temps)
        
        for j = 1:length(decIdx)
            tmp = d.ewma(:,decIdx(j));
            tmpH = d.haul(:,decIdx(j));
            tmpS = d.size(:,decIdx(j));
            Ltmp = zeros(6,6);
            price = p*((tmpH>=15)*1.2+(tmpH<15)); % current policy...
                                                  %price = ones(6,1);
            
            for k = 1:6
                x = log((price./tmp)./(price(k)/tmp(k))); % BEST
                %x = ((price./tmp)./(price(k)/tmp(k))).^2;
                %x = log10((price./tmp)./(price(k)/tmp(k)));
                %x = (price./tmp)./(price(k)/tmp(k));
                %x = price./tmp - price(k)/tmp(k);
                %x = tmp(k) - tmp;
                %x = 1./tmp - 1/tmp(k);
                %x = tmpH./tmpH(k);
                %x = (price.*tmpS./tmp)./(price(k)*tmpS(k)/tmp(k));
                %x = (tmpH./tmp)./(tmpH(k)/tmp(k));
                Ltmp(k,:) = x;         

            end
            
            Ltmp = Ltmp + diag(NaN.*ones(1,6));% remove diagonal...
            
            probs = exp(Ltmp/temps(i))/sum(nansum(exp(Ltmp/temps(i))));
            
            if ~isnan(probs(d.from(decIdx(j)), d.to(decIdx(j))))
                likelihood(i) = likelihood(i) + probs(d.from(decIdx(j)), ...
                                                      d.to(decIdx(j)));
            end
        end
    end
    
    if length(prices) == 1
        figure, plot(temps, likelihood, 'b.')
        title(sprintf('Total likelihood: %2.2f',p))
        xlabel('Temperature');
        ylabel('Sum of likelihoods');
    end
    
    [mL(cnt),i] = max(likelihood);
    bestT(cnt) = temps(i);
    disp(sprintf(['P: %2.2f; Best T: %2.2f; Total Likelihood: %2.2f, Relative to ' ...
                  'Random: %2.2f'], p, temps(i),mL(cnt),mL(cnt)/4.2));
    cnt = cnt + 1;
    
    
end

if length(prices) > 1
    figure, subplot(211), plot(prices, bestT, 'b.');
    subplot(212), plot(prices, mL, 'b.');
end