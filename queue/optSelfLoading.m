% This script is a quick proof of concept for using queueing theory
% to make TDSR design and policy decisions regarding the use of
% self-loading trucks in the debris removal mission
%
% 12 Feb 2013
% J.Brooks

close all;
clear all;
pRange = [0:0.01:1];
MOVIE = 0;
L_RANGE = [1:30];                       % incoming rate (per hour)
C = 3;                                  % number of simulataneous
                                        % self-loading spots
mu1range = [1:10];                      % service rate of self-loaders
mu2range = [1:20];                      % service rate of chipper
% Particular case:
%mu1range = 10
%mu2range = 10
prop = NaN*zeros(length(mu1range),length(mu2range),length(L_RANGE));
equality =NaN*zeros(length(mu1range),length(mu2range),length(L_RANGE));

i = 1;

for mu_s = mu1range
    j = 1;
    for mu_c = mu2range
        idx = 1;
        for lambda = L_RANGE
            stable = find((((1-pRange)*lambda/mu_c)<1 & ...
                          (pRange*lambda/(C*mu_s))<1) == 1);
            if length(stable) < 2
                continue;
            end
            waitC = 1./(mu_c - (1-pRange)*lambda);

            % from Kleinrock1975, p.102, 404
            rho = pRange*lambda/mu_s;
            for k = 1:length(rho)
                summation(k) = sum(rho(k).^[0:C-1]./ ...
                                   factorial([0:C-1]));
            end
            p0 = 1./( summation + ...
                      (rho.^C./factorial(C)).*(1./(1-rho/C)) );
            waitS = (1./(C*mu_s - pRange*lambda)).*...
                    (rho.^C./(factorial(C)*(1-rho/C))).*...
                    p0 + 1/mu_s;        % from Buzacott1993, p.79
                                        % note that \rho = \rho/C
            wait = (1-pRange).*waitC + pRange.*waitS;
            %            figure, plot(pRange,p0), hold on, plot(pRange, summation, 'r'), ...
            %                title('diagnostics');
            %            figure, plot(pRange(stable),waitS(stable), 'r');
            %            hold on;
            %            plot(pRange(stable),waitC(stable), 'g');
            %            plot(pRange(stable),wait(stable), 'b');
            [minWait, idxT] = min(wait(stable));

            prop(i,j,idx) = pRange(stable(idxT));
            equality(i,j,idx) = waitC(stable(idxT))/waitS(stable(idxT));
            idx = idx + 1;
        end
        j = j + 1;
    end
    i = i + 1;
end

if MOVIE
    aviobj = avifile('tdsrOptSL.avi', 'FPS', 1);
end

for i = 1:length(L_RANGE)
    figure, surfc(mu2range,mu1range,prop(:,:,i))
    title(['Arrival Rate (\lambda): ', sprintf('%d trucks/hr', L_RANGE(i))]);
    xlabel('Chipper Rate (trucks/hr)');
    ylabel('Self-Unloading Rate (trucks/hr)');
    zlabel('Optimal Proportion of Self-Load');
    axis([0,max(mu2range),0,max(mu1range),0,1]);
    %    print('-deps', sprintf('graphics/optSelf/prop%d.eps', i));
    if MOVIE
        F = getframe(gcf);
        aviobj = addframe(aviobj,F);
    end
end

if MOVIE
    aviobj = close(aviobj);

    aviobj = avifile('tdsrOptSL_EQ.avi', 'FPS', 1);
end

for i = 1:length(L_RANGE)
    figure, surfc(mu2range,mu1range,equality(:,:,i))
    title(['Arrival Rate (\lambda): ', sprintf('%d trucks/hr', L_RANGE(i))]);
    xlabel('Chipper Rate (trucks/hr)');
    ylabel('Self-Unloading Rate (trucks/hr)');
    zlabel('Optimal Proportion of Self-Load');
    axis([0,max(mu2range),0,max(mu1range),0,5]);
    %    print('-deps', sprintf('graphics/optSelf/equal%d.eps', i));
    if MOVIE
        F = getframe(gcf);
        aviobj = addframe(aviobj,F);
    end
end

if MOVIE
    aviobj = close(aviobj);
end


for i = [1:max(mu1range)]
    figure;
    for j = [1:max(mu2range)]
        plot(L_RANGE, reshape(prop(i,j,:),1,30)), hold on;
        xlabel('Arrival Rate');
        ylabel('Optimal Proportion of Self-Load');
    end
    title(sprintf('Self-loader rate: %d', mu1range(i)));
    %    print('-deps', sprintf('graphics/optSelf/flat%d.eps', i));
end