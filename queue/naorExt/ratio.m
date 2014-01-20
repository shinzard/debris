Tvec = linspace(0,20,200);
Pvec = linspace(0.01,0.99,20);
r = zeros(length(Pvec),length(Tvec));
i = 1;
for p = Pvec
    %    disp(sprintf('p = %f',p));
    j = 1;
    for T = Tvec
        %        disp(sprintf('T = %f',T));
        [tmp,tmp,tmp,tmp,r(i,j)] = explore(1,5,T,p);
        j = j + 1;
    end
    i = i + 1;
end

mesh(Tvec, Pvec, r)


