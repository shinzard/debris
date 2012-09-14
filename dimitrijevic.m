function dimitrijevic

days = unique(floor(t601)');
for i = days
    idx = find(floor(t601)==i);
    
    if length(idx)>150
        
        % Segment lengths (in hours)
        S = diff(t601(idx)*24);
        
        % Arrivals per segment
        k = zeros(1,length(idx));
        
        % Queue length right before ith departure
        q = zeros(1,length(idx));
       
        kr = k;
        kr1 = ;
        kr(1) = 1;
        for r = 2:length(idx)
           
           
        end
    end
end


function prob = p(k, S)
    n = length(S);
    prod = 1;
    
    for i = 1:n
        prod = prod * S(i)^k(i)/factorial(k);
    end
    prob = prod*factorial(n)/sum(S)^n;
    
function l = F(k, S)
    l = sum(log(factorial(k))-k.*log(S));
        
function difference = delF(k, i, S)
    kp = k;
    kp(i) = kp(i) + 1;
    difference = F(k,S) - F(k,S)
    log(kp(i)+1) - log(S(i))