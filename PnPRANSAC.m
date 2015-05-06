function [C, R] = PnPRANSAC(X, x, K)

M = 1000;
er = 10;
N = length(X(:,1));

n_rands = 6;

n = 0;
Sprev = [];

C = [0 0 0]';
R = eye(3);

for i=1:M
    
    inds = zeros(1,n_rands);
    
    j = 0;
    
    while j < n_rands
        
        r = floor(rand*N);
        
        if ~ismember(r, inds)
            inds(j+1) = r;
            j = j + 1;
        end
        
    end
    
    Xd = zeros(n_rands, 3);
    xd = zeros(n_rands, 2);
    
    for j=1:n_rands
        Xd(j,:) = X(inds(j),:);
        xd(j,:) = x(inds(j),:);
    end
    
    [C1,R1] = LinearPnP(Xd, xd, K);
    
    S = [];
    
    P = K*[R1 -R1*C1];
    
    for j = 1:N
        
        Xhom = [X(j,:) 1]';
        
        e = (x(j,1) - (P(1,:)*Xhom)/(P(3,:)*Xhom))^2 + (x(j,2) - (P(2,:)*Xhom)/(P(3,:)*Xhom))^2;
        
        if e < er
            S = [S j];
        end
        
    end
        
    if n < length(S)
        C = C1;
        R = R1;
        n = length(S);
        Sprev = S;
    end
        
end