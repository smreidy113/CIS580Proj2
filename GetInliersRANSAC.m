function [y1, y2, idx] = GetInliersRANSAC(x1, x2)

M = 1000;
N = length(x1);
e = 0.005;

n_rands = 8;

n = 0;
Sprev = [];

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
    
    x1d = zeros(n_rands,2);
    x2d = zeros(n_rands,2);
    
    for j=1:n_rands
        x1d(j,:) = x1(inds(j),:);
        x2d(j,:) = x2(inds(j),:);
    end
    
    %disp(size(x1d));
    %disp(size(x2d));
    
    F = EstimateFundamentalMatrix(x1d, x2d);
    
    S = [];
    
    for j=1:N
        if norm([x2(j,:) 1] * F * [x1(j,:)'; 1]) < e
            S = [S j];
        end
    end
    
    if n < length(S)
        n = length(S);
        Sprev = S;
    end
    
end

y1 = zeros(n,2);
y2 = zeros(n,2);

for i=1:n
    y1(i,:) = x1(Sprev(i),:);
    y2(i,:) = x2(Sprev(i),:);
    idx = Sprev(i);
end