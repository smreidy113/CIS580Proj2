function F = EstimateFundamentalMatrix(x1, x2)

n = length(x1);

A = zeros(n,9);

for i=1:n
    A(i,:) = [x2(i,1)*x1(i,1) x2(i,1)*x1(i,2) x2(i,1) x2(i,2)*x1(i,1) x2(i,2)*x1(i,2) x2(i,2) x1(i,1) x1(i,2) 1];
end

[~,~,V] = svd(A);
fcol = V(:,end);

F = [fcol(1) fcol(2) fcol(3);
     fcol(4) fcol(5) fcol(6);
     fcol(7) fcol(8) fcol(9)];
 
[U,D,V] = svd(F);
D(3,3) = 0;
F = U*D*V';
