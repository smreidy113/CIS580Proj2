function [C, R] = LinearPnP(X, x, K)

A = zeros(length(X(:,1))*3, 12);

for i=1:length(X(:,1))
   
    Xhom = [X(i,:) 1];
    
    A((i-1)*3+1,:) = [0 0 0 0 Xhom x(i,2).*Xhom];
    A((i-1)*3+2,:) = [Xhom 0 0 0 0 -x(i,1).*Xhom];
    A((i-1)*3+3,:) = [-x(i,2).*Xhom x(i,1).*Xhom 0 0 0 0];
    
end

[~,~,V] = svd(A);

Vcol = V(:,end)./V(end,end);
P = [Vcol(1) Vcol(2) Vcol(3) Vcol(4);
     Vcol(5) Vcol(6) Vcol(7) Vcol(8);
     Vcol(9) Vcol(10) Vcol(11) Vcol(12)];

 Rt = K\P;

R = Rt(1:3,1:3);
t = Rt(:,4);

[U,D,V] = svd(R);

R = U*V';

if det(R) <= 0
    R = -1.*R;
end

C = -R'*t;