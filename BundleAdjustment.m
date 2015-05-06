function [Cset, Rset, Xset] = BundleAdjustment(K, Cset, Rset, X, ReconX, V, Mx_bundle, My_bundle)

opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'MaxIter', 1e64, 'Display', 'iter');

I = length(Cset);
J = length(X(:,1));

params = zeros(I*7+3*J,1);

for i=1:I
    params((i-1)*7+1) = Cset{i}(1);
    params((i-1)*7+2) = Cset{i}(2);
    params((i-1)*7+3) = Cset{i}(3);
    
    q = R2q(Rset{i});
    
    params((i-1)*7+4) = q(1);
    params((i-1)*7+5) = q(2);
    params((i-1)*7+6) = q(3);
    params((i-1)*7+7) = q(4);
end

for i=1:J
    params(I*7+(i-1)*3+1) = X(i,1);
    params(I*7+(i-1)*3+2) = X(i,2);
    params(I*7+(i-1)*3+3) = X(i,3);
end

optparams = lsqnonlin(@(params) func(params, K, Mx_bundle, My_bundle, V, I, J), params, [], [], opts);

[Cset, Rset, Xset] = decodeParams(optparams, I, J);

function f = func(params, K, Mx_bundle, My_bundle, V, I, J)

[Cset, Rset, X] = decodeParams(params, I, J);

f = zeros(I*J*2,1);

for i = 1:I
    
    P = K*Rset{i}*[eye(3) -1.*Cset{i}];
    
    for j = 1:J
        Xhom = [X{j}' 1]';
        f(((i-1)*I+j-1)*2+1:((i-1)*I+j-1)*2+2) = [V(j,i)*(Mx_bundle(j,i) - (P(1,:)*Xhom)/(P(3,:)*Xhom));
                                                  V(j,i)*(My_bundle(j,i) - (P(2,:)*Xhom)/(P(3,:)*Xhom))];
    end
end

function [Cset, Rset, X] = decodeParams(p, I, J)

Cset = cell([1 I]);
Rset = cell([1 I]);
X = cell([1 J]);

for i=1:I
    Cset{i} = p((i-1)*7+1:(i-1)*7+3);
    Rset{i} = q2R(p((i-1)*7+4:(i-1)*7+7));
end

for j=1:J
    X{i} = p(7*I+(j-1)*3+1:7*I+(j-1)*3+3);
end

function r = quat2mat(q)

r = [1-2*q(3)^2-2*q(4)^2 2*q(2)*q(3)-2*q(4)*q(1) 2*q(2)*q(4)+2*q(3)*q(1);
     2*q(2)*q(3)+2*q(4)*q(1) 1-2*q(2)^2-2*q(4)^2 2*q(3)*q(4)-2*q(2)*q(1);
     2*q(2)*q(4)-2*q(3)*q(1) 2*q(3)*q(4)+2*q(2)*q(1) 1-2*q(2)^2-2*q(3)^2];