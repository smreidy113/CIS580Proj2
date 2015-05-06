function [C, R] = NonlinearPnP(X, x, K, C, R)

opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'MaxIter', 1e64, 'Display', 'iter');

cq0 = [C' R2q(R)'];

cqo = lsqnonlin(@(cq) func(cq, K, X, x), cq0, [], [], opts);

C = cqo(1:3)';
q = cqo(4:7)';
R = q2R(q);

function f = func(cq, K, X, x)

C = cq(1:3)';
q = cq(4:7)';
if numel(q) == 1
    q = q';
end
R = q2R(q);

j = length(X(:,1));

f = zeros(j*2,1);

P = K*[R C];

for i=1:j

    Xhom = [X(j,:) 1]';
    f((i-1)*2+1:(i-1)*2+2) = [x(j,1)-(P(1,:)*Xhom)/(P(3,:)*Xhom);
         x(j,2)-(P(2,:)*Xhom)/(P(3,:)*Xhom)];

end

function r = quat2mat(q)

r = [1-2*q(3)^2-2*q(4)^2 2*q(2)*q(3)-2*q(4)*q(1) 2*q(2)*q(4)+2*q(3)*q(1);
     2*q(2)*q(3)+2*q(4)*q(1) 1-2*q(2)^2-2*q(4)^2 2*q(3)*q(4)-2*q(2)*q(1);
     2*q(2)*q(4)-2*q(3)*q(1) 2*q(3)*q(4)+2*q(2)*q(1) 1-2*q(2)^2-2*q(3)^2];