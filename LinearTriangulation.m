function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)

P1 = K*R1*[eye(3) -1.*C1];
P2 = K*R2*[eye(3) -1.*C2];

X = zeros(length(x1), 3);

for i=1:length(x1)

    x1p = [x1(i,:) 1]';
    x2p = [x2(i,:) 1]';
    
    skew1 = Vec2Skew(x1p);
    skew2 = Vec2Skew(x2p);
    
    A = [skew1*P1; skew2*P2];
    [~,~,v] = svd(A);
    X(i,:) = v(1:3,end)./v(end,end);

end

function skew = Vec2Skew(v)
skew = [0 -v(3) v(2);
        v(3) 0 -v(1);
        -v(2) v(1) 0];