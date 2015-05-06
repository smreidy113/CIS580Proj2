function [Cset Rset] = ExtractCameraPose(E)

W = [0 -1 0;
     1 0 0;
     0 0 1];
 
[U, D, V] = svd(E);

C1 = U(:,3);
C2 = -1.*U(:,3);
C3 = U(:,3);
C4 = -1.*U(:,3);

R1 = U*W*V';
if det(R1) < 0
    R1 = -1.*R1;
end
R2 = U*W*V';
if det(R2) < 0
    R2 = -1.*R2;
end
R3 = U*W'*V';
if det(R3) < 0
    R3 = -1.*R3;
end
R4 = U*W'*V';
if det(R4) < 0
    R4 = -1.*R4;
end


Cset = {C1, C2, C3, C4};
Rset = {R1, R2, R3, R4};