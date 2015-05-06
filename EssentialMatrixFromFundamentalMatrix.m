function E = EssentialMatrixFromFundamentalMatrix(F, K)

E = K'*F*K;

[U,~,V] = svd(E);

Dcorr = [1 0 0;
         0 1 0;
         0 0 0];
     
E = U*Dcorr*V';