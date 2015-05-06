function [C, R, X0] = DisambiguateCameraPose(Cset, Rset, Xset)

numSat = zeros(1,4);

for i=1:4
    Cpot = Cset{i};
    Rpot = Rset{i};
    
    Xpts = Xset{i};
    
    for j=1:length(Xpts)
        X = Xpts(j,:);
        if (Rpot(3,:)*(X - Cpot')' >= 0)
            numSat(i) = numSat(i) + 1;
        end
    end
end

[~,ind] = max(numSat);

C = Cset{ind};
R = Rset{ind};
X0 = Xset{ind};