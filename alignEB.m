function [mG, stdG, mR, stdR ] = alignEB(dat1, dat2, n1, n2)

%datG = flipud(dat.GROIaveMax); %we're numbering from 1 at L9 to 18 at R9 on Tanya's diagram
%datR = flipud(dat.RROIaveMax); %this means CCW rotation (+ve vrot) gives higher numbers



Gal = [];
Ral = [];

s = size(dat1);

for i = 1:s(2);

    G = dat1(:,i);
    R = dat2(:,i);
    newG = G;
    newR = R;
    
    [val idx] = max(G);
    shift = 9-idx;
    
    for j = 1:16
        k = j+shift;
        
        if k > 16
            newG(k-16) = G(j);
            newR(k-16) = R(j);      
        elseif k < 1
            newG(k+16) = G(j);
            newR(k+16) = R(j);
        else    
            newG(k) = G(j);
            newR(k) = R(j);
        end    
    end
    
    Gal(:,i) = newG;
    Ral(:,i) = newR;
end


mG = [];
stdG = [];
mR = [];
stdR = [];


for i = 1:16
    mG = [mG mean(Gal(i,:))];
    stdG = [stdG std(Gal(i,:))];
    mR = [mR mean(Ral(i,:))];
    stdR = [stdR std(Ral(i,:))];
end