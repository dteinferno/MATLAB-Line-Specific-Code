function [mG, stdG, mR, stdR ] = alignEB(dat1, dat2)
%given two single trials for single flies, aligns trial 1 and trial 2 to the
%max of trial 1 (at index 9) and returns the mean and standard deviation at
%each position



s = size(dat1);

Gal = zeros(s(1),s(2));
Ral = zeros(s(1),s(2));

for i = 1:s(2); %go over each timepoint

    G = dat1(:,i);%get the vectors
    R = dat2(:,i);
    newG = G;
    newR = R;
    
    [~, idx] = max(G);
    shift = 9-idx;%this is how much we need to shift stuff
    
    for j = 1:16 %shift each element of the green and red vectors
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
    
    Gal(:,i) = newG; %update data
    Ral(:,i) = newR;
end


mG = mean(Gal, 2)'; %mean green intensity
stdG = zeros(1,16);
mR = mean(Ral, 2)';
stdR = zeros(1,16);


for i = 1:16
    stdG(i) = std(Gal(i,:)); %get std at each position
    stdR(i) = std(Ral(i,:));
end