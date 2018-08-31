function meanvecmag = meanvecmag(vector, window)

vecmags = zeros(1,length(vector)-window);

for i = 1:length(vector)-(window-1)
    
    [mag, direc] = circvecavg(vector(i:i+window-1)); %get the vector average of this set of orientations
    
    vecmags(i) = mag;
    
end

meanvecmag = mean(vecmags);