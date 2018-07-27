function newvec = Smooth(vec, range)
%make range an odd number; the number of cells we average over

newvec = [];

for i = (range+1)/2:length(vec)-(range-1)/2
   m = mean( [vec(i-1), vec(i), vec(i+1)] );
   newvec = vertcat(newvec, m);
    
end