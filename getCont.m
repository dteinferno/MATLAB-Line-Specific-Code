function vec = getCont(vec, thresh);

for i = 1:length(vec)-1;
   if abs(vec(i+1)-vec(i)) > thresh;
       vec(i+1) = nan;
   end
end