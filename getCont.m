function vec = getCont(vec, thresh);
%if two consecutive points are more than tresh apart, set the first one to
%nan

for i = 2:length(vec);
   if abs(vec(i)-vec(i-1)) > thresh;
       vec(i-1) = nan;
   end
end