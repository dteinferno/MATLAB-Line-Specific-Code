function offset = getOffset(vec1s, vec2s, circumference)
%%given two list of directions of vectors at the same timepoints, returns a
%%vector of the CCW offsets in radians. Also requires a circumference argument as the
%%offset is always considered to be the shortest path along the circle.
%%Offset FROM vec2s TO vec1s (i.e. vec1s-vec2s)


offset = zeros(length(vec1s), 1);

for i = 1:length(vec1s);
    a = vec1s(i);
    b = vec2s(i);
    
    dif = abs(a - b);
    
    if dif <= circumference/2 %this is already shortest path, return simple difference
        offset(i) = a-b;
        
    else %need to go the other way around the circle
        if b > a %positive offset the other way around
            offset(i) = a+circumference-b;
            
        else %negative offset the other way around
            offset(i) = a-(b+circumference);
    
        end
    end
    
end

offset = transpose(offset);

