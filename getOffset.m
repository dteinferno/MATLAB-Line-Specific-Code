function offset = getOffset(vec1s, vec2s, diameter)
%%given a list of directions of vectors at timepoints, finds the offset in
%%radians


offset = [];

for i = 1:length(vec1s);
    a = vec1s(i);
    b = vec2s(i);
    
    dif = abs(a - b);
    
    if dif <= diameter/2
        offset = [ offset a-b ];
        
    else
        if b > a
            offset = [ offset a+diameter-b ];
            
        else
            offset = [ offset a-(b+diameter) ];
    
        end
    end
    
end

offset = transpose(offset);

