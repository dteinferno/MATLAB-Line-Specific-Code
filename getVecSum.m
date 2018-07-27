function [direc, mag] = getVecSum(list);

x = 0.0;
y = 0.0;

L = length(list);

for i = 1:L;
    
    xi = list(i) * cos( (2*pi / L) * i);
    yi = list(i) * sin( (2*pi / L) * i);
    
    %display([i, list(i), xi, yi]);
    
    x = x + xi;;
    y = y + yi; ;
end
    
mag = sqrt(x^2 + y^2);

if y > 0
    direc = acos( x / mag );
else
    direc = 2*pi - acos( x / mag );
end 