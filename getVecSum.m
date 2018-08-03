function [direc, mag] = getVecSum(list)
%given a vector of intensities, assumes these to be circularly arranged and
%returns the direction and magnitude of the vector sum of these evenly
%spaced constituent vectors

x = 0.0;
y = 0.0;

L = length(list);

for i = 1:L; %divide into L segments
    
    xi = list(i) * cos( (2*pi / L) * i); %x coordinate of vector number i
    yi = list(i) * sin( (2*pi / L) * i);
    
    %display([i, list(i), xi, yi]);
    
    x = x + xi;
    y = y + yi;
end
    
mag = sqrt(x^2 + y^2);

if y > 0 %let the directionrun from 0 to 2pi
    direc = acos( x / mag );
else
    direc = 2*pi - acos( x / mag );
end 