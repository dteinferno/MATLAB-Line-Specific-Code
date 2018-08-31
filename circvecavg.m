function [mag,direc] = circvecavg(vector)
%given a list of unit vectors, calculates  the mean direction and magnitude
%of their average


x = 0.0;
y = 0.0;

L = length(vector);

for i = 1:L; %divide into L segments
    
    xi = cos(vector(i)); %x coordinate of vector number i
    yi = sin(vector(i));
    
    %display([i, list(i), xi, yi]);
    
    x = x + xi;
    y = y + yi;
end
    
mag = sqrt(x^2 + y^2); %normalize magnitude

if y > 0 %let the directionrun from -pi to pi
    direc = acos( x / mag );
else
    direc = - acos( x / mag );
end

mag = mag / L;