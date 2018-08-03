function [magsL, magsR, dirsL, dirsR] = getPBVec(array);
%%Return lists of magnitudes and directions

%initialize
magsL = zeros(length(array), 1);
magsR = zeros(length(array), 1);
dirsL = zeros(length(array), 1);
dirsR = zeros(length(array), 1);

s = size(array);
for i = 1:s(2); %go through time points
    
   vec = array(:,i); 
    
   [dirl, magl] = getVecSum( vec(1:9) ); %split into left and right and get magnitude and direction of vector sum
   
   [dirr, magr] = getVecSum( vec(10:18) );
      
   dirl = 9*dirl / (2*pi); %convert to radians
   dirr = 9 + 9*dirr / (2*pi);
   
   magsL(i) = magl; %add to vectors
   magsR(i) = magr;
   
   dirsL(i) = dirl;
   dirsR(i) = dirr;
   
end