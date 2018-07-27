function [magsL, magsR, dirsL, dirsR] = getPBVec(array);
%%Return lists of magnitudes and directions

magsL = [];
magsR = [];
dirsL = [];
dirsR = [];

s = size(array);

for i = 1:s(2);
    
   vec = array(:,i); 
   vec(1:9);
    
   [dirl, magl] = getVecSum( vec(1:9) );
   
   [dirr, magr] = getVecSum( vec(10:18) );
      
   dirl = 9*dirl / (2*pi);
   dirr = 9 + 9*dirr / (2*pi);
   
   magsL = [magsL magl];
   magsR = [magsR magr];
   
   dirsL = [dirsL dirl];
   dirsR = [dirsR dirr];
   
end

magsL = transpose(magsL)
magsR = transpose(magsR)
dirsL = transpose(dirsL)
dirsR = transpose(dirsR)
    
    
end