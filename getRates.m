function rates = getRates(vecs, ts, circumference)
%Given a list of PVA directions and timepoints (length N), returns a list of the rate of
%change of the directions at all times (length N-1). Also requires a
%'diameter' argument as the distance between timepoints is always assumed
%to be the shortest distance around the circle.

rates = zeros(1, length(vecs)-1) ;

for i = 1 : length(vecs)-1;
    offset = getOffset(vecs(i+1), vecs(i), circumference);
    dt = ts(i+1)-ts(i);
    rates(i) = offset/dt;
    
end

rates = transpose(rates); %return column vector