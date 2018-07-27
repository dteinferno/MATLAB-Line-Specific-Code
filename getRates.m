function rates = getRates(vecs, ts, diameter);

rates = [];

for i = 1 : length(vecs)-1;
    offset = getOffset(vecs(i+1), vecs(i), diameter);
    dt = ts(i+1)-ts(i);
    rates = [ rates offset/dt ];
    
end

rates = transpose(rates);