function maxarray = maxk(array, k, dimension)


minus_base = false;

if minus_base
    array = array - repmat( min(array, [], dimension), [16, 1]);
end

%returns k highest values along dimension (1 or 2)

array = sort(array, dimension, 'descend');

if dimension == 1
    maxarray = array(1:k, :);
elseif dimension == 2
    maxarray = array(:, 1:k);
else
    display('dimension must be 1 or 2')
end
