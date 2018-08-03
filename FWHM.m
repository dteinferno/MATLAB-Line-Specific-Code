function [xs, HM] = FWHM(means, baseline)
%Given a vector (means), returns the x values of the full-width at half
%maximum as well as the half-maximum value. The baseline argument tells
%whether to get FWHM from 0 or from the minimum value of the data set
%('min')

[maxval, idx] = max(means);
minval = min(means);

if strcmp(baseline, 'min')
    HM = 0.5*(maxval + minval); %find FWHM relative to baseline
else
    %actually et reative to zero since wider bump gives higher acitvity on
    %opposite side of EB
    HM = 0.5*maxval;
end

n = 100; %resolution of the script

means = {means(1:idx), means(idx:end)}; %divide into left and right of maximum

xs = {[] []};

for i = 1:2; %find where y = HM in left and right half individually

    points = zeros( 1, (length(means{i})-1)*n + 1 );

    for j = 1:length(means{i})-1

        points(n*(j-1)+1:n*j+1) = linspace(means{i}(j),means{i}(j+1),n+1); %construct grid

    end

    difs = abs(points - HM); %get deviation from HM
    
    [~, HMidx] = min(difs); %this is the HM point
    
     %convert back to the indexing of the input vector
    try
        if i == 1
            x = ( HMidx - 1 ) / n + 1; %our index starts at 1
        else
            x = ( HMidx - 1 ) / n + length(means{1}); %add the other half
        end
    catch %if baseline is higher than HM, we just return the min and max x values
        if i == 1
            x = 1
        else
            x = length(means{1})+length(means{2}) %if we don't reach the HM value
        end
    end
        
    xs{i} = x; %add to list of xs
    
end