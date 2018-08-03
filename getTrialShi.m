function [DF, data] = getTrialShi(cond, period, i, k, smoothF, smoothV)
%Given a container with data for a given experiment, as well as a period, flyID and
%trialID of interest, returns DF with the raw imaging data (possibly
%smoothened) as well as a cell 'data' containing tPts, vRot, vF, intensity,
%|PVA|, PVA direction for RT and 30C (in that order; this is a pretty silly
%data structure and I should probably change it at some point)

%i is flyID, k is trialID (1, 2 for RT AND 30C)
%returns data where 1-6 is RT, 7-12 30C
%1 is tPts, 2 is vRot, 3 is vF, 4 is int, 5 is |PVA|, 6 is direction


data = { [] [] [] [] [] [] [] [] [] [] [] [] };


conds = {'All' 'All_30C'};
    
for j = 1:2; %iterate over conditions

    if j == 1
        ind = k;
    else
        ind = k+1; %don't use first trial from 30C. Use 2-3 instead of 1-2
    end

    [tPts,darkPer, OLPer, CLPer, CWPer, CCWPer, DF, heading, headingPlt, vRot, vF, stripePos, stripePosPlt, stripeJumps]...
        = extractShiData(cond, i, conds{j}, ind, smoothF, smoothV); %get data

    int = mean(DF,1); %use mean rather than sum since this actually gives a dF/F like value

    %getmagnitude of PVAs
    m = zeros(1, length(int));
    dirs = zeros(1, length(int));
    for ind = 1:length(int);
        vec = DF(:, ind);

        [di, mi] = getVecSum( vec );
        m(ind) = mi/sum(abs(vec)); %use normalized PVA amplitude
        dirs(ind) = di;
    end

    vF = transpose(vF); %use row vectors
    vRot = transpose(vRot);
    tPts = transpose(tPts);

    %pick out indices corresponding to our period of interest
    if strcmp(period, 'dark')
        per = darkPer;
    elseif strcmp(period, 'OL')
        per = OLPer;
    elseif strcmp(period, 'CL')
        per = CLPer;
    else
        display('per not recognized, please select dark or OL')
    end

    %this is pretty silly, sorry to anyone who has to deal with it
    data{6*j-5} = [data{6*j-5} tPts(per)];
    data{6*j-4} = [data{6*j-4} vRot(per)];%vRot
    data{6*j-3} = [data{6*j-3} vF(per)];%forward vel
    data{6*j-2} = [data{6*j-2} int(per)];%sum of intensity
    data{6*j-1} = [data{6*j-1} m(per)];%magnitude of PVA
    data{6*j} = [data{6*j} dirs(per)];
    
    DF = DF(:, per);

end