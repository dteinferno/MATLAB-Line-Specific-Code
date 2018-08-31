function [DFs, data, flyinds] = getData(cond, period, smoothF, smoothV, nWedge)
%Given a container with data for a given experiment, as well as a period of interest, returns DF with the raw imaging data (possibly
%smoothened) as well as a cell 'data' containing tPts, vRot, vF, intensity,
%|PVA|, PVA direction for RT and 30C (in that order; this is a pretty silly
%data structure and I should probably change it at some point). Pooled
%across flies and trials. Only considers green channel.
%nWedge specifies the number of max intensity wedge to calculate mean
%intensity over

%i is flyID, k is trialID (1, 2 for RT AND 30C)
%returns data where 1-6 is RT, 7-12 30C
%1 is tPts, 2 is vRot, 3 is vF, 4 is int, 5 is |PVA|, 6 is direction



data = { [] [] [] [] [] [] };

DFs = [];

flyinds = [];

for i = 1:length(cond{1}.allFlyData); %iterate over flies
    flyinds = [flyinds length(data{4})+1];
    
    try%we may not have consistent nomenclature
        L = length(cond{1}.allFlyData{i}.Dark);
        label = 'Dark';
    catch
        L = length(cond{1}.allFlyData{i}.All);
        label = 'All';
    end
   
    for k = 1:L %iterate over trials


        [tPts,darkPer, OLPer, CLPer, CWPer, CCWPer, DF, heading, headingPlt, vRot, vF, stripePos, stripePosPlt, stripeJumps]...
            = extractShiData(cond, i, label, k, smoothF, smoothV); 

        
        int = mean( maxk(DF, nWedge, 1), 1); %specify how many wedges to calculate the mean over

        %getmagnitude of PVAs
        m = zeros(1, length(int));
        dirs = zeros(1, length(int));
        for ind = 1:length(int);
            vec = DF(:, ind);

            [di, mi] = getVecSum( vec );
            m(ind) = mi/sum(abs(vec)); %normalized magnitude of PVA
            dirs(ind) = di;
        end

        vDF = transpose( getRates(dirs, tPts, 2*pi) ); %work with row vectors
        vF = transpose(vF);
        vRot = transpose(vRot);
        tPts = transpose(tPts);

        %get indices of period of interest
        if strcmp(period, 'dark')
            per = darkPer;
        elseif strcmp(period, 'OL')
            per = OLPer;
        elseif strcmp(period, 'CL')
            per = CLPer;
        else
            display('per not recognized, please select dark or OL')
        end

        DF = DF(:, per);
        
        
        size(DF);
        
        DFs = [DFs DF];
        try
            data{1} = [data{1} tPts(per)];
            data{2} = [data{2} vRot(per)];
            data{3} = [data{3} vF(per)];
            data{4} = [data{4} int(per)]; %mean fluorescent intensity
            data{5} = [data{5} m(per)]; %PVA magnitudes
            data{6} = [data{6} dirs(per)]; %PVA directions
        catch %reaching end of index
            data{1} = [data{1} tPts(per(1:end-1))];
            data{2} = [data{2} vRot(per(1:end-1))];
            data{3} = [data{3} vF(per(1:end-1))];
            data{4} = [data{4} int(per(1:end-1))]; %mean fluorescent intensity
            data{5} = [data{5} m(per(1:end-1))]; %PVA magnitudes
            data{6} = [data{6} dirs(per(1:end-1))]; %PVA directions
        end
    end

end