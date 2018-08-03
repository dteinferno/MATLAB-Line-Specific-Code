function [data, vRots] = align1C(cond, period)
%Given a container with data for a given experiment and a period of interest, realigns to maximum
%intensity at 5 for each timeframe and pools over flies and trials. Returns
%'data' where data{1} is RT, data{2} is 30C, and vRots with same structure.
%Changing smoothing requires hardcording at the moment.

conds = {'All' 'All_30C'};

%let's not smooth since we're averaging over time anyways and I guess we shouldn't
%mess with the data if we don't have to
smoothV = 0;
smoothF = 0;

data = {[], []}; %data{1} is RT data{2} is 30C. Pool across flies and trials
vRots = {[], []}; %same
    
for j = 1:2 %iterate over conditions
    
    Gals = []; %initialize data
    
    for i = 1:length(cond{1}.allFlyData); %iterate over flies
   
        for k = 1:2 %iterate over trials
            
            if j == 1
                ind = k;
            else
                ind = k+1; %don't use first trial from 30C; only use 2-3 (given by k=1-2)
            end
        
            [tPts,darkPer, OLPer, CLPer, CWPer, CCWPer, DF, heading, headingPlt, vRot, vF, stripePos, stripePosPlt, stripeJumps]...
                    = extractShiData(cond, i, conds{j}, ind, smoothF, smoothV);
                
            
            %find period of interest    
            if strcmp(period, 'dark')
                per = darkPer;
            elseif strcmp(period, 'CL')
                per = CLPer;
            elseif strcmp(period, 'OL')
                per = OLPer;
            else
                display('period not recognized')
            end
                
            DF = DF(:,per); %get dF/F data for period
            
            vRots{j} = [vRots{j} vRot(per)']; %add vRot data to the pooled data
            
            s = size(DF);
            
            Gal = zeros(16, s(2));

            for index = 1:s(2); %for each timepoint
                G = DF(:,index); %get vector for that timepoint
                newG = G;
                [~, idx] = max(G);
                shift = 9-idx;

                for jnew = 1:16 %realign so maximum is at 9
                    knew = jnew+shift;

                    if knew > 16
                        newG(knew-16) = G(jnew);     
                    elseif knew < 1
                        newG(knew+16) = G(jnew);
                    else    
                        newG(knew) = G(jnew);
                    end    
                end

                Gal(:,index) = newG; %add new data to Gal
            end
            Gals = [Gals Gal]; %add re-aligned data to pooled re-aligned data

        end
    end
    data{j} = Gals; %add to data. Data now contains full re-aligned dataset
end
