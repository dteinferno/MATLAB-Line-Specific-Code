%go through directories with EB data and run analyses

dirs = {'~/Documents/Imaging/Data_Dan/EB/PEN2_G_PEG_R_EB/',...
    '~/Documents/Imaging/Data_Dan/EB/PEN1_G_EPG_R_EB/',...
    '~/Documents/Imaging/Data_Dan/EB/PEN1_R_EPG_G_EB/',...
    '~/Documents/Imaging/Data_Dan/EB/EPG_G_GE_R_EB/',...
    '~/Documents/Imaging/Data_Dan/EB/EPG_R_GE_G_EB/',...
    '~/Documents/Imaging/Data_Dan/EB/PEN2_G_EPG_R_EB/',...
    '~/Documents/Imaging/Data_Dan/EB/PEN2_R_EPG_G_EB/'};

greens = {'PEN2', 'PEN1', 'EPG', 'EPG', 'GE', 'PEN2', 'EPG'};

reds = {'PEG', 'EPG', 'PEN1', 'GE', 'EPG', 'EPG', 'PEN2'};

for i = 1:length(dirs)
    
    %display(dirs{i})
    
    try
        from_file = load(strcat(dirs{i}, 'cont'), 'alldata');
        cond=from_file.alldata;    
    catch
        cond = FlyDatLoad(1);
        save(strcat(dirs{i}, 'cont'), 'alldata');
    end
    
    getDecayDelayEBPB(strcat(dirs{i}, 'decayDelay/'), cond, 7, 10, greens{i})
    %EB_analysis(dirs{i}, greens{i}, reds{i})    
    %scatterVelEB(dirs{i}, greens{i}, reds{i})
    %alignbyVelEB(dirs{i}, greens{i}, reds{i})

    
end