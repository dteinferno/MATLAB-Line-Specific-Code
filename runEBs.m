
dirs = {'~/Documents/Imaging/Data_Dan/PEN2_G_PEG_R_EB/',...
    '~/Documents/Imaging/Data_Dan/PEN1_G_EPG_R_EB/',...
    '~/Documents/Imaging/Data_Dan/PEN1_R_EPG_G_EB/',...
    '~/Documents/Imaging/Data_Dan/EPG_G_GE_R_EB/',...
    '~/Documents/Imaging/Data_Dan/EPG_R_GE_G_EB/',...
    '~/Documents/Imaging/Data_Dan/PEN2_G_EPG_R_EB/',...
    '~/Documents/Imaging/Data_Dan/PEN2_R_EPG_G_EB/'}

greens = {'PEN2', 'PEN1', 'EPG', 'EPG', 'GE', 'PEN2', 'EPG'}

reds = {'PEG', 'EPG', 'PEN1', 'GE', 'EPG', 'EPG', 'PEN2'}

for i = 1:length(dirs)
    display(dirs{i})
    EB_analysis(dirs{i}, greens{i}, reds{i})    
    scatterVelEB(dirs{i}, greens{i}, reds{i})
    alignbyVelEB(dirs{i}, greens{i}, reds{i})

    
end