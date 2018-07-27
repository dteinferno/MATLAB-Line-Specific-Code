dirs = { '~/Documents/Imaging/Data_Dan/PEN1_G_EPG_R_PB/',...
    '~/Documents/Imaging/Data_Dan/PEN1_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/PEN2_G_PEG_R_PB/',...
    '~/Documents/Imaging/Data_Dan/PEN2_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/D7_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/D7_G_EPG_R_PB/' };

greens = { 'PEN1' 'EPG' 'PEN2' 'EPG' 'EPG' 'D7' };

reds = { 'EPG' 'PEN1' 'PEG' 'PEN2' 'D7' 'EPG' };

for i = 1:length(dirs)
    display(dirs{i})
    PB_analysis(dirs{i}, greens{i}, reds{i})
    scatterVelPB(dirs{i}, greens{i}, reds{i})
    alignByVelPB(dirs{i}, greens{i}, reds{i})
    
end
