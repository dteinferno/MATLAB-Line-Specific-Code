clear


dirs = { '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/EPG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN1/'};


for i = 1:length(dirs)-1;
    dir = dirs{i}
    
    try
        from_file = load(strcat(dir, 'cond'), 'cond');
        cond=from_file.cond;    
    catch
        cond = FlyDatLoad(1);
        save(strcat(dir, 'cond'), 'cond');
    end
    
    data = scatterVelEB1C(strcat(dir, 'intensities/'), cond, 'OL');
    
    %cor = correlationVel(strcat(dir, 'correlation/'), cond);
    
    %[decays, delays] = getDecayDelay(strcat(dir, 'decayDelay/'), cond, 3, 10);
    %dir, cond, length of stretch that must be zero, length of stretch to
    %consider for re-emergence of fluorescence
    
    %alignByVel1C(strcat(dir, 'aligned/'), cond)
    
end
