function Conditions = SeparateLaserFrequencies(Trig, fs)
offset = 0;
c = 2;
promptStrings = {'no. of Intensities','no. of Stimulation Patterns'};
defInputs = {'1', '3'};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
if isempty(answ)
    fprintf(1,'Cancelling...\n')
    return
else
    n = str2num(answ{1,1});
    m = str2num(answ{2,1});
end

for a = 1:n
    promptStrings = {'Laser Intensity: mW'};
    defInputs = {'1'};
    answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
    if isempty(answ)
        fprintf(1,'Cancelling...\n')
        return
    else
        intensityName = [num2str(answ{1}), 'mW'];
    end
    for b = 1:m
        laser = Trig(2).triggers{1,a};
        if b == 1
            laserStim = [false, false; diff(laser) > 0.15*fs & diff(laser) < 1.1*fs];
            patternName = '1Hz';
        elseif b == 2
            laserStim = [false, false; diff(laser) < 0.15*fs];
            patternName = '10Hz';
        end
        
        On = diff(laserStim(:,1)) == 1;
        On = [On; false];
        Off = diff(laserStim(:,1)) == -1;
        Off = [Off; false];
        if b == 3
            laserStim = diff(laser') > 3*fs;
            On = laserStim;
            Off = On;
            patternName = 'Continuous_Pulse';
        end
        Conditions(c).name = ['Laser_', patternName, '_', intensityName];
        Conditions(c).Triggers = laser(On,1) + offset;
        Conditions(c).Triggers(:,2) = laser(Off,2) + offset;
        c = c + 1;
    end
    offset = offset + Trig(2).offset(a);
end
Conditions(1).name = 'Laser_All';
Conditions(1).Triggers = Conditions(2).Triggers;
for a = 3: length(Conditions)
    Conditions(1).Triggers = [Conditions(1).Triggers; Conditions(a).Triggers];
end
Conditions(1).Triggers = sort(Conditions(1).Triggers, 'ascend');
end

