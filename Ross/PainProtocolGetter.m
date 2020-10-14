% PainProtocolGetter
expt_part = 1;
a = expt_part;
Trig(1).name = head69.title;
Trig(2).name = head71.title;
Trig(3).name = head66.title;
Trig(4).name = head67.title;
Trig(5).name = head68.title;
Trig(1).info{a} = chan69;
Trig(1).offset(a) = length(chan69);
Trig(2).info{a} = chan71;
Trig(2).offset(a) = length(chan71);
Trig(3).info{a} = chan66;
Trig(3).offset(a) = length(chan66);
Trig(4).info{a} = chan67;
Trig(4).offset(a) = length(chan67);
Trig(5).info{a} = chan68;
Trig(5).offset(a) = length(chan68);
%% Getting Trigger data
for b = 1:2
    for a  = 1:length(Trig(b).info)
        Obj = StepWaveform(Trig(b).info{a}, fs);
        Trig(b).triggers{a} = Obj.subTriggers;
    end
end
%% Concatenating Triggers
for b = 1:2
    Trig(b).AllTriggers = Trig(b).triggers{1};
    lngth = length(Trig(b).triggers);
    offset = Trig(b).offset(1);
    for a = 2:lngth
        offInd = a;
        Trig(b).AllTriggers = [Trig(b).AllTriggers; offset + Trig(b).triggers{a}];
        offset = offset + Trig(b).offset(offInd);
    end
end
%% Creating Trigger Struct Data
lngth = length(Trig(1).info);
Triggers.Whisker = [Trig(1).info{1}];
Triggers.Laser = [Trig(2).info{1}];
Triggers.CED = [Trig(3).info{1}];
Triggers.LFP = [Trig(4).info{1}];
Triggers.Resp = [Trig(5).info{1}];
for a = 2: lngth
    Triggers.Whisker = [Triggers.Whisker; Trig(1).info{a}];
    Triggers.Laser = [Triggers.Laser; Trig(2).info{a}];
    Triggers.CED = [Triggers.CED; Trig(3).info{a}];
    Triggers.LFP = [Triggers.LFP; Trig(4).info{a}];
    Triggers.Resp = [Triggers.Resp; Trig(5).info{a}];
end
%% Separating Conditions

Conditions(1).name = 'MechALL'; Conditions(1).Triggers = Trig(1).AllTriggers;
Conditions(2).name = 'LaserAll'; Conditions(2).Triggers = Trig(2).AllTriggers;
save(fullfile(dataDir,[expName,'analysis.mat']),'Conditions','Triggers','-v7.3');
