% Transform Conditions from the ProtocolGetter to LTP 
NTa = size(Conditions(1).Triggers,1);
[biGaps, gapSubs] = sort(diff(Conditions(1).Triggers(:,1)), 'descend');
gapFlag = abs(zscore(biGaps)) > 3;
mainGaps = sort(gapSubs(gapFlag), 'ascend');
Ng = numel(mainGaps);
Ncond = Ng+1;
condFlags = false(NTa, Ncond);
initSub = 1;
for cg = 1:Ng
    condFlags(initSub:mainGaps(cg),cg) = true;
    initSub = mainGaps(cg) + 1;
end
condFlags(initSub:NTa, Ncond) = true;
Na = sum(condFlags,1);
%% Creating the LTP Condition structure
NewConditions = struct('name','Control','Triggers',...
    Conditions(1).Triggers(condFlags(:,1),:));
NewConditions(2) = Conditions(2); NewConditions(2).name = 'Induction';
NewConditions(3) = struct('name','AfterInduction','Triggers',...
    Conditions(1).Triggers(condFlags(:,2),:));
NewConditions(4) = Conditions(1);
NewConditions(5) = Conditions(2);
Conditions = NewConditions;
