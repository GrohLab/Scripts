for a = 1:length(ISIhist)
    IsiStruct(a).name = ISIhist(a).name;
    IsiStruct(a).Stacks(1).name = 'Baseline';
    IsiStruct(a).Stacks(1).ISIs = ISIhist(a).Vals(1).TriggeredIsI;
    IsiStruct(a).Stacks(2).name = 'Evoked';
    IsiStruct(a).Stacks(2).ISIs = ISIhist(a).Vals(2).TriggeredIsI;
   
    % Unit ReturnMap loop (Baseline)
    clear new;
    for cu = 1:length(IsiStruct(a).Stacks(1).ISIs(:,1,1))
        firstTrialISIs = IsiStruct(a).Stacks(1).ISIs(cu,:,1);
        ftInd = firstTrialISIs ~= 0;
        firstTrialISIs = firstTrialISIs(ftInd);
        current = firstTrialISIs(1:end-1);
        next = firstTrialISIs(2:end);
        isiReturnMap = [current', next'];
        
        
        for ct = 2:length(IsiStruct(a).Stacks(1).ISIs(cu,1,:))
            TrialISIs = IsiStruct(a).Stacks(1).ISIs(cu,:,ct);
            tInd = TrialISIs ~= 0;
            TrialISIs = TrialISIs(tInd);
            current = TrialISIs(1:end-1);
            next = TrialISIs(2:end);
            new = [current', next'];
            isiReturnMap = [isiReturnMap; new];
        end
        IsiStruct(a).Stacks(1).ReturnMap(cu).Unit = isiReturnMap;
    end
    
    % Condition ReturnMap loop (Baseline)
     IsiStruct(a).Stacks(1).ConditionMap =  IsiStruct(a).Stacks(1).ReturnMap(1).Unit;
     for cu = 2:length(IsiStruct(a).Stacks(1).ReturnMap)
         IsiStruct(a).Stacks(1).ConditionMap = [IsiStruct(a).Stacks(1).ConditionMap; IsiStruct(a).Stacks(1).ReturnMap(cu).Unit];
     end
    
     
     
     
     % Unit ReturnMap loop (Evoked)
     clear new;
     for cu = 1:length(IsiStruct(a).Stacks(2).ISIs(:,1,1))
        firstTrialISIs = IsiStruct(a).Stacks(2).ISIs(cu,:,1);
        ftInd = firstTrialISIs ~= 0;
        firstTrialISIs = firstTrialISIs(ftInd);
        current = firstTrialISIs(1:end-1);
        next = firstTrialISIs(2:end);
        isiReturnMap = [current', next'];
        
        
        for ct = 2:length(IsiStruct(a).Stacks(2).ISIs(cu,1,:))
            TrialISIs = IsiStruct(a).Stacks(2).ISIs(cu,:,ct);
            tInd = TrialISIs ~= 0;
            TrialISIs = TrialISIs(tInd);
            current = TrialISIs(1:end-1);
            next = TrialISIs(2:end);
            new = [current', next'];
            isiReturnMap = [isiReturnMap; new];
        end
        IsiStruct(a).Stacks(2).ReturnMap(cu).Unit = isiReturnMap;
    end
    
    % Condition ReturnMap loop (Evoked)
     IsiStruct(a).Stacks(2).ConditionMap =  IsiStruct(a).Stacks(2).ReturnMap(1).Unit;
     for cu = 2:length(IsiStruct(a).Stacks(2).ReturnMap)
         IsiStruct(a).Stacks(2).ConditionMap = [IsiStruct(a).Stacks(2).ConditionMap;IsiStruct(a).Stacks(2).ReturnMap(cu).Unit];
     end
        
end
save(fullfile(dataDir,[expName,'_ISI_Struct.mat']), 'IsiStruct', '-v7.3');
