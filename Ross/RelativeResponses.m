% Relative responses

b = Nccond;
for a = 1 : Nccond
   
  
     sigDif = Results(b).Activity(1).Pvalues < 0.05;
     ind = find(sigDif);
     RelativeResp(a).name = consCondNames{1,a};
     RelativeResp(a).Resp = mean(Counts{a,2}(ind,:)') - mean(Counts{a,1}(ind,:)');
     b = b + Nccond - a;
end
save(fullfile(dataDir,[expName,'RelativeResponses.mat']),'RelativeResponses','-v7.3');
