% getDirection
% given you've run the statstest function to give Counts and Results
% we can get say  the no. of clusters that increase or decease activity
% for conditon effect a = 1, for MR a = 2


%%   Significant mechanical responses per condition
b = Nccond;
for a = 1 : Nccond
   
    
    
    
    fprintf([ consCondNames{1,a},' nMR', ' = ',  num2str(sum(Results(b).Activity(1).Pvalues < 0.05)), '  \n']);
    dC = mean(Counts{a,2}') - mean(Counts{a,1}');
    fprintf([ consCondNames{1,a},' Increased MR', ' = ',   num2str(sum((Results(b).Activity(1).Pvalues < 0.05) & dC' > 0)), '  \n']);
    fprintf([ consCondNames{1,a},' Decreased MR', ' = ',  num2str(sum((Results(b).Activity(1).Pvalues < 0.05) & dC' < 0)), '  \n']);
    
    b = b + Nccond - a;
    
end
%% Comparing Across Conditions (i.e. CNO effect on unevoked and evoked activity)
     sZ = size(Counts);
     d = Nccond;
     e = 1;
     f = Nccond;
for a = 1:Nccond - 1
    for b = (a + 1): Nccond
        
        for c = 1: sZ(1,2)
        
        dC = mean(Counts{b,c}') - mean(Counts{a,c}');
        fprintf([ consCondNames{1,a},' vs ', consCondNames{1,b}, ' ', Results(e).Activity(c).Type, ' differences = ', num2str(sum(Results(e).Activity(c).Pvalues < 0.05)), '  \n']);
        fprintf([ consCondNames{1,a},' vs ', consCondNames{1,b}, ' ', Results(e).Activity(c).Type, ' increases = ', num2str(sum((Results(e).Activity(c).Pvalues < 0.05) & dC' > 0)), '  \n']);
        fprintf([ consCondNames{1,a},' vs ', consCondNames{1,b}, ' ',Results(e).Activity(c).Type, ' decreases = ', num2str(sum((Results(e).Activity(c).Pvalues < 0.05) & dC' < 0)), '  \n']);
        end
        
        e = e + 1;
        
        if e == f
            e = e + 1;
        end
    end
    d = d - 1;
    f = f + d;
end
