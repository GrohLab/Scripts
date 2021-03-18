% PercentageDifferences


percentdif = zeros(length(Counts{1,1}),1);
sp1 = mean(Counts{1,1}')'/(responseWindow(2) - responseWindow(1));
sp2 = mean(Counts{2,1}')'/(responseWindow(2) - responseWindow(1));


for a = 1:length(Counts{1,1})
percentdif(a,1) = ((sp2(a,1) - sp1(a,1))/sp1(a,1))*100;
end

% Filter for FRs significantly altered by CNO  
percentdif(find(Results(1).Activity(1).Pvalues >= 0.05)) = [];
fprintf([' no. Total = ', num2str(numel(percentdif)), '  \n']);
% Inf is 'silent' cell before CNO that becomes active
fprintf([' no. Inf = ', num2str(numel(percentdif(percentdif == inf))), '  \n']);
fprintf([' no. Decreased FR = ', num2str(sum(percentdif < 0)), '  \n']);
fprintf([' no. Increased FR = ', num2str(sum(percentdif > 0)), '  \n']);
percentdif(percentdif == inf) = [];
figure; bar(percentdif);
