for a = length(Results)
fprintf(1, [Results(a).name,' Kolgorov-Smirnov Test = ', num2str(Results(a).kstest), ' \n \n'])
end