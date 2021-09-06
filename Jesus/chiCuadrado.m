tbl = [68, 250; 33, 218];
sumCol = sum(tbl)';
sumFil = sum(tbl,2)';
multSum = (sumCol * sumFil)';
esperado = multSum / sum(tbl, "all");
resChi = sum(((tbl - esperado).^2)./esperado, "all");
figure; plot(0:100, chi2pdf(0:100, 1))
hold on; scatter(resChi, chi2pdf(resChi,1))