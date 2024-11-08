tbl = [68, 33, 217;...
       33, 35, 180];
sumCol = sum(tbl)';
sumFil = sum(tbl,2)';
multSum = (sumCol * sumFil)';
esperado = multSum / sum(tbl(:));
resChi = ((tbl - esperado).^2)./esperado;
resChi = sum(resChi(:));
figure; plot(0:100, chi2pdf(0:100, 2))
hold on; scatter(resChi, chi2pdf(resChi,2))