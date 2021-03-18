function findCCs(corrs)
N = size(corrs{1},2);
xh = (1:(N-1)/2)/fs;
x = [-flip(xh), 0 , xh];


for a  = randi(length(corrs),1,25)
figure; plot(x, corrs{a}(2,:))
hold on
[sp, mdls] = fitSpline(x, corrs{a}(2,:), 1, 5e-2, 2/3);
plot(x, sp)
end