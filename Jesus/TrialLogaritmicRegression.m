% Create random data
x = y3dsampleTodos(:,1);
y = y3dsampleTodos(:,2);
%y(1:50) = x(1:50) > 0.3; % To avoid perfect seperation
% Fit modelx,
mdl = fitglm(x, y,'Distribution','binomial');
xnew = (linspace(0,1281,1281)'); % test data
ynew = predict(mdl, xnew);
scatter(x, y);
hold on;
%plot(xnew(:,2), ynew);
plot(xnew,ynew)
%%
% Create random data
x = rand(100, 1);
y = x > 0.5;
y(1:50) = x(1:50) > 0.3; % To avoid perfect seperation
% Fit model
mdl = fitglm(x, y,"quadratic", "Distribution", "binomial");
xnew = linspace(0,1,1000)'; % test data
ynew = predict(mdl, xnew);
scatter(x, y);
hold on;
plot(xnew, ynew);
%%
x = xycombined;
y = xy3DTodosFin(:,3);
%y(1:50) = x(1:50) > 0.3; % To avoid perfect seperation
% Fit modelx,
mdl = fitglm(x, y, "Distribution", "binomial");
xnew = linspace(0,50000,3000)'; % test data
ynew = predict(mdl, xnew);
scatter(x, y);
hold on;
plot(xnew, ynew);