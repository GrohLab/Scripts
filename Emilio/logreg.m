function y_hat = logreg(Xtrain, ytrain, Xtest)
mdl = lassoglm( Xtrain, ytrain, 'Distribution', 'binomial', 'Link', 'logit' );
y_hat = predict(mdl, Xtest);
end