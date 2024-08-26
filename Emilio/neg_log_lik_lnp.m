function nllik = neg_log_lik_lnp(theta, X, y)
lambda = exp( X * theta );
log_lik = y' * log( lambda ) - sum( lambda );
nllik = -log_lik;
end