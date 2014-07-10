function [mu,sigma] = logninvstat(m,v)

%notation from the Matlab documentation for the log-normal pdf. given the 
%mean (m) and variance (v) or a lognormal distribution, calculate
%the mean (mu) and standard deviation (sigma) of the associated normal
%distribution. so if x is log-normally distributed, m is the mean of x, and
%v is the variance of x, and mu is the mean of log(x) and sigma is the
%standard deviation of log(x).

mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

return