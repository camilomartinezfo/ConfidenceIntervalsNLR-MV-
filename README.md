# Confidence Intervals for a nonlinear regression model 


Confidence intervals for the average response of a nonlinear regression model.

This repository has three methods for calculate confidence intervals: Delta method, Bates & Watts (1988, p. 59) and a Monte Carlo simulation (Spiess, 2018). 


### Monte Carlo Simulation

It consists in error propagation of the predictor variables. This approach uses the variance and covariance matrix (vcov) of the parameters fitted of the model and predictor variable $\mu_{i}$ where for each $i$, n samples are created from a multivariate normal distribution using as input the vcov matrix (Spiess, 2013).  Hence, the assumption is that for each i there is a normal distribution of the dependent variable. After, with each sample simulated, the function of response variable is calculated. Finally, some statistics are calculated like mean, standard deviation, and quantiles of confidence intervals (Spiess, 2013).  

