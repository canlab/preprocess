function [L] = ilogit(t)
L = exp(t)./(1+exp(t));