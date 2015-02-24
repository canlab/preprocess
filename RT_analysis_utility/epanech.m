function k = epanech(t)

% EPANECH Epanechnikov kernel estimate of the hazard function
% from Trisha van Zandt

% The Epanechnikov kernel takes arguments between -sqrt(5) and
% sqrt(5). All other values return 0.

% 3 parts to routine: functions epanech and Iepanech and the function
% 'hazard', which returns the hazard function estimate
% input argument (t) is the vector of points at which the estimate is
% desired

k = (abs(t)<sqrt(5)).*(.75.*(1.-0.2.*t.^2)/sqrt(5));