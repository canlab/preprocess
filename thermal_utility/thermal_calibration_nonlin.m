function [p,errval,fit,linkfun,fhan,target_temps, ramp_rate] = thermal_calibration_nonlin(ratings,temps)
%  thermal_calibration(ratings,temps)
%  1. choose temp, 2. calibrate 8 places, choose 3 most similar; test 3
%  temps on those three spots 
%  "temps"= [t1;t2;t3] 
%  "ratings" = [avg_t1;avg_t2;avg_t3]
%  choose target ratings
% output will be target temps to match those ratings, and ramp rate to set
% for 1.5s increase from baseline temp of 32.

target_ratings = [1;3;6;8];


%  P = polyfit(ratings,temps,1); %that was okay for linear.
%  int = P(2); slope = P(1);
  
%  target_temps = slope * target_ratings + int;

linkfun = @(p, x) pain_predict_function(p, x);
[p,errval,fit,linkfun,fhan] = nonlin_fit(ratings,temps, 'link', linkfun, 'start', [3 1 2], 'plot', 'verbose');

target_temps = (target_ratings / (p(2))) .^ (1/(p(3))) + 32 + p(1);
  
ramp_rate = (target_temps - 32) / 1.5;
  
  
hold on;
plot(target_temps, target_ratings, 'ko', 'MarkerFaceColor', 'r');

xlabel('Temperature'); ylabel('Rating');

end

