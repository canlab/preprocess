function [target_temps, ramp_rate] = thermal_calibration(ratings,temps)
%  thermal_calibration(ratings,temps)
%  1. choose temp, 2. calibrate 8 places, choose 3 most similar; test 3
%  temps on those three spots 
%  "temps"= [t1;t2;t3] 
%  "avg_ratings" = [avg_t1;avg_t2;avg_t3]
%  output will be temps that correspond to ratings of 20, 50, and 80 using
%  a linear polynomial fit by least squares, and then the rate to go from
%  32 deg to target temp in 1.5s

figure; plot(temps, ratings, 'ko', 'MarkerFaceColor', [.5 .5 .5]);

target_ratings = [2;5;8];


  P = polyfit(ratings,temps,1);
  int = P(2); slope = P(1);
  
  target_temps = slope * target_ratings + int;
  
  ramp_rate = (target_temps - 32) / 1.5;
  
  
  figure; plot(temps, ratings, 'ko', 'MarkerFaceColor', [.5 .5 .5]);
refline
hold on;
plot(target_temps, target_ratings, 'ko', 'MarkerFaceColor', 'r');

xlabel('Temperature'); ylabel('Rating');


end

% linkfun = @(p, x) pain_predict_function(p, x);
% [p,errval,fit,linkfun,fhan] = nonlin_fit(ratings,temps, 'link', linkfun, 'start', [3 1 2], 'plot', 'verbose');