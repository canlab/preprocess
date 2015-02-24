function x = nlhrf2_TR1(params,xdata)
% function x = nlhrf2_TR1(params,xdata)
%
% params are hrf parameters [height delay intercept]
% xdata     is the delta function to convolve with the HRF
%
% fixed for TR of 1 s!
%
% Basic use:
% beta = lsqcurvefit(@nlhrf,[1 6 0],x,y')
%
%
% To generate and test with sample data:
% ---------------------------
% x1 = .1:.2:10*6.28;y = sin(x) + rand(1,length(x));x = x(1:30);y = y(1:30);
% x = 1:30; % x should be integer values to index time since trial onset (TRs).
% x = [x x x x x x x x x x];
% y = sin(x1(1:length(x))) + rand(1,length(x));  
% beta = lsqcurvefit(@nlhrf,[1 6 1],x,y');fit = nlhrf(beta,x);figure;plot(y);hold on;plot(fit,'r')
%
% To get confidence intervals:
% ---------------------------
% [beta,resnorm,residual,exitflag,output,lambda,J]= lsqcurvefit(@nlhrf,[1 6 1],x,y');
% ci = nlparci(beta,residual,J);
%
% To get confidence intervals AND use Levenberg-Marquardt method:
% ---------------------------
% myoptions = optimset('LevenbergMarquardt','on');
% [beta,resnorm,residual,exitflag,output,lambda,J]= lsqcurvefit(@nlhrf,[1 6 1],x,y',-Inf,Inf,myoptions);
% ci = nlparci(beta,residual,J);    
%
% to have the function consist of MORE than one impulse response
% and convolve the IR function with a delta function before returning
% enter IR stick function as third argument.
%
%

global TR

len = length(xdata);
%TR = xdata(end);
xdata = xdata(1:end-1);

height = params(1);
delay = 6;
%delay = params(2);
udelay = 16;
disp = 1;
udisp = 1;
rtou = 6;
%onset = 0;
onset = params(2);
klength = 32;

normhrf = spm_hrf(TR,[delay udelay disp udisp rtou onset klength]);
normhrf = normhrf ./ max(normhrf);

%whos normhrf             

funct = height * normhrf;
funct = conv(xdata,funct);
if size(funct,2) > size(funct,1), funct = funct';,end

x = funct(1:len);

x = x + params(3);

return      
