function y = nlhrf(params,x,varargin)
% function y = nlhrf(params,x,[opt] stick function to convolve)
%
% params are, right now, [height delay dispersion]
% x should be index values of time bin at each data point.
% y should be a column vector of observed values.
%
% Basic use:
% beta = lsqcurvefit(@nlhrf,[1 6 1],x,y')
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

mylen = length(x);

height = params(1);
delay = params(2);
udelay = 16;
disp = params(3);
udisp = 1;
rtou = 6;
onset = 0;
klength = 32;

%whos x

normhrf = spm_hrf(.5,[delay udelay disp udisp rtou onset klength]) ./ ...
    max(spm_hrf(.5,[delay udelay disp udisp rtou onset klength]));

%whos normhrf             

funct = height * normhrf;

if nargin > 2  % convolve
	stick = varargin{1};
    funct = conv(stick,funct);
    if size(funct,2) > size(funct,1), funct = funct';,end
end

y = funct(x);

%whos x

return      
