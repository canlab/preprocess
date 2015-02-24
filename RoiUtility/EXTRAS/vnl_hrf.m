function hrf = vnl_hrf(TR,x)
% hrf = vnl_hrf(TR,x)
% 
% HRF modified for nonlinearity, depending on the stimulus position (x) in
% a series.  Values derived for series of stimuli 1 s apart.
%
% Uses spm_hrf.m from SPM's Wellcome group at the FIL / UCL
% x of 1 is the first occurrence of a stimulus, followed by x = 2 is the
% second in a series, etc.
%
% tor wager, 10/04
%
% example:
% figure('Color','w'); set(gca,'FontSize',18)
% col = {'r' 'r--' 'g' 'g--' 'b' 'b--' 'k'};
%xx = [1 2 5 6 10 11 30];  % positions
%for x = 1:7,
%   hrf = vnl_hrf(.1,xx(x));hold on; plot(hrf,col{x},'LineWidth',2);
%end
% legend({'1st' '2nd' '5th' '6th' '10th' '11th' '30th'})
%


% equations for height, onset delay, time to peak as a function of stimulus
% position x

heighteq = inline('1.7141.*(exp(-2.1038.*x)) + 0.4932.*(exp(-0.0770.*x))');
delayeq = inline('-13.4097.*(exp(-1.0746.*x)) + 4.8733.*(exp(-0.1979.*x))');
peakeq = inline('37.5445.*(exp(-2.6760.*x)) + -3.2046.*(exp(-0.2120.*x)) + 5.6344');

m = heighteq(x);    % height
d = delayeq(x);     % delay
p = peakeq(x);      % time to peak

% SPM gamma function, normalized to area = 1
%gamm = inline('el.^h .* x.^(h-1) .* exp(-el .* x) ./ gamma(h)','x','h','el');
%res = TR./16;                       % resolution, bins per sec
%x     = [0:(32/res)] - d/res;    % x - d
%hrf = gamm(x,p,res);         % gamma function

% Use SPM HRF instead because it avoids roundoff error and edge effects
% where the function is invalid at different delays 
%delay = d; 
%peak = p; 

uonset = 16;
disp = 1; 
udisp = 1; 
rtou = 6; 
klength = 32;

p = [p uonset disp udisp rtou d klength];

hrf = spm_hrf(TR);


hrf = m .* hrf ./ max(hrf);     % normalize height to % of max height


return

