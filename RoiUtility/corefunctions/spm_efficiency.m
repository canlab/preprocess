function O = spm_efficiency(varargin)
% function O = spm_efficiency([SPM.mat file name],[vec. of contrasts of interest])
% Tor Wager  9/30/02
%
% Loads an SPM.mat file (prompts if no inputs entered)
% Calculates efficiency for contrast estimation
% And for HRF shape estimation (given stick f's in Sess)
% Prints an output table.
%
% Functions Called:
% spm_get and other spm functions
% iofun toolbox
% efficiency.m
%   calcEfficiency.m
%   designvector2model.m
%       sampleInSeconds.m
%       getPredictors.m
%       resample.m    (toolbox; free version available from Sourceforge.net)
%       modelSaturation.m
% tor_make_deconv_mtx3.m
% fullpath.m
%
% contrasts should be a matrix of row vectors for each contrast


if ~isempty(varargin)
    Pspm = varargin{1};
    Pspm2 = getfullpath(Pspm);
    if isempty(Pspm2), Pspm2 = which(Pspm); end
    if isempty(Pspm2), Pspm2 = which(Pspm); end
    Pspm = Pspm2;
else
    Pspm = spm_get(1,'SPM.mat','Choose SPM.mat file');
end

if isempty(Pspm), error(['File ' Pspm ' not found.']),end

disp(['Loading: ' Pspm])
eval(['load ' Pspm])

% for SPM2/5
if ~exist('xX', 'var'), xX = SPM.xX; end
if ~exist('Sess', 'var'), Sess = SPM.Sess; end

switch spm('Ver')
    case {'SPM2', 'SPM5', 'SPM8'}
        O.X = xX.X;
        
    otherwise
        O.X = xX.xKXs.X;
end
X = O.X;
nrows = size(X, 1);

[d f e] = fileparts(Pspm);
xconname = [d filesep 'xCon.mat'];

% Load contrasts
if exist(xconname, 'file')
    disp(['Loading: ' d filesep 'xCon.mat'])
    eval(['load ' d filesep 'xCon.mat'])

    disp('Contrasts are:')
    str2mat(xCon.name)

    if length(varargin) > 1, myint = varargin{2}; else myint = 1:length(xCon); end
    fprintf(1,'\nContrasts of interest are:\n')
    str2mat(xCon(myint).name)
    mycons = cat(2,xCon(myint).c);  % contrasts are columns in
    for i = 1:length(myint), mynames{i} = xCon(myint(i)).name; end

else
    disp('No contrasts. Using regressors of interest');
    
    n_of_interest = size(SPM.xX.iC, 2);
    mycons = [eye(n_of_interest); zeros(size(X, 2) - n_of_interest, n_of_interest)];
    
    mynames = SPM.xX.name(SPM.xX.iC);
end


% ---------------------------------------------------------
% * power in design stuff (how to choose a HP filter cutoff for an ER design)
% ---------------------------------------------------------
mm = [];

x = abs(fft(X(:,1:max(xX.iC)))); x = x(1:round(length(x) ./ 2),:);  % this one for all predictors
x = cumsum(x);
for i = 1:size(x,2), x(:,i) = x(:,i) ./ max(x(:,i)); end   % don't use repmat - no toolboxes

% nyquist = TR / 2;  
% fft elements are frequencies, last one is 2*nyquist = TR; divide el index by (TR*len) so last el = TR
if isfield(xX, 'RT')
    % old SPM
    TR = xX.RT;
else
    % new spm: spm5, 2?
    TR = SPM.xY.RT;
end

f = 1:nrows; f = f * 1/(TR * nrows);
create_figure('spm_efficiency', 3, 2);
plot(f(1:length(x)),x); set(gcf,'Color','w'); title('Cumulative power of all predictors','FontSize',14);
for i = 1:size(x,2); m(i) = sum(x(:,i) < .2); end
axis([0 f(max(m)) 0 1]); xlabel('Hz','FontSize',14)

%for i = 1:size(x,2); m2(i) = sum(x(:,i) < .1);, end
%for i = 1:size(x,2); m3(i) = sum(x(:,i) < .05);, end
%for i = 1:size(x,2); m4(i) = sum(x(:,i) < .01);, end
clear m

% now get contrasts of interest, which is critical for choosing cutoff
% -----------------------------------------------------------------------
conX = X * mycons;

x = abs(fft(conX)); x = x(1:round(length(x) ./ 2),:);
x = cumsum(x);
for i = 1:size(x,2), x(:,i) = x(:,i) ./ max(x(:,i)); end   % could use repmat instead

f = 1:nrows; f = f * 1/(TR * nrows);
create_figure('spm_efficiency', 3, 2, 1);
subplot(3, 2, 2);
plot(f(1:length(x)),x); 
title('Cumulative power of contrasts of interest','FontSize',14);
for i = 1:size(x,2); m(i) = sum(x(:,i) < .2); end
axis([0 f(max(m)) 0 1]); xlabel('Hz','FontSize',14)
if length(mynames) < 10, legend(mynames); end

for i = 1:size(x,2); m2(i) = sum(x(:,i) < .1); end
for i = 1:size(x,2); m3(i) = sum(x(:,i) < .05); end
for i = 1:size(x,2); m4(i) = sum(x(:,i) < .01); end

% m is how many elements fall below the power cutoff level 
% use it to find frequency value, since this is a cumulative distribution function 
% m = 0 if all power, or enough power, is at the lowest frequency!

m(m == 0) = 1; m2(m2 == 0) = 1; m3(m3 == 0) = 1; m4(m4== 0) = 1;    % avoid looking for 0th element.
mm = [f(min(m)) f(min(m2)) f(min(m3)) f(min(m4))];
mm = mm'; mm(:,2) = mm;  mm(:,1) = [20 10 5 1]'; mm(:,3) = 1 ./ mm(:,2); mm = mm';

fprintf(1,'\nEfficiency for design stored in %s\n',Pspm); 
fprintf(1,'___________________________________________________\n')
fprintf(1,'Power in regressors\n')
fprintf(1,'Less than %3.0f%% power in all predictors at %3.3f Hz, %3.1f s\n',mm)
fprintf(1,'\n')

%fprintf(1,'\n5%% and 1%% Power-loss Levels for the Lowest-Power Predictor by Session (in s)\n')
%fprintf(1,'(suggested cutoff values)\n')

%clear mm; cumnr = 1;mystr = [];

%for i = 1:length(Sess)
%    nr = (Sess{i}.col); 
%    mall = [f(min(m(nr))) f(min(m2(nr))) f(min(m3(nr))) f(min(m4(nr)))];
%    mm(i,:) = mall(3:4);    % take the 5% value and the 1% value
%    mystr = [mystr '%3.1f\t'];
%end
%mm = mm'; 
%mm = 1 ./ mm;   % convert to s
%mm = [[5;1] mm];
%mm = mm';

%eval(['fprintf(1,''%3.0f%% max power loss: ' mystr '\n'',mm)'])
%fprintf(1,'\n')  


% ---------------------------------------------------------
% * recommended cutoff values with a-optimality criterion
% the idea is to find the optimal tradeoff between signal
% and noise eliminated from the design by the filter
% ---------------------------------------------------------
% this way of doing autocorr does not take TR into account
%load hiautocorr
%myscannerxc = myscannerxc - myscannerxc(end);

% bottom line: how efficient it is to filter really depends on the 
% amount of noise coloration!

f = inline('A * exp(-a * x)','A','a','x');  % exponential function
xa = 0:TR:30;
myscannerxc = f(.9475,.5323,xa);    % params at 3T from vnl experiment, n = 10

create_figure('spm_efficiency', 3, 2, 1);
subplot(3, 2, 3);
plot(xa,myscannerxc), title('Assumed autocorrelation'),xlabel('Time (s)'),drawnow
disp('Finding recommended HP cutoff values with A-optimality criterion')
disp('Getting autocorrelation matrix')
V = getv('make',myscannerxc,size(xX.X,1));

it = round(10):10:round(mm(3,3));
it = round(10:10:90);
fprintf(1,'Iterating over: %s ', num2str(it))
fprintf(1,'\n')
for i = 1:length(it)
    fprintf(1,'%3.0f.',it(i))
    [S] = use_spm_filter(TR,nrows,'none','specify',it(i));
    xtxitx = pinv(S * X); svi = S * V;
    [e,ev] = calcEfficiency(ones(1, size(mycons,2)),mycons',xtxitx,svi);
    eff(i) = e(1); effv(:,i) = ev;
    
    svi = S * eye(nrows);
    [e2,ev2] = calcEfficiency(ones(1, size(mycons,2)),mycons',xtxitx,svi);
    eff2(i) = e2(1); effv2(:,i) = ev2;
end
fprintf(1,'\n')
create_figure('spm_efficiency', 3, 2, 1)
subplot(3, 2, 4);
plot(it,eff,'b','LineWidth',3), hold on; plot(it,eff2,'r','LineWidth',3),
legend({'Efficiency - Colored noise' 'Efficiency - white noise'})
xlabel('HP filter cutoff (Sec)');
ylabel('Efficiency');

a = it(eff==max(eff)); a = a(1); b = it(eff2==max(eff2)); b = b(1);
fprintf(1,'\nRecommended cutoff is %3.0f s (colored) and %3.0f (white scanner noise)\n',a,b);

% ---------------------------------------------------------
% * Efficiency of params and intercept
% ---------------------------------------------------------

[O.peff]= efficiency(X,O);  % before defining contrasts

if iscell(Sess)
    % Old SPM
    O.pname = Sess{1}.name;
else
    % New SPM (5, 2?)
    nregs = length(Sess(1).Fc);
    O.pname = cell(1, nregs);
    for ii = 1:nregs, O.pname{ii} = Sess(1).Fc(ii).name; end
end

O.pint = O.peff(end-length(Sess)+1:end);    % intercepts

%SS edit 12/13/2010
for n = 1:numel(Sess)
    ndesign(n) = numel(Sess(n).Fc);
end
for n = 1:numel(Sess)-1
    nequal(n) = ~isequal(ndesign(1),ndesign(n+1));
end
nequal = sum(nequal);
if nequal > 0
    idx = 0;
    O.pname = cell(1,numel(Sess));
    for n = 1:numel(Sess)
        for ii = 1:ndesign(n), O.pname{n}{ii} = Sess(n).Fc(ii).name; end
        tempeff{n} = O.peff(idx+1:idx+numel(Sess(n).Fc));
        idx = idx + numel(Sess(n).Fc);
    end
    O.peff = tempeff;
else
    O.peff = reshape(O.peff(1:length(Sess)*length(O.pname)),length(O.pname),length(Sess));  %this breaks when there are unequal component numbers across sessions in the design
end
%end edit




% ---------------------------------------------------------
% * contrast and HRF setup stuff
% ---------------------------------------------------------

O.contrasts = mycons;
O.names = char(mynames{:});  %str2mat(xCon(2:end).name);
O.HRFtime = 30;     % time in seconds to estimate FIR model (deconvolution)
O.Vi = [];
O.ISI = TR;
O.TR = TR;


% ---------------------------------------------------------
% * HRF and contrast calculations
% for each session
% ---------------------------------------------------------

if iscell(Sess)
    % Old SPM
    allsf = Sess{i}.sf;
else
    % New SPM (5, 2?)
    allsf = full(cat(2, Sess(1).U(:).u));
    allsf = mat2cell(allsf, size(allsf, 1), ones(1, length(Sess(1).U)));
end

for i = 1:length(Sess)
    [DX,sf] = tor_make_deconv_mtx3(allsf,round(O.HRFtime ./ TR),16); % sf is resampled at TR
    O.delta = cat(2,sf{:});
    sf2{i} = O.delta;
    [ceff,hrfeff{i}]= efficiency(conX,O);
    hrfintercept(i) = hrfeff{i}(end);
    hrfeff{i} = mean(reshape(hrfeff{i}(1:end-1),round(O.HRFtime ./ TR), nregs));
    
    % cols are conditions, rows are time points
end

O.delta = cat(1,sf2{:});            % overall stick function
[ceff,hrfeff{end+1}]= efficiency(conX,O);
hrfintercept(end+1) = hrfeff{end}(end);
hrfeff{end} = mean(reshape(hrfeff{end}(1:end-1),round(O.HRFtime ./ TR), nregs));

O.hrfeff = hrfeff;
O.ceff = ceff;

% ---------------------------------------------------------
% * table
% ---------------------------------------------------------


% HRF
% ---------------------------------------------------------

fprintf(1,'Regressor Beta Parameter and HRF Efficiency Estimates\n')
fprintf(1,'\tCondition\t');
for j = 1:length(Sess), fprintf(1,'Sess %3.0f\t\t',j); end
fprintf(1,'Overall\t')
fprintf(1,'\n')
fprintf(1,'\t\t');
for j = 1:length(Sess), fprintf(1,'Param\tHRF\t');,end
fprintf(1,'\n')

%SS edit 12/13/2010
if nequal > 0
    for i = 1:length(O.pname)
%         fprintf(1,'\t%s\t',O.pname{i});
        for j = 1:numel(O.pname{i})
            fprintf(1,'\t%s\t',O.pname{i}{j})
            fprintf(1,'%3.3f\t',O.peff{i}(j))
        end
        for j = 1:length(Sess)
            fprintf(1,'%3.3f\t',hrfeff{j}(i))
        end
        
        fprintf(1,'%3.3f\t',hrfeff{end}(i))
        
        fprintf(1,'\n')
    end
else
    for i = 1:length(O.pname)
        fprintf(1,'\t%s\t',O.pname{i});
        
        for j = 1:length(Sess)
            fprintf(1,'%3.3f\t',O.peff(i,j))
            fprintf(1,'%3.3f\t',hrfeff{j}(i))
        end
        
        fprintf(1,'%3.3f\t',hrfeff{end}(i))
        
        fprintf(1,'\n')
    end
end
%end edit

fprintf(1,'\tIntercept\t');
    for j = 1:length(Sess)
        fprintf(1,'%3.3f\t%3.3f\t',O.pint(j),hrfeff{j}(end))
    end
fprintf(1,'%3.3f\t',hrfeff{end}(end))
  
    
    
% Contrasts
% ---------------------------------------------------------

fprintf(1,'\n')
fprintf(1,'Contrast Efficiency Estimates\n')
fprintf(1,'\tContrast\tEfficiency\tStd. Err (Design component)\tEff / Intercept \n');
for i = 1:size(O.names,1)
    fprintf(1,'\t%s\t%3.3f\t%3.3f\t%3.3f\n',O.names(i,:),ceff(i),(1./ceff(i)).^.5,ceff(i)./sum(O.pint))
end


return




