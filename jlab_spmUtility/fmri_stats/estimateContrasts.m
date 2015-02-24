function estimateContrasts(c,names,type,SelectedContrast,varargin)
% function estimateContrasts(c,names,type,SelectedContrast,myContrastName)
%
% c is cell array of contrasts
% SelectedContrast is a flag for estimating only a selected contrast (vs all)
% names is a cell array of contrast names
% type is a cell array of contrast types - 'T' or 'F'
%
% Tor Wager, adapted from KUL scripts

if length(varargin)>0, myContrastName = varargin{1};,end

% *******************************************************************
% 'ATTENTION : Make sure that you are in the correct results directory!!!'
% *******************************************************************

% ===================================================================
% function 'f99_GetCont_first'
% ===================================================================
% * is a function that helps you to define easily all your contrasts
%   in one individual ('fixed effects analysis')
% * the main contrasts for a 'second level analysis' (on population
%   level) can also be determined this way (i.e. one con???.img per
%   condition and per subject (not session))
% PARAMETERS
% - matrix containing the string names of the contrasts
% - matrix with the values for each defined contrast(one session)
% - the number of sessions (runs) for this subject
% ===================================================================
%f99_GetCont_first(...
%    ['sid-fix';'sid-dim';'lin_pos'],...
%    [ 1  1 1 1 1  0 -5;...
%      1  1 1 1 1 -5  0;...
%     -2 -1 0 1 2  0  0],...
%     8)

% ===================================================================
% function 'f99_GetCont_second'
% ===================================================================
% is a function that helps you to define the same contrasts in all
% subjects for further use in a 'RFX-analysis', if you want one
% con???.img per condition and per session
% PARAMETERS
% - matrix containing the string names of the contrasts(one session)
% - matrix with the values for each defined contrast
% - the number of sessions (runs) for this subject
% ===================================================================

disp(['	...Making xCon.mat'])
   disp('	-----------------------------------------------')
   % eval(['makecon_model' model])
   
	load SPM.mat
	try
        % spm 99 version
        sX = xX.xKXs;
    catch
        % SPM2 analysis
        N = fieldnames(SPM);
        for i = 1:length(N), eval([N{i} ' = SPM.' N{i} ';']);,end
        try
            sX = xX.xKXs;
        catch
            disp('No smoothed/filtered model.  Is this a basic model? Assuming it is.');
            sX = xX;
        end
        
    end

	% load xCon.mat
	% =================================
if exist(fullfile('.','xCon.mat'),'file'), 
	load('xCon.mat'), 
	lxCon = length(xCon);
else, 
   disp('	...Starting new xCon matrix.')
   xCon = spm_FcUtil('FconFields');
	lxCon = 0;
end

% estimate contrasts

for CONS = 1:size(c,2)
   disp(['	...making contrast ' num2str(CONS)])
   % add zeros to end
 	c{CONS} = [c{CONS} zeros(size(c{CONS},1),size(sX.X,2)-size(c{CONS},2))];
     
    % SPM 99 version
	%contrast = spm_FcUtil('Set',names{CONS},type{CONS},'c',c{CONS}',sX);

    % SPM2
    contrast = spm_FcUtil('Set',names{CONS},type{CONS},'c',c{CONS}',sX.X);
    
	%iFc2 = spm_FcUtil('In', contrast, sX, xCon);
    iFc2 = spm_FcUtil('In', contrast, sX.X, xCon);
    
   if ~iFc2, 
      %xCon(length(xCon)+1) = contrast;
      xCon(lxCon+1) = contrast;
   else 
      %- 
      fprintf('\ncontrast %s (type %s) already in xCon', names{CONS}, type{CONS});
   end
   lxCon = length(xCon);
end

disp(['	...saving xCon.mat in ' pwd])
disp('	-----------------------------------------------')

try
	save('xCon.mat','xCon')
catch
	str = ['Can''t write xCon.mat to the results directory: ' swd];
	warning(str);
%	status.str = str;
%	status.err = 3;
end

disp(['	...Estimating contrasts'])
disp('	...-----------------------------------------------')   

% ===================================================================
% function 'f99_DoCont'
% ===================================================================
% * this function computes contrast and statistical images
% * there are 3 options for argunment 1:
%     - 'all'      : all contrasts in xCon.mat will be computed
%     - 'selected' : those contrasts defined in arg.2 will be computed
%                    the exact names of the contrasts should be used!!
% PARAMETERS
% - option ('all' or 'selected')
% - no 2nd argument in case of 'all'
% - contrast names (as defined in f99_GetCont_...)
% ===================================================================
if ~(exist('SelectedContrast') == 1), SelectedContrast = 0;, end

if SelectedContrast
	f99_DoCont(...
	'selected',...
	[myContrastName] ...
	)
else
	f99_DoCont(...
   	'all' ...
   	)

end

	%f99_DoCont(...
	%'selected',...
	%[names{1};names{2};names{3};names{4};names{5};names{6};names{7}] ...
	%)
 
 %['s1 Session 1';...
  %   's2 Session 3';...
   %  's3 Session 6';...
    % 'sid-fix     ';...
     %'sid-dim     '])
     