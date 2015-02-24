
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

disp('===========================================================')
for C = 1:length(ss)
   ResDir = TOR.allresdirs{C};
   
   disp('Changing to subject directory')
   disp('-----------------------------------------------')
   str = (['cd ' ResDir])
   eval(str)
   pwd
   
   if removeoldstat
   		disp('Removing old analysis files')
   		disp('-----------------------------------------------')
   		!rm *T*; rm*F*;rm con*;rm ess*;rm ncon*;rm sncon*;rm xCon.mat;rm *.ps
   end
      
   disp(['Making xCon.mat: Model' model])
   disp('-----------------------------------------------')
   % eval(['makecon_model' model])
   
	load SPM.mat
	sX = xX.xKXs;


	% load xCon.mat
	% =================================
if exist(fullfile('.','xCon.mat'),'file'), 
	load('xCon.mat'), 
	lxCon = length(xCon);
else, 
   disp('		Starting new xCon matrix.')
   xCon = spm_FcUtil('FconFields');
	lxCon = 0;
end

% estimate contrasts

for CONS = 1:size(c,2)
   disp(['Subject ' num2str(ss(C)) ' model ' model ': making contrast ' num2str(CONS)])
   % add zeros to end
 	c{CONS} = [c{CONS} zeros(size(c{CONS},1),size(sX.X,2)-size(c{CONS},2))];
     
	contrast = spm_FcUtil('Set',names{CONS},type{CONS},'c',c{CONS}',sX);

	iFc2 = spm_FcUtil('In', contrast, sX, xCon);
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
disp('-----------------------------------------------')

try
	save('xCon.mat','xCon')
catch
	str = ['Can''t write xCon.mat to the results directory: ' swd];
	warning(str);
%	status.str = str;
%	status.err = 3;
end

disp(['Estimating contrasts for subject ' num2str(ss(C)) ', model ' model])
disp('-----------------------------------------------')   

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
f99_DoCont(...
   'all' ...
   )

%'selected',...
%   [names{1};names{2};names{3};names{4};names{5};names{6};names{7}] ...
%    )
 
 %['s1 Session 1';...
  %   's2 Session 3';...
   %  's3 Session 6';...
    % 'sid-fix     ';...
     %'sid-dim     '])
     
     
end % end loop thru subjects
