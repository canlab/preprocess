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
f99_GetCont_first(...
    ['sid-fix';'sid-dim';'lin_pos'],...
    [ 1  1 1 1 1  0 -5;...
      1  1 1 1 1 -5  0;...
     -2 -1 0 1 2  0  0],...
     8)

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
f99_GetCont_second(...
    ['s1';'s2';'s3';'s4';'s5';'di'],...
    [ 1 0 0 0 0 0 -1;...
      0 1 0 0 0 0 -1;...
      0 0 1 0 0 0 -1;...
      0 0 0 1 0 0 -1;...
      0 0 0 0 1 0 -1;...
      0 0 0 0 0 1 -1],...
     8)


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
    'selected',...
    ['s1 Session 1';...
     's2 Session 3';...
     's3 Session 6';...
     'sid-fix     ';...
     'sid-dim     '])

