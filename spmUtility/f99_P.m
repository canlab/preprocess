function [P] = f99_P(DIR,EXP)
%----------------------------------------------------------------------------
% Makes a matrix P which contains the filenames of the images
%============================================================================
global fmriDIR

%hold the position of the current working directory
cwd_old = pwd;

% Change the working directory to the base fmriDIR
%-------------------------------------------------
cwd = [deblank(fmriDIR)];
str = ['cd ' cwd];
eval(str);
CWD = pwd;

[nd, tmp] = size(DIR);
P = [];
for i = 1:nd
  [files, tmp] = spm_list_files(DIR(i,:),EXP(i,:));
  [nf, tmp] = size(files);
  lfs = '';
  for j = 1:nf
    lf = [DIR(i,:) filesep files(j,:)];
    lfs = str2mat(lfs,lf);
  end
  P = str2mat(P,lfs(2:nf+1,:));
end
[np, tmp] = size(P);

%go back to old current working directory 
str = ['cd ' cwd_old];
eval(str);
CWD = pwd;

P = P(2:np,:);

