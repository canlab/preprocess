function myFiles = tor_list_files(dirs,imgwildcard)
% function myFiles = tor_list_files(dirs,imgwildcard)
%
% input:
%	dirs: cell array of directory names
%		one directory = one session per cell
%
%	imgwildcard: wildcard of image files to get.
%		example - 'sravol*img'
%
% output: matrix of all files, padded with zeros
% 	  separate cells of output are separate sessions
%
% uses spm_list_files.m from SPM99 distribution
% returns full path names
% 
% Tor Wager, 10/17/01


% get the list of actual file names (P in SPM)
% ====================================================================== 
if ~iscell(dirs), d{1} = dirs;, dirs = d;,end
startdir = pwd;

myFiles = [];
for j = 1:length(dirs)

	[FNames,dummy] = spm_list_files(dirs{j},imgwildcard);

	% ...and add the directory name to get the full path
 	a = repmat([startdir filesep dirs{j} filesep],size(FNames,1),1);
	FNames = [a FNames];

    FNames = str2mat(FNames);
    FNames(FNames == '/') = filesep;
    FNames(FNames == '\') = filesep;
    
	myFiles{j} = FNames;
    
    
    
end

return