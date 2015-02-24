function [P,nimgs,imnumbers] = getfiles(input,varargin)
% [P,nimgs,imnumbers] = getfiles(input,nimgs [opt])
% 
% Gets list of image names
%
% Inputs can be:
%		name with wildcard, e.g. 'sra*img'
%		basename, w/o numbers, e.g. 'sravol_e3303_100_'
%		single number, to list all images up to that number
%		row vector of image numbers
%		combinations of the above: e.g., basename and vector of images to select.
%
% Outputs:
%		P = list of images in cell array
%		nimgs = how many images in the list
%		imnumbers = the numbers attached to the images selected.
%
% Search path:
%		The function first tries to list each file in the current directory
%		Then tries to list imgs in subdirectories called 'scan#' (# is a number)
%		Then prompts you for the subdirectory name of the images.
%
% Also try using a directory name as input instead of a basename.
%
% by Tor Wager.  last edit, 5/15/01

% figure out what type of input
% ========================================================================================= 
if isstr(input), 
   for i = 1:size(input,2)
   		if input(i) == '*', Wild = 1;,end
   end
   if ~exist('Wild'),input = [input '*img'];,Wild = 0;,end
        
   if nargin > 1,
      nimgs = varargin{1};
      switch size(varargin{1},2)
      case 1
         Selection = 'singlenumber';
      otherwise Selection = 'numberlist';
      end
      if isempty(nimgs),Selection = 'no selection';,end
   else Selection = 'none';
   end

else 
   nimgs = input;
   if isempty(input), Selection = 'no selection';
   elseif size(input,2) == 1, Selection = 'singlenumber';
   else Selection = 'listofnumbers';
   end
   input = ['*.img'];		% so that no wildcard applies in list
end

% =========================================================================================  
	% separate the directory from the filename wildcard.
   lastchar = 0;
   for i = 1:size(input,2)
   		if input(i) == '/' | input(i) == '\', lastchar = i;,end
   end
   if lastchar, directory = input(1:lastchar);,input = input(lastchar+1:end);
   else directory = [];
   end
   

	% list all files in directory with wildcard, or '*.img' if no wildcard
   disp(['!ls ' directory input ' > file.txt'])
   eval(['!ls ' directory input ' > file.txt'])
   !chmod 777 file.txt

   P = textread('file.txt','%s');
   
	% if empty, then list files in scan* directories
   if isempty(P), 
      %disp(['ls ' directory 'scan*/' input ' > file.txt'])
      %eval(['ls ' directory 'scan*/' input ' > file.txt'])
      % modify to work around the argument list too long problem in ls
      %James temporary fix 6/15 for use with 10 scan vnl data. Should either calculate this from existing data or add box for numscans.
      for i = 1:10	% max 20 scans
	disp(['!ls ' directory 'scan' num2str(i) filesep input ' >> file.txt'])
         eval(['!ls ' directory 'scan' num2str(i) filesep input ' >> file.txt'])
      end
      P = textread('file.txt','%s');
   end
   
   % or look in other user-entered directory directories
   if isempty(P), 
      %directory = inputdlg('?>','Can''t find scan directory or wrong wildcard. Enter directory, minus indexing #?')
      %input = inputdlg('?>','Enter filename wildcard?')
      %scandir = input('Enter name of scan subdirectory in single quotes: ');
      %eval(['ls ' directory '/*/' input ' > file.txt'])
      %P = textread('file.txt','%s');
   		return   
   end
   
   % if still empty, then error.
   if isempty(P), error(['Couldn''t find specified files in ' scandir '1/' input]),end
   
     % build imnumbers
     Q = str2mat(P);
     Q = str2num(Q(:,end-7:end-4))';
  
   switch Selection
   case 'singlenumber'
      P = P(1:nimgs,:);
      Q = Q(1:nimgs);
      nimgs = size(P,1);
   case 'numberlist'
      index = 1;
      for i = 1:size(nimgs,2),
         imgexists = sum(Q == nimgs(i));
         if imgexists, Pnew{index,1} = P{Q == nimgs(i)};
            index = index + 1;
         end
      end
      P = Pnew;
      %re-build numbers and nimgs
      nimgs = size(P,1);
      Q = str2mat(P);
     	Q = str2num(Q(:,end-7:end-4))';
   end
   
   nimgs = size(P,1);
   imnumbers = Q;

   
return
   
   
   
   
   
   
   
   
   
% old stuff   
DUMMY = 'noneoftheabove';  
   
switch DUMMY
      
case 'basename_and_number'
% ========================================================================================= 
basename = input;
nimgs = varargin{1};
if nimgs < 9
   for i = 1:nimgs
      str{i} = ([basename '000' num2str(i) '.img']);
	end
else
	for i = 1:9
   	    str{i} = ([basename '000' num2str(i) '.img']);
	end
   if nimgs > 99
		for i = 10:99
   		    str{i} = ([basename '00' num2str(i) '.img']);
		end
   else
      for i = 10:nimgs
   		    str{i} = ([basename '00' num2str(i) '.img']);
		end
   end
	if nimgs > 999
		for i = 100:999
   		    str{i} = ([basename '0' num2str(i) '.img']);
   	    end
   	    for i = 1000:nimgs
   		    str{i} = ([basename num2str(i) '.img']);
		end
	else
   	    for i = 100:nimgs
   		    str{i} = ([basename '0' num2str(i) '.img']);
		end
   end
end
P = str;

      
case 'basename_and_numberlist'
% =========================================================================================  
% for fixed list of image numbers
	basename = input;
	nimgs = varargin{1};
      disp('		Timeseries: fixed list of images')
      for i = 1:size(nimgs,2)
         if nimgs(1,i) < 10
            str{i} = ([basename '000' num2str(nimgs(i)) '.img']);
         elseif nimgs(1,i) < 100, str{i} = ([basename '00' num2str(nimgs(i)) '.img']);
         elseif nimgs(1,i) < 1000, str{i} = ([basename '0' num2str(nimgs(i)) '.img']);
         elseif nimgs(1,i) < 10000, str{i} = ([basename num2str(nimgs(i)) '.img']);
         else error('image numbers from 1 to 10000 only.')
         end
      end
   P = str;
   
case 'singlenumber'
% =========================================================================================  

end

   


   
   

