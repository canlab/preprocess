function ss=rnfield(s,oldname,newname)
%RNFIELD Rename Structure Fields.
% RNFIELD(S,OldName,NewName) returns the structure S with the field name
% denoted by the string OldName changed to NewName. Oldname must exist in S.
%
% If OldName and NewName are cell arrays of strings of equal length, each
% field name in OldName is changed to the corresponding element in NewName.
%
% Examples:
%
% RNFIELD(S,fieldnames(S),lower(fieldnames(S))) makes all field names lower
% case as long as they are unique when lowercase
%
% If S.a=pi; S.b=inf; then 
% RNFIELD(S,{'a' 'c'},{'c' 'd'}) produces S.d=pi; S.b=inf;
%
% See also ORDERFIELDS, RMFIELD, ISFIELD, FIELDNAMES.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-02-08, 2006-02-09

if nargin~=3
   error('rnfield:IncorrectNumberofInputArguments',...
         'Three Input Arguments Required.')
end
if ~isstruct(s)
   error('rnfield:IncorrectInputArgument',...
         'First Argument Must be a Structure.')
end
if ~(ischar(oldname)||iscellstr(oldname)) &&...
      ~(ischar(newname)||iscellstr(newname))
   error('rnfield:IncorrectInputArgument',...
         'Last Two Arguments Must be Strings or Cell Strings.')
end
if ischar(oldname)   % convert to cell
   oldname=cellstr(oldname);
end
if ischar(newname)   % convert to cell
   newname=cellstr(newname);
end
nold=length(oldname);
if nold~=length(newname)
   error('rnfield:IncorrectInputArgument',...
         'OldName and NewName Must Have the Same Number of Elements.')
end

fnames=fieldnames(s); % get field names of input structure

for k=1:nold          % do the work and perform error checking
   if ~isequal(oldname{k},newname{k}) % no work if oldname==newname
      if ~isvarname(newname{k})
         error('rnfield:NewFieldNameNotValid',...
               'Invalid New Field Variable Name: %s',newname{k})
      end
      old=strcmp(fnames,oldname{k});
      if ~any(old)
         error('rnfield:OldFieldNameDoesNotExist',...
               'Structure Does Not Contain the Field: %s',oldname{k})
      end
      if any(strcmp(fnames,newname{k}))
         error('rnfield:NewFieldNameAlreadyExists',...
               'Structure Already Contains the Field: %s',newname{k})
      end
      fnames(old)=newname(k); % replace old field name with new one
   end
end
ssize=size(s);
c=struct2cell(s);                          % convert structure to cell
ss=reshape(cell2struct(c,fnames,1),ssize); % rebuild with revised fields