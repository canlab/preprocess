% read_edat_output
%
% Script that creates variables in the workspace that correspond to columns
% of the edat file output (saved in text tab delimited "excel" format)
%
% Tor Wager
% 
% You must define the filename and the number of columns before running
% fname = 'myfile.txt'; numc = 148; (where 148 is a placeholder for the number of cols)

% ------------------------------------------------------------------------
% read the database initially to get all column names
% number of columns must be predefined
% ------------------------------------------------------------------------

addone = 0;

%numc = 148;
[d] = textread(fname,'%s','headerlines',1,'delimiter','\t');

if strcmp(d{1},'This file contained edited data.'),
    % add an extra input row!
    disp('edited data: extra row.')
    addone = 1;
    [d] = textread(fname,'%s','headerlines',2,'delimiter','\t');
end

names = d(1:numc)'; disp(['Last name column read is ' names{end}])

% ------------------------------------------------------------------------
% replace illegal characters that will cause Matlab problems
% ------------------------------------------------------------------------

for i = 1:length(names), names{i}(findstr('.',names{i})) = '_';, end
for i = 1:length(names), names{i}(findstr('[',names{i})) = '';, end
for i = 1:length(names), names{i}(findstr(']',names{i})) = '';, end
for i = 1:length(names), names{i}(findstr(':',names{i})) = '';, end
for i = 1:length(names), names{i}(findstr('=',names{i})) = '';, end

% ------------------------------------------------------------------------
% create formatting string with number of columns and output var names
% ------------------------------------------------------------------------

fmt = repmat('%s',1,numc);
fmt = [fmt '%*[^\n]'];

outs = ['[' names{1}];
for i = 2:numc, outs = [outs ',' names{i}];, end
outs = [outs ']'];

% works, but requires an extra step
%outs = ['[out{1}'];
%for i = 2:numc, outs = [outs ',out{' num2str(i) '}'];, end
%outs = [outs ']'];

% ------------------------------------------------------------------------
% read the database again with proper formatting and output names
% ------------------------------------------------------------------------

if addone
    str = [outs ' = textread(''' fname ''',fmt,''headerlines'',3,''delimiter'',''\t'');'];
else
    str = [outs ' = textread(''' fname ''',fmt,''headerlines'',2,''delimiter'',''\t'');'];
end

eval(str)

warning off
ww = whos('*RT*'); ww = ww(end).name; eval(['ww = isempty(str2num(' ww '{1}));'])
warning on
if ww,
    str = [outs ' = textread(''' fname ''',fmt,''headerlines'',3,''delimiter'',''\t'');'];
    eval(str)
end

% not necessary if all filenames are OK
%for i = 1:length(outs),
%    eval([names{i} ' = out{i};'])
%end


% ------------------------------------------------------------------------
% convert blanks to NaN's, to avoid losing placeholders, and
% convert columns with numeric information to numeric vectors
% ------------------------------------------------------------------------
for i = 1:numc
    
    % The replacing blanks part
    eval(['wh = strmatch('' '',str2mat(' names{i} '{:}));'])
    eval([names{i} '(wh) = {''NaN''};'])
    
    % the conversion to numbers part
    eval(['a = str2num(str2mat(' names{i} '{:}));'])
    if isempty(a) | sum(isnan(a)) == length(a)
        % leave it alone; it's text
    else
        eval([names{i} ' = a;']), 
    end
end

