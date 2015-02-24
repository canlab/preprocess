function DATA=scaledata(DATA)
% function DATA=scaledata(DATA)
%
% Purpose: to mean-center and windsorize data within conditions
%          recursive trimming, 3 iterations
%          And to perform high-bass filtering if necessary
%
% Inputs:
% DATA, data structure from the mvroi tool
%       requires: DATA.DATA.dat field, (conditions x time) within 
%       {subjects x regions} cell array
%
% DATA.SPEC.trim, number of standard deviations to trim to
%       zero or empty for no trimming
%
% filter options: define in DATA.SPEC
% see filterAdjust for filter options and structure definition
%
% experiment-critical fields for filtering, with example values:
% DATA.SPEC.HP = 120;
% DATA.SPEC.TR = 1.5;
%
% for no filtering, leave out HP field.
%
% Outputs:
% DATA.DATA.filtered_dat
%
% calls: scale (row normalization)
%        filterAdjust (high pass filtering & trimming)


% newer vector of cells for each subject format
numsub=length(DATA.DATA.dat);    %number of subjects
numreg=size(DATA.DATA.dat{1},2);  %number of ROIs

if ~isfield(DATA.SPEC,'trim'),
    trim = [];
else
    trim = DATA.SPEC.trim;
end

% set trimming default to 3sd
if trim > 0 & ~isempty(trim)            
	fprintf(1,'Removing session means . Windsorizing, %3.0f std . ',trim);   
end

% here we column-center data, making sure the timeseries from each region has
% the same mean.  Unnecessary, but defines the output structure.

for s=1:numsub
    for r=1:numreg
         dat = DATA.DATA.dat{s}(:,r);
         dat = scale(dat,1);
         DATA.DATA.filtered_dat{s}(:,r) = dat;
     end    
end

% Add high-pass filtering and trimming

% defaults; also removes session means

y = ones(size(DATA.DATA.filtered_dat{s}(:,r)));  % fake data to set up filter

[y,I,S] = hpfilter(y,DATA.SPEC.TR,DATA.SPEC.HP,DATA.SPEC.spersess);

% old way, doesn't do sessions of different lengths
%O = struct('filtertype','spm','doHP',1,'HP',DATA.SPEC.HP,'nruns',DATA.SPEC.nruns,'scanadjust',1, ...
%'percent',0,'plot',0,'trimts',trim,'TR',DATA.SPEC.TR,'verbose',0,'cyclecorrection2',0);

% add options from DATA.SPEC, replacing defaults
N = fieldnames(DATA.SPEC);
for i = 1:length(N), eval(['O.' N{i} ' = DATA.SPEC.' N{i} ';']), end

if ~isfield(O,'HP'), 
    % no filtering specified in DATA.SPEC, return
    return
end

fprintf(1,'High-pass filtering, %3.0f s . \n',O.HP);


% here we do the actual filtering.  All the above was just set-up and
% extras.

for s=1:numsub
    if length(DATA.SPEC.nruns)>1;   %if different numbers of runs for each subject
        nruns=DATA.SPEC.nruns(s);
        %O.nruns=nruns;
    end
    for r=1:numreg       
        %O.y = DATA.DATA.filtered_dat{s}(:,r);   
        %dat = filterAdjust(O)';
        dat = hpfilter(DATA.DATA.filtered_dat{s}(:,r),[],S,DATA.SPEC.spersess,I);
        [dat,DATA.DATA.ntrimmed(s,r)] = trimts(dat,3,[],1,3);      % 3 std, 3 times, with spike correction
        
        DATA.DATA.filtered_dat{s}(:,r) = dat;              
    end
end


return








return

