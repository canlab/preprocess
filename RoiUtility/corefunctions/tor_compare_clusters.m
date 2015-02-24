function [clusters1,clusters2,test,SPM1,SPM2,VOL1,VOL2] = tor_compare_clusters(varargin)
% function [clusters1,clusters2,test,SPM1,SPM2,VOL1,VOL2] = tor_compare_clusters(SPM1,VOL1,SPM2,VOL2)
% arguments optional.
% This function does a t-test at each value for a set of clusters structures.
% must use either 0, 2 or 4 input arguments
% inputs must be structures with XYZ and XYZmm pointlists

try
	% everything

if nargin > 0, SPM1 = varargin{1};, VOL1 = varargin{2};,end
if nargin > 2, SPM2 = varargin{3};, VOL2 = varargin{4};,end

test.data{1} = [];
test.data{2} = [];

% ------------------------------------------------------------------------------------
% * get SPM structure, image names, and point list for all voxels
% ------------------------------------------------------------------------------------

% > first contrast
% --------------------------------------------------------

if ~(exist('SPM1') == 1)
    disp('Select SPM.mat file for contrast 1')
    [SPM1,VOL1] = spm_getSPM;
end

if isfield(SPM1,'imnames')
    test.imnames{1} = SPM1.imnames;
else
    test.imnames{1} = spm_get([1 1000],'*.img','Select images for 1st contrast',pwd,0);
end

SPM1.imnames = test.imnames{1};

% > second contrast
% --------------------------------------------------------

if ~(exist('SPM2') == 1)
    disp('Select SPM.mat file for contrast 2')
    [SPM2,VOL2] = spm_getSPM;
end

if isfield(SPM2,'imnames')
    test.imnames{2} = SPM2.imnames;
else
    test.imnames{2} = spm_get([1 1000],'*.img','Select images for 2nd contrast',pwd,0);
end

SPM2.imnames = test.imnames{2};


% > point list
% --------------------------------------------------------
test.XYZ = [SPM1.XYZ SPM2.XYZ];




% ----------------------------------------------------------------------------------
% get the timeseries for all voxels in combined map
% ----------------------------------------------------------------------------------
fprintf(1,'Extracting timeseries values for contrast 1.')
O.coords = test.XYZ';
ts = timeseries2('multi',test.imnames{1},O);
test.ts{1} = ts.avg;
test.data{1} = ts.indiv;

fprintf(1,'...contrast 2\n')
ts = timeseries2('multi',test.imnames{2},O);
test.ts{2} = ts.avg;
test.data{2} = ts.indiv;

% ----------------------------------------------------------------------------------
% do t-test
% ----------------------------------------------------------------------------------
if size(test.data{1},2) ~= size(test.data{2},2), error('Cl timeseries not same length.'),end

for i = 1:size(test.data{1},2)
	est = test.data{1}(:,i) - test.data{2}(:,i);
	[h,test.p(i),ci,test.t(i),ser] = t_test2(est);
   	%
	% THIS IS WRONG.
	% [test.t(i), test.p(i)] = t_test(test.data{1}(:,i),test.data{2}(:,i));
end


SPM1.t = test.t(1:size(SPM1.XYZ,2));
SPM1.p = test.p(1:size(SPM1.XYZ,2));
SPM2.t = test.t(size(SPM1.XYZ,2)+1:end);
SPM2.p = test.p(size(SPM1.XYZ,2)+1:end);


% ----------------------------------------------------------------------------------
% define cluster structure for each cluster
% ----------------------------------------------------------------------------------

clusters1 = define_cluster(SPM1,VOL1);
clusters2 = define_cluster(SPM2,VOL2);




% ----------------------------------------------------------------------------------
% print table
% ----------------------------------------------------------------------------------
myTit = [SPM1.title '_' SPM2.title];
myTit(myTit == ' ' | myTit == ':') = '_';

eval(['diary ' myTit '.txt'])

disp(['Contrast 1: ' SPM1.title '	Contrast 2:' SPM2.title ])
disp(['Contrast 1: ' SPM1.imnames(1,:)])
disp(['Contrast 2: ' SPM2.imnames(1,:)])
disp(['_____________________________________________________________'])
disp([SPM1.title ': ' num2str(length(clusters1)) ' clusters.'])
for i = 1:length(clusters1)
	disp(['	Cluster ' num2str(i) ': ' num2str(size(clusters1(i).XYZ,2)) ' total / ' num2str(size(clusters1(i).con1_XYZ,2)) ' ' SPM1.title ' / ' num2str(size(clusters1(i).con2_XYZ,2)) ' ' SPM2.title ' / ' num2str(size(clusters1(i).fuzzy,2)) ' undetermined / cluster lev. t = ' num2str(clusters1(i).cl_t) ', p = ' num2str(clusters1(i).cl_p)])
end

disp(['_____________________________________________________________'])
disp([SPM2.title ': ' num2str(length(clusters2)) ' clusters.'])
for i = 1:length(clusters2)
	disp(['	Cluster ' num2str(i) ': ' num2str(size(clusters2(i).XYZ,2)) ' total / ' num2str(size(clusters2(i).con1_XYZ,2)) ' ' SPM1.title ' / ' num2str(size(clusters2(i).con2_XYZ,2)) ' ' SPM2.title ' / ' num2str(size(clusters2(i).fuzzy,2)) ' undetermined / cluster lev. t = ' num2str(clusters2(i).cl_t) ', p = ' num2str(clusters2(i).cl_p)])
end


diary off

eval(['save ' myTit])

tor_image_compare_clusters(clusters1,clusters2) 


catch
	disp('Problem in tor_compare_clusters')
	lasterr
end

return





        