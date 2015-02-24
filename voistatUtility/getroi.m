% written by Tor Wager
% last modified 5/21/01 by Tor Wager



% sample input:
%
% ====== for montage =======
%anatname = 's8t1'; %'t1' %'meanavol_e3424_11_16_100_0001';
%orient = 'ax';
%range = 'all';
%tmaskname = 'spmT_0003'; %'ravol_e3424_11_16_100_0008';
%ROI.t_thresh = 2;
%tcolor = 'rs';
% ====== for getvoxels ========
%numareas = 2;
%usemask = 0;
% ====== for timeseries =======
%ROI.imgname = 'ravol_e3424_11_16_100_';  
%nimages = 200;
% ====== for window plot ======
%avgover = 10;
% ====== for voxel surfer =====
%meanf = 'meanavol_e3424_11_16_100_0001';
%voxsize = [3.12 3.12 4];
%HChoice = 40;
%TR = 2;

% USER ENTERS BEFORE RUNNING
% workspace variable funct:
%	'setupinput'
%	'newmontage'
%	'newplot'
%	'getregion'
%	'plotregion'
%	'initwindows'
%	'surf'
%	'roistats'
%	'voxstats'
%	'all'

if ~(exist('funct') == 1)
	disp(	'	setupinput')
	disp(	'	newmontage')
	disp(	'	newplot')
	disp(	'	getregion')
	disp(	'	plotregion')
	disp(	'	initwindows')
	disp(	'	surf')
	disp(	'	all')
	funct = inputdlg('Enter the function you want to run: ')
elseif isempty(funct)
	disp(	'	setupinput')
	disp(	'	newmontage')
	disp(	'	newplot')
	disp(	'	getregion')
	disp(	'	plotregion')
	disp(	'	initwindows')
	disp(	'	surf')
	disp(	'	all')
	funct = inputdlg('Enter the function you want to run: ')
end


	
switch funct
   
   
% ========================================================================================
% * newmontage *
% ========================================================================================

case 'newmontage'
% new montage
% --------------------------------------------------------------------------
if isempty(anatname),
   [filename,pathname] = uigetfile('*.img','pick an ANATOMY img file.');
   anatname = [pathname filename(1:end-4)];
   H = findobj('Tag','E1'); set(H,'String',anatname);
end

if isempty(tmaskname),
   [filename,pathname] = uigetfile('*.img','pick an OVERLAY img file.');
   tmaskname = [pathname filename(1:end-4)];
   H = findobj('Tag','E4'); set(H,'String',tmaskname);
end

T = montage4(anatname,orient,range,tmaskname,ROI.t_thresh,tcolor);






% ========================================================================================
% * newplot *
% ========================================================================================

case 'newplot'
%new plot
% --------------------------------------------------------------------------
T.range = range;
T.color{2} = tcolor;
redT = montage4(T,'mask');





% ========================================================================================
% * initwindows *
% ========================================================================================

case 'initwindows'
% init windows and region
% --------------------------------------------------------------------------
clear voifigH
clear vox
vox = voxelsurf(meanf,ROI.imgname,nimages,voxsize,avgover,HChoice,TR);





% ========================================================================================
% * getregion *
% ========================================================================================

case 'getregion'
% get voxels
% --------------------------------------------------------------------------
if usemask
	[ROI.voxels,mask,h] = getvoxels('poly',numareas,redT,2);
else
	[ROI.voxels,mask,h] = getvoxels('poly',numareas,redT);
end

ROI.ts = timeseries('multi',ROI.imgname,nimages,ROI.voxels);
%ROI.dstats = descriptives(ROI.ts.indiv);





% ========================================================================================
% * gettsfrommask *
% ========================================================================================

case 'gettsfrommask'
% use existing mask to get ROI data; not in The Works
% --------------------------------------------------------------------------
ROI.ts = timeseries('multi',ROI.imgname,nimages,ROI.voxels);
%ROI.dstats = descriptives(ROI.ts.indiv);
funct = 'plotregion';getroi





% ========================================================================================
% * plotregion *
% ========================================================================================
case 'plotregion'

% plot ROI data
% --------------------------------------------------------------------------
if ~exist('voifigH'),voifigH = [];,end

if (ishandle(voifigH) & ~useselective & ~usecontinuous)
	[voifigH,ROI.avgdata,ROI.adjustedy] = voistat('voifig',voifigH,ROI.ts.avg,[],avgover,HChoice,TR,numscans);
elseif (~useselective & ~usecontinuous)
	[voifigH,ROI.avgdata,ROI.adjustedy] = voistat('voifig',[],ROI.ts.avg,[],avgover,HChoice,TR,numscans);
elseif ishandle(voifigH)
	ROI.adjustedy = voistat('adjusty',ROI.ts.avg,HChoice,TR,numscans,usercovariates);
	[voifigH] = voistat('voifig',voifigH,ROI.adjustedy);
else
	ROI.adjustedy = voistat('adjusty',ROI.ts.avg,HChoice,TR,numscans,usercovariates);
	[voifigH] = voistat('voifig',[],ROI.adjustedy);
end
axes(voifigH(2)); title('Filtered average timeseries for ROI')

if useselective   % do selective averaging with events as cond function.
	if isempty(events),warning('No condition function found'),end
	[ROI.avgdata,ROI.avgste,ROI.ntrials] = voistat('trialavg',ROI.adjustedy,avgover,scolors,voifigH,vnl_condf,groups,basepoints,dosteplot,outthresh);
end

if usecontinuous
   disp('Getting continuous avg.')
   axes(voifigH(4));cla
   if isempty(events),error('events cell array is empty. specify event times.'),end
   [ROI.avg.time,ROI.avg.resp,ROI.avg.avg,ROI.avg.ste,ROI.avg.bincenter] = ...
      continuousavg(ROI.adjustedy,window,subjevents,groups,binres,'voifig',voifigH,'nopoints');
end

%[TMAP,hdr] = readim2(tmaskname,'t');
%checkvoxel(ROI.voxels,redT,TMAP,ROI.ts.indiv(:,1),ROI,[]);






% ========================================================================================
% * surf *
% ========================================================================================

case 'surf'
% surf a voxel
% --------------------------------------------------------------------------
%if isobject(vox.pointh),delete(vox.pointh),end
vox = voxelsurf(redT,vox);
axes(vox.voifigh(2)); title(['Filtered timeseries for voxel [' num2str(vox.coord) ']'])
[TMAP,hdr] = readim2(tmaskname,'t');
vox.z = TMAP(vox.coord(2), vox.coord(1), vox.coord(3));
checkvoxel(vox.coord,redT,TMAP,vox.ts,ROI,vox.pointh);
vox





% ========================================================================================
% * roistats *
% ========================================================================================

case 'roistats'
   % stats on the ROI ts
   if ~exist('X'),error('define X matrix prior to stats.'),
   elseif isempty(X) | ~size(size(X),2),error('X matrix is empty or not design matrix.')
   end 
   %ROI.gls = genls(X,ROI.adjustedy);
   ROI.gls = genls_ridge(X,ROI.adjustedy,.01,.05);

   format long
  	disp('T-values for X:')
	disp([num2str(ROI.gls.tbetas)])
	disp('P-values for X:')
	disp([num2str(ROI.gls.pbetas)])
   
   
   
   
   
% ========================================================================================
% * voxstats *
% ========================================================================================  
   
case 'voxstats'
   % stats on the ROI ts
   if ~exist('X'),error('define X matrix prior to stats.'),
   elseif isempty(X) | ~size(size(X),2),error('X matrix is empty or not design matrix.')
   end
   vox.gls = genls(X,vox.ts);
   format long
  	disp('T-values for X:')
	disp([num2str(vox.gls.tbetas)])
	disp('P-values for X:')
	disp([num2str(vox.gls.pbetas)])
	
   
   
   
   
   
% ========================================================================================
% * all *
% ========================================================================================   
   
case 'all'
	funct = 'setupinput'; getroi;
	funct = 'newmontage'; getroi;
	funct = 'newplot'; getroi;
	funct = 'initwindows'; getroi;
	funct = 'getregion'; getroi;
   	funct = 'plotregion'; getroi;
   	%funct = 'roistats'; getroi;
	%funct = 'surf'; getroi;
   
   
   
   
   
% ========================================================================================
% * batch *
% ========================================================================================
   
case 'batch'
   disp('===== starting batch processing ======')
   funct = 'setupinput';getroi;funct = 'batch';
   numsubs = size(imgdirarray,2);
   disp(['	... ' num2str(numsubs) ' subjects'])
   colors = {'rs-' 'b^-' 'go-' 'mv-' 'kx' 'c+'};
   batchfh = figure;
   
   for i = 1:numsubs
      disp(['	...Starting subject ' num2str(i)])
      
      % ---- get the timeseries ----
      ROI.batch.imnames{i} = [imgdirarray{i} '/' origimgname];
      ROI.batch.masknames{i} = [resdirarray{i} '/' origmaskname];
      ROI.batch.ts{i} = timeseries('multi',ROI.batch.imnames{i},nimages,ROI.voxels);
      
      % ---- get and plot the average and adjusted y ----
      if exist('voifigH')
			[voifigH,ROI.batch.trialavg1{i},ROI.batch.ts{i}.adjustedy] = voistat('voifig',voifigH,ROI.batch.ts{i}.avg,[],avgover,HChoice,TR);
		else
			[voifigH,ROI.batch.trialavg1{i},ROI.batch.ts{i}.adjustedy] = voistat('voifig',[],ROI.batch.ts{i}.avg,[],avgover,HChoice,TR);
		end
		axes(voifigH(2)); title('Filtered average timeseries for ROI')
      
      disp(['		...Getting trial average.'])
      % ---- get the continuous trial average, if possible ----
      if usecontinuous
         axes(voifigH(4));cla
   			if isempty(events),error('events cell array is empty. specify event times.'),end
   			[dummy,dummy2,ROI.batch.trialavg2{i}.avg,ROI.batch.trialavg2{i}.ste,ROI.batch.trialavg2{i}.bincenter] = ...
      		continuousavg(ROI.batch.ts{i}.adjustedy,window,events{i},groups,binres,'voifig',voifigH,'nopoints');
         
         % ---- plot the continuous timeseries on multiplot ----
         figure(batchfh)
      		subplot(4,4,i),hold on,
      		for j = 1:size(ROI.batch.trialavg2{i}.avg,2)
         		plot(ROI.batch.trialavg2{i}.avg{j},colors{j})
      		end
      
       else
          % ---- ...or plot the scan timeseries on multiplot ----

          figure(batchfh)      		
          subplot(4,4,i),hold on,
      		for j = 1:size(ROI.batch.trialavg1{i}.avg,2)
         		plot(ROI.batch.trialavg1{i}.avg{j},colors{j})
      		end
       end 
    end
    
   disp('Done TS for all subjects.')
   disp('Checking voxels for each subject')
   for i = 1:numsubs
      %P = getfiles(ROI.batch.masknames{i});
		[TMAP,hdr] = readim2(tmaskname,'t');
		checkvoxel(ROI.voxels,redT,TMAP,ROI.batch.ts{i}.indiv(:,1),ROI,[],ROI.batch.imnames{i});
   	end
   
   
   
      
      
      
% ========================================================================================
% * setup input *
% ========================================================================================
	
case 'setupinput'
clear H
H = findobj('Tag','E1');
if H

% ====== for montage =======
H = findobj('Tag','E1');
anatname = get(H,'String');

H = findobj('Tag','E2');
orient = get(H,'String');

H = findobj('Tag','E3');
range = get(H,'String');
if ~strcmp(range,'all'),range = str2num(range);,end

H = findobj('Tag','E4');
tmaskname = get(H,'String');

H = findobj('Tag','E5');
ROI.t_thresh = get(H,'String');
ROI.t_thresh = str2num(ROI.t_thresh);

H = findobj('Tag','E6');
tcolor = get(H,'String');

% ====== for getvoxels ========
H = findobj('Tag','E7');
numareas = get(H,'String');
numareas = str2num(numareas);

H = findobj('Tag','E8');
usemask = get(H,'String');
usemask = str2num(usemask);

% ====== for timeseries =======
H = findobj('Tag','E9');
ROI.imgname = get(H,'String');

H = findobj('Tag','E10');
nimages = get(H,'String');
nimages = str2num(nimages);

% ====== for window plot ======
H = findobj('Tag','E11');
avgover = get(H,'String');
avgover = str2num(avgover);

% ====== for voxel surfer =====
H = findobj('Tag','E12');
meanf = get(H,'String');

H = findobj('Tag','E13');
voxsize = get(H,'String');
voxsize = str2num(voxsize);

H = findobj('Tag','E14');
HChoice = get(H,'String');
HChoice = str2num(HChoice);

H = findobj('Tag','E15');
TR = get(H,'String');
TR = str2num(TR);

H = findobj('Tag','E16');
Xname = get(H,'String');

H = findobj('Tag','E17');
Xvar = get(H,'String');


% ====== for continuous average =====
H = findobj('Tag','ContinuousRadio');
usecontinuous = get(H,'Value');

H = findobj('Tag','E21');
binres = str2num(get(H,'String'));

H = findobj('Tag','E22');
window = str2num(get(H,'String'));

H = findobj('Tag','E23');
str = get(H,'String');
try
	eval(['events = ' str ';'])
catch
	disp(['loading event times / cond function from file.'])
	disp(['		(var name and file name should be the same for this to work.)'])
	eval(['load ' str])
	eval(['events = ' str ';'])
end

H = findobj('Tag','E24');
groups = str2num(get(H,'String'));


% ====== for selective average =====
H = findobj('Tag','SelectiveRadio');
useselective = get(H,'Value');

H = findobj('Tag','E27');
str = get(H,'String');
scolors = eval(str);

H = findobj('Tag','E28');
path2usercov = get(H,'String');
if ~(strcmp(path2usercov(end),filesep)),path2usercov = [path2usercov filesep];,end

H = findobj('Tag','E29');
usercovname = get(H,'String');

try
	if ~isempty(usercovname)
		eval(['load ' path2usercov usercovname ' -ASCII'])
		eval(['usercovariates = ' usercovname ';'])
	end
catch	
	disp('Warning: user covariate or mvmt parameter not found.')
	disp(['Command: load ' path2usercov usercovname])
end

H = findobj('Tag','E30');
dosteplot = get(H,'Value');

H = findobj('Tag','E31');
str = get(H,'String');
numscans = str2num(str);

H = findobj('Tag','E32');
str = get(H,'String');
basepoints = str2num(str);

H = findobj('Tag','E33');
str = get(H,'String');
outthresh = str2num(str);



% ====== for batch mode and directory specification =====
H = findobj('Tag','E18');
str = get(H,'String');
if ~(strcmp(str(end),filesep)),str = [str filesep];,end
eval(['imgdirarray{1} = ''' str ''';'])


H = findobj('Tag','E19');
str = get(H,'String');
if ~(strcmp(str(end),filesep)),str = [str filesep];,end
eval(['resdirarray{1} = ''' str ''';'])

H = findobj('Tag','E20');
str = get(H,'String');
if ~(strcmp(str(end),filesep)),str = [str filesep];,end
eval(['anatdirarray{1} = ''' str ''';'])


H = findobj('Tag','E25');
str = get(H,'String');
if ~isempty(str),
   	disp('Batch entry overrides entry for image directory.  Using batch.')
	if ~(strcmp(str(end),filesep)),str = [str filesep];,end
      	eval(['imgdirarray = ' str ';'])
      	if isempty(imgdirarray), disp('Warning: image directory list is empty.'),end
end

H = findobj('Tag','E26');
str = get(H,'String');
if ~isempty(str),
   	disp('Batch entry overrides entry for results directory.  Using batch.')
	if ~(strcmp(str(end),filesep)),str = [str filesep];,end
      	eval(['resdirarray = ' str ';'])
      	if isempty(resdirarray), disp('Warning: results directory list is empty.'),end
end

if ~isempty('imgdirarray')
   origimgname = ROI.imgname;
   ROI.imgname = [imgdirarray{1} ROI.imgname];
   meanf = [imgdirarray{1} meanf];
end
if ~isempty('resdirarray')
   origmaskname = tmaskname;
   tmaskname = [resdirarray{1} tmaskname];
   if ~isempty(Xname),Xname = [resdirarray{1} Xname];,end
end
if ~isempty('anatdirarray')
	anatname = [anatdirarray{1} anatname];
end


% ====== load the design matrix, if any ======
if ~isempty(Xname)
	fid = fopen(Xname);
	if fid > -1
      fclose(fid)
      load(Xname);
      eval(['X = ' Xvar ';'])
		disp('design matrix loaded from file.')
	else disp('cannot find X matrix file by the name specified.')
	end
end

% set to events{1}, unless running through batch option.
% =======================================================
if usecontinuous, subjevents = events{1};, end
      


else
   warning('GUI not found - Using default values for input!')	
% defaults
% ====== for montage =======
%anatname = 's8t1'; %'t1' %'meanavol_e3424_11_16_100_0001';
orient = 'ax';
range = 'all';
%tmaskname = 'spmT_0003'; %'ravol_e3424_11_16_100_0008';
ROI.t_thresh = 4;
tcolor = 'rs';
% ====== for getvoxels ========
numareas = 2;
usemask = 0;
% ====== for timeseries =======
ROI.imgname = 'ravol_e3424_11_16_100_';  
nimages = 200;
% ====== for window plot ======
avgover = 10;
% ====== for voxel surfer =====
meanf = 'meanavol_e3424_11_16_100_0001';
voxsize = [3.12 3.12 4];
HChoice = 40;
TR = 2;

end

end % end switch

	
