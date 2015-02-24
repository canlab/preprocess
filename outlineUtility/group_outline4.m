function [gvol,clusters] = group_outline3(mysbs,u,k,wcon,varargin)
% function [gvol,clusters] = group_outline3(mysbs,u,k,wcon,[opt] bkground img name, cell array of colors)
%
% You must start in the directory above the individual subject results directories!
% Background overlay image must be in same space as individual subjects!
% Background image is hard-coded in this script!
%
% mysbs		string of individual subject directory codes
% u		height threshold
%		corr. method is hard-coded in ind_subject_outline
% k		extent threshold
% wcon		index number of which contrast to show in results
%
% gt is the threshold for how many subjects must activate for the voxel to be plotted
%	if gt is a vector, plot will be separated in colors based on multiple levels
%	order: 'g' 'b' 'r'.  maximum of 3 levels recommended but not necessary.
%	elements of gt should be in ASCENDING ORDER!
% 	THIS VERSION AUTOMATICALLY PICKS GT based on HISTOGRAM OF GVOL (range of subject counts
%	of significance in each voxel)
%
% Optional inputs:
% background image name; try using spm_get to get a filename of an analyze img file
% cell array of colors: e.g., {'g' 'b' 'r' 'y'}
%
% gvol		group mask, with value in elements indicating number of subjects 
%		activating that voxel.
%
% Tor Wager, 2/15/02    Improved version over group_outline, but may be quite a bit slower.
% This version does NOT plot individual outlines, as the previous group_outline does.
%
% group_outline4 does not rely on loading SPM results
%
% svol is the individual subject mask of results for each subject0
% gvol is the sum of individual subject masks
%
%
% Uses montage_clusters.m
% Example:
% Po = spm_get(1,'*img','Choose anatomical overlay image');
% [gvol,clusters] = group_outline3(EXPT.subjects,.05,5,5,Po);
% [gvol,clusters] = group_outline3(EXPT.subjects,.05,5,5,Po,{'b' 'r' 'y'});

docolbar = 1;

if length(varargin) > 0, 
    bgimg = varargin{1};, 
else
    bgimg = spm_get(1,'*img','Choose anatomical overlay image.'); % '/usr/private/spm99/templates/scalped_avg152T1.img';
end
bgV = spm_vol(bgimg);

if length(varargin) > 1,
    mycol = varargin{2};
else
    mycol = {'g' 'b' 'r' 'y' 'm' 'c' 'k'};
end

% ---------------------------------------------------------------
% * loop through subjects, get results, make outline, plot
% ---------------------------------------------------------------

index = 1;
for snum = mysbs
	
	SubjCode = snum{1};
  	if index == 1, cd([SubjCode]), else  cd(['../' SubjCode]), end

    
    
    % -------------------------------------------------------------
    % * process for individual subject
    % -------------------------------------------------------------

    % get the results for this subject
    % - - - - - - - - - - - - - - - - - - - - - - - - - -
	TOR.u = u;
	TOR.k = k;
	TOR.Ic = wcon;
	TOR.resdir = pwd;
	TOR.maskWOtherContrasts = 0;
	TOR.multCompCorrect = 'uncorrected'; %'FDR' 'FWE' 'uncorrected'; 
	[hReg,SPM,VOL,xX,xCon,xSDM] = tor_spm_results_ui(TOR);
    svol = voxel2mask(SPM.XYZ',VOL.DIM');
    mytit{index} = SPM.title;
     
    % check VOL mat files - all subs should be the same
    if index > 1,
        if any(VOL.M - oldVOL.M), 
            error([snum{1} ': Subject VOL.M field differs from previous subject!'])
        end
    end
    oldVOL = VOL;
    

	sxyz{index} = SPM.XYZ;
	sz{index} = SPM.Z;

	% add subject to group mask	
	% --------------------------------------------------------------
	if index == 1, gvol = svol;, else,  gvol = gvol + svol;, end
	svols{index} = svol;	% save for making table later
	index = index + 1;

end

cd ..

% -------------------------------------------------------------
% * check titles to make sure they're all the same.
% -------------------------------------------------------------
testdt = str2mat(mytit);
if length(strmatch(testdt(1,:),testdt)) < size(testdt,1)
    warning('Contrast titles do not all match!!!')
end


% -------------------------------------------------------------
% * table of number subjects by number of voxels
% -------------------------------------------------------------
if exist('ind_outlines','dir') == 0, mkdir ind_outlines, end
diary(['ind_outlines' filesep 'group_outline3_results.txt'])

fprintf(1,'\nTotal subjects = %3.0f, u = %3.2f, k = %3.0f, contrast = %3.0f\n',length(mysbs),u,k,wcon)
fprintf(1,'\nCorrection method: %s', TOR.multCompCorrect)
fprintf(1,'\nTotal voxels = %6.0f',sum(sum(sum(gvol(~isnan(gvol))))))
fprintf(1,'\nSubjects\tProp. Sub\tNo. voxels')
gogo = 1; gtmean = [];
for myg = 0:length(mysbs)
	mycount = sum(sum(sum(gvol >= myg)));
	fprintf(1,'\n%3.0f\t%3.2f%%\t%5.0f\t',myg,myg/length(mysbs),mycount)
	if mycount < 3000 & gogo, gtmean = myg;, gogo = 0;, end
end
fprintf(1,'\n')
if isempty(gtmean), disp(['Sig. voxel counts above 3000 for all subjects! Exiting...']), return, end

% -------------------------------------------------------------
% * find gt (thresholds) based on the number of found subjects
% -------------------------------------------------------------
gtmax = max(max(max(gvol))); disp(['Max subjects activating any voxel is: ' num2str(gtmax)])
if gtmax == 0, disp(['No results: Exiting']),return, end
% gtmean = max(round(mean(mean(mean(gvol(gvol>0))))),round(length(mysbs)*.4));
% define based on how many voxels are found - see table code above.

if length(varargin) > 1,
	gt = input('Enter number of subjects for each threshold value (e.g., [3 5 7]): ');
else
	gt = [gtmean round(gtmean + ((gtmax-gtmean)./2)) gtmax];
	gt = unique(gt);
end

disp(['gt threshold values are: ' num2str(gt) ' = ' num2str(round(100*gt ./ length(mysbs))) '% of sample'])
if any(gt - sort(gt)), error('gt must be in ascending order.'),end


% -------------------------------------------------------------
% * table of individual subject voxels at each threshold gt level
% -------------------------------------------------------------
% header
fprintf(1,'\nTable of individual subject activation summaries')
fprintf(1,'\nSubject\tTotal Vox\t')
for j = 1:length(gt),
	fprintf(1,'Area def. by %3.0f subs\t\t',gt(j))
end
fprintf(1,'\n\t\t')
for j = 1:length(gt),
	fprintf(1,'Vox in area\tProp this sub\t')
end
fprintf(1,'\n')
% body
for i = 1:length(svols)
	subcoverall = sum(sum(sum(svols{i} > 0)));
	fprintf(1,'%s\t%3.0f\t',mysbs{i},subcoverall)
	for j = 1:length(gt)
		% find out how many voxels in group mask were activated by this subject
		mysubvol = sum(sum(sum(svols{i} .* gvol >= gt(j))));
		fprintf(1,'%3.0f\t%3.2f\t',mysubvol,mysubvol./sum(sum(sum(gvol >= gt(j)))));
	end
	fprintf(1,'\n')
end

diary off

% -------------------------------------------------------------
% * Make clusters out of mask, for compatibility with montage_clusters
% -------------------------------------------------------------
gt(end+1) = Inf;
for i = 1:length(gt)-1
    
    [x,y,z] = ind2sub(size(gvol),find(gvol >= gt(i) & gvol < gt(i+1)));
    XYZ = [x y z]';
    clusters{i}.XYZmm = voxel2mm(XYZ,VOL.M);
    clusters{i}.title = mytit{1};
    clusters{i}.threshold = gt(i);
    
end
  
% -------------------------------------------------------------
% * Write image file of gvol (borrowing hdr info from an existing con img)
% -------------------------------------------------------------

if i < 10, myz = '000';, elseif i < 100, myz = '00'; else, myz = '0';,end
V = spm_vol(fullfile(TOR.resdir,['con_' myz num2str(i) '.img']));
V.fname = ['ind_outlines' filesep 'group_con' num2str(i) '.img'];
V.descrip = ['u=' num2str(TOR.u) ' k=' num2str(TOR.k) ' corr=' TOR.multCompCorrect ' dir=' TOR.resdir];

spm_write_vol(V,gvol);

% -------------------------------------------------------------
% * Image montage with montage_clusters.m
% -------------------------------------------------------------

str = ['montage_clusters(bgimg,clusters{1}'];
for i = 2:length(gt)-1
    if ~isempty(clusters{i}.XYZmm)
    	str = [str ',clusters{' num2str(i) '}'];
    end
end
% 0 is for suppressing overlap plotting, mycol is color specification
str = [str ',0,mycol);'];

if ~isempty(clusters{1}.XYZmm)
	eval(str)
else
	disp(['Con ' num2str(wcon) ': No voxels significant for ' num2str(gt(1)) ' participants.'])

end

% -------------------------------------------------------------
% * Compare subjects on a slice
% -------------------------------------------------------------
v = spm_read_vols(bgV);
bgV.M = bgV.mat;


ss = 1;
while ~isempty(ss)
	ss = input(['Pick a slice (1:' num2str(size(gvol,3)) ')']);
	if isempty(ss), break, end

	figure('Color','w'); rc = sqrt(length(svols))+1;

	ssmm = voxel2mm([0 0 ss]',VOL.M);
	sso = mm2voxel(ssmm,bgV); sso = sso(3);

	subplot(rc,rc,1), imagesc(v(:,:,ss)'); set(gca,'YDir','normal');
        hold on; axis image; axis off;colormap gray

	% color map
	h1 = (0:1/99:1)';
	h2 = ones(size(h1)); 
	h3 = zeros(size(h1));
	h = [h1 h3 h3; h2 h1 h3; h2 h2 h1];
	h(1:50,:) = [];
	% in new matlab: h = colormap(hot(300));

	% determine overall z-score range
	zrange = cat(2,sz{:}); zrange = [min(zrange) max(zrange)];
	
	zh = zrange(1):(zrange(2)-zrange(1))./249:zrange(2);
	zh = round(zh*100);

	subplot(rc,rc,1);imagesc(v(:,:,sso)'); hold on; colormap gray
	set(gca,'YDir','normal'),axis image, axis off
	
	gslice = gvol(:,:,ss);
	[x,y] = ind2sub(size(gslice),find(gslice > 0)); 
	myz = gslice(find(gslice > 0));
	z = ss * ones(size(x));

	myxyzmm = voxel2mm([x y z]',VOL.M);
	xyz = mm2voxel(myxyzmm,bgV,1)';
	clear h2
	for i = 1:length(gt)-1
		myxyz = xyz(:,myz >= gt(i));
		for j = 1:length(myxyz)
			h2(j) = plot(myxyz(1,j),myxyz(2,j),'Color',mycol{i},'MarkerSize',2,'MarkerFaceColor',mycol{i});
		end
	end
	if exist('h2') == 1, set(h2,'Marker','square'),end	
	title('Group'),xlabel(['rgb = ' num2str(gt(1:end-1))])

	hold on
	for i = 1:length(sxyz)
		subplot(rc,rc,i+1), imagesc(v(:,:,sso)'); hold on; colormap gray
		set(gca,'YDir','normal'),axis image, axis off

		myxyz = sxyz{i}(:,sxyz{i}(3,:) == ss);
		myz = sz{i}(sxyz{i}(3,:) == ss);

		myxyzmm = voxel2mm(myxyz,VOL.M);
		myxyz = mm2voxel(myxyzmm,bgV,1)';
		clear h2
		for j = 1:size(myxyz,2)
			tmp = find(round(myz(j)*100) == zh);
			if isempty(tmp), 
				tmp = find((zh-round(myz(j)*100)).^2 == min((zh-round(myz(j)*100)).^2));
			end

			wh(j) = tmp(1);

			h2(j) = plot(myxyz(1,j),myxyz(2,j),'Color',h(wh(j),:),'MarkerSize',2,'MarkerFaceColor',h(wh(j),:));

		end
		if exist('h2') == 1, set(h2,'Marker','square'),end
		title(mysbs{i})
	drawnow
	end
	
    % -------------------------------------------------------------
	% color scale bar - we must create by hand
    % -------------------------------------------------------------
    
    % does not work well for Matlab 5.3!
    %cc = colormap(gray); cc(1:10,:) = repmat([1 1 1],10,1);    % [0 0 .3] for dark blue
    %colormap(cc)

    if docolbar % only on the 1st time thru
    % weird bug in this; make sep fig
	%subplot(rc,rc2,length(sxyz)+1),hold on
    figure('Color','w'); hold on;
	zh2 = zh./100;
    axis([0 .3 zh2(1) zh2(end)]),hold on
	for i = 1:size(h,1), plot([0 1],[zh2(i) zh2(i)],'Color',h(i,:));, end
	%set(gca,'YLim',[zh(1) zh(end)])
	set(gca,'XTickLabel',''); % ylabel('Z-score')
	%h3 = get(gca,'Position');
	%set(gca,'Position',[h3(1:2) h3(3)./3 h3(4)])
	h3 = get(gcf,'Position');
    set(gcf,'Position',[h3(1:2) h3(3)*.3 h3(4)*.5])
    docolbar = 0;
    
    end
	
end	% end loop


return
