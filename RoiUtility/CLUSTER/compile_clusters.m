function [cl_summary,clusters,beta_clusters] = compile_clusters(O)
% cl_summary = compile_clusters(O)
% 
% by Tor Wager
%
% Use this function to compile sub*_clusters (timeseries data) 
% extracted with batch_multisubject_cluster_ui
%
% DX: deconvolution matrix
% number of columns = # time pts. estimated * # of trial types to estimate + 1 (for intercept)
%
% if you use batch_multisubject_cluster_ui to make original DX's (saved in sub#_clusters.mat)
% then each run is estimated separately, so number of trial types = trial types * # of runs
% e.g., for hi switch and lo switch, with 6 runs, # trial types = 12, and DX should have 12 * 16 + 1 columns
% if estimating 16 timepoints.

fprintf(1,'\n*compile_clusters.m\n________________________________________________________________\n')

% try to load saved betas, if we want them
% -----------------------------------------------------------------------
if O.extractbetas
        try
            eval(['load .' filesep 'beta_clusters'])
            disp(['Loaded beta_clusters.mat from current directory.'])
        catch
            disp('Did not find beta_clusters.mat in pwd: extracting betas from beta*img files')
            beta_clusters = [];
        end
end
            
pwd
fprintf(1,'Subject ')

%11 subjects in this expt
for mySub = O.numsubs

	fprintf(1,'%d ',mySub)

	eval(['load sub' num2str(mySub) '_clusters;'])
    
    % load the correct xX by loading spm.mat
    str = (['load ' O.an_files{mySub}]);
    %str = (['load c:\tor_documents\currentexperiments\biman5\' O.an_files{mySub}(end-24:end)]);
    %str(str == '/') = '\';
    
    try
        eval(str)
    catch
        try
            disp(['Cannot find: ' O.an_files{mySub}])
            str = (['load ' O.expdir filesep O.an_files{mySub}(end-O.numcharsback:end)]);
            str(str == '/') = '\';
            disp(['Trying: ' str]) 
            eval(str)
        catch
            
            disp(['looking for: ' O.an_files{mySub}])
            disp('Cannot find SPM.mat file for this subject.  Please locate.')
            P = spm_get(1,'*',['Choose SPM.mat or an def for sub' num2str(mySub)]);
            eval(['load ' P])
        end
        
    end
    
	if O.radius > 0
		% break into sub-clusters and re-save 
		% --------------------------------------------------------
		[clusters,orig_clusters,sphere,SPM,VOL] = tor_get_spheres2(clusters,O.radius,SPM,VOL);

		%xX.RT = O.TR;
		%[clusters,HP] = analyze_cluster_rois(clusters,xX,O.HP);

		eval(['save sub' num2str(mySub) '_subclusters;'])
	end



	if O.plot
		figure;
		axisx = ceil(sqrt(length(clusters)));
		axisy = floor(sqrt(length(clusters)));
		while axisx * axisy < length(clusters)
			axisy = axisy + 1;
		end
	end

    
    % hack: baseline time shift, to avoid averaging at beginning of baseline period
    % -----------------------------------------------------------------------
    if isfield(O,'basegroup') & isfield(O,'baseshift')
        disp(['Shifting over DX functions for baseline group.'])
        if length(O.groups) ~= (length(O.adjustfor) ./ O.nruns)
            error('O.groups must be same length as number of within session regressors.')
        end
        
        myloc = repmat(O.groups,1,O.nruns);
        myloc = find(myloc == O.basegroup);     % get which columns to shift
        myz = zeros(O.baseshift,1);             % zeros to append before column
        for mycol = myloc
            DX(:,mycol) = [myz; DX(1:end-length(myz),mycol)];
        end
    end
        
        

	% create custom adjustment matrix from subject's xX
    % so that these effects are filtered out with filterAdjust
	% -----------------------------------------------------------------------
    if any(O.adjustfor) & O.regressOut
        O.adjustmatrix = xX.xKXs.X(:,find(O.adjustfor));    % smoothed regressors
    end
    
	% collapse sf into specified groups
	% and re-make DX to estimate an arbitrary number of time points
	% -----------------------------------------------------------------------
	sf = DX(:,1:O.origtp:size(DX,2) - 1);
	for j = 1:max(O.groups)
		sf2(:,j) = sum(sf(:,O.groups == j),2);
	end
	% original DX matrix with O.tp timepoints, for making myDelta, etc.
	DX = tor_make_deconv_mtx3(sf,O.tp,1);
	% DX collapsed into groups.
	DX3 = tor_make_deconv_mtx3(sf2,O.tp,1,O.tbefore);

	% make DX2, which collapses all conditions
	% -----------------------------------------------------------------------
	%sf = sum(DX(:,1:O.tp:size(DX,2) - 1),2);
	sf = sum(sf,2);
	DX2 = tor_make_deconv_mtx3(sf,O.tp,1,O.tbefore);

    
	% get title for figures, if necessary
    % -----------------------------------------------------------------------
    %if ~isfield(clusters,'title'),	
	    for j = 1:length(clusters),
		    clusters(j).title = O.mytitle;
	    end
        %end

	% eliminate clusters of size 1 voxel; they cause problems
    % -----------------------------------------------------------------------
    elim = cat(1,clusters.numVox);
    if any(elim == 1), disp(['Warning!  Omitted ' num2str(sum(elim == 1)) ' clusters that had voxel sizes of 1.'])
        disp(['These clusters were: ' num2str(find(elim == 1)') '; cl. numbers have been adjusted accordingly.'])
    end
    clusters(elim == 1) = [];
    
    
    
    % special: make CLU from clusters and VOL.M (needs SPM.mat loaded)
    % xM.VM.mat contains affine transformation mat - as in VOL.M
    % don't need CLU any more
    %if ~exist('CLU') == 1
    %    CLU = clusters2CLU(clusters,xM.VM.mat);
    %end
        
    % extract betas from beta*img, if we have the CLU structure
    % -----------------------------------------------------------------------
    if O.extractbetas & length(beta_clusters) < length(O.numsubs)

                    % get image names for all betas from current subject
                    str = (['load ' O.expdir filesep O.an_files{mySub}(end-O.numcharsback:end)]);
                    str(str == '/') = '\';
                    subDir = str(6:end-7);
                    imloc = [subDir 'beta*img'];
                    try
                        imnames_2 = getfiles2(imloc);
                        imnames = [imnames_2(end,:); imnames_2(1:end-1,:)];
                    catch
                        disp(['Looking for: ' O.expdir filesep O.an_files{mySub}(end-24:end)]);
                        subDir = spm_get(-1,'*',['Choose sub directory for ' num2str(mySub)]);
                        imloc = [subDir filesep 'beta*img'];
                        imnames_2 = getfiles2(imloc);
                        imnames = [imnames_2(end,:); imnames_2(1:end-1,:)];
                    end
                    
                    disp(['Getting beta values from beta*imgs'])
                    % extract betas from images
                    beta_clusters{mySub} = tor_extract_rois_cls(imnames,clusters);
                    
    end
     
    
    % get all selective avgs and deconvolution estimates for each cluster
    % -----------------------------------------------------------------------
    
	for j = 1:length(clusters)
        
		if O.plot,subplot(axisy,axisx,j),end

		% ----------------------------------------------------------------
		% get deconv hrfs for this cluster
		% ----------------------------------------------------------------
		O.y = clusters(j).timeseries;
 		[y,O] = filterAdjust(O);
		clusters(j).adjustedy = y;

        %figure;hold on; subplot(2,1,1)
        %plot(clusters(j).timeseries,'b','LineWidth',1.5)
        %subplot(2,1,2)
        % plot(clusters(j).adjustedy,'r','LineWidth',1.5)
        % pause(4)
        % close
         
		clusters(j).DXb = pinv(DX3) * clusters(j).adjustedy; 
		clusters(j).oDXb = pinv(DX2) * clusters(j).adjustedy; 

		index = 1;
		for k = 1:O.tbefore+O.tp:length(clusters(j).DXb) - 1
			cl_summary(j).hrf_ind{index}(mySub,:) = clusters(j).DXb(k:k+O.tbefore+O.tp-1);
			index = index + 1;
		end

        % ----------------------------------------------------------------
		% if getting regression betas from O.adjustfor, get them
        % they are stored in O.custombeta as a column vector
		% ----------------------------------------------------------------
		if any(O.adjustOfInterest)
            if isfield(O,'custombeta')
                cbeta = O.custombeta(find(O.adjustOfInterest))';
                % this was wrong in previous versions - but should be ok of same # of runs and vars of interest
                clusters(j).custombeta = mean(reshape(cbeta,length(cbeta)/O.nruns,O.nruns)');    
                % saves the mean of each beta of interest across scans
                cl_summary(j).custombeta(mySub,:) = clusters(j).custombeta;
            end
            
            % converting to % in filterAdjust changes the betas in weird ways; avoid...???
            % tested 3/25/02, betas w/o percent conversion are almost exactly betas from beta*img files.
            O2 = O;
            if isfield(O2,'adjustmatrix'),O2 = rmfield(O2,'adjustmatrix');,end
            O2.percent = 0; O2.y = clusters(j).timeseries;
            newy = filterAdjust(O2);
            cbeta2 = pinv(xX.xKXs.X) * newy;
            cbeta2 = cbeta2(find(O.adjustfor));
            clusters(j).custombeta2 = mean(reshape(cbeta2,length(cbeta2)/O.nruns,O.nruns)'); 
            cl_summary(j).custombeta2(mySub,:) = clusters(j).custombeta2;
            
            % betas for 1st voxel only, for check
            O2.y = clusters(j).all_data(:,1);
            newy = filterAdjust(O2);
            cbeta4 = pinv(xX.xKXs.X) * newy;
            cbeta4 = cbeta4(find(O.adjustfor));
            clusters(j).custombeta4 = mean(reshape(cbeta4,length(cbeta4)/O.nruns,O.nruns)'); 
            cl_summary(j).custombeta4(mySub,:) = clusters(j).custombeta4;
            
            % now get the actual betas from the image files
            if exist('beta_clusters') == 1
                    
                    % average betas of interest
                    cbeta3 = beta_clusters{mySub}(j).timeseries(find(O.adjustfor));
                    clusters(j).custombeta3 = mean(reshape(cbeta3,length(cbeta3)/O.nruns,O.nruns)'); 
                    cl_summary(j).custombeta3(mySub,:) = clusters(j).custombeta3; 
                    
                    myc = corrcoef(cbeta2,cbeta3);
                    myc2 = corrcoef(cbeta4,beta_clusters{mySub}(j).all_data(find(O.adjustfor),1));
                    disp([O.mytitle ' cl ' num2str(j) ' fitted w/extracted betas: r = ' num2str(myc(1,2)) ' average, r = ' num2str(myc2(1,2)) ' 1st voxel.']) 
                     
            end
            
        end
        
		% ----------------------------------------------------------------
		% plot selective averages for individual
		% ----------------------------------------------------------------
		if j == 1
			O.title = ['Sub ' num2str(mySub) ' contrast ' clusters(j).title(end-3:end) ' cluster ' num2str(j) ' - tmax = ' num2str(max(clusters(j).Z))]; 
		else
			O.title = ['cluster' num2str(j)];
		end		

		myDelta = DX(:,[1:O.tp:size(DX,2)-1]);
		
		[ROI,clusters(j).sel_avg_all_conds] = trialavg2(y,myDelta,[O.mytrialstart O.tp-1],'options',O);
		clusters(j).sel_avg_ind = ROI.grpavg;		

		if O.plot
			legend off
			set(gcf,'Position',[ 224    75   870   871])
			set(gcf,'Color','w')
		end

		% ----------------------------------------------------------------
		% save group selective averages and deconv betas
		% ----------------------------------------------------------------

		for myCond = 1:length(clusters(j).sel_avg_ind)
			cl_summary(j).sel_avg_ind{myCond}(mySub,:) = clusters(j).sel_avg_ind{myCond}(1,:);

			cl_summary(j).DXb_ind(mySub,:) = clusters(j).DXb'; 
			cl_summary(j).oDXb_ind(mySub,:) = clusters(j).oDXb';
			
		end
		
		if O.plot,drawnow,end

	end

	%if O.plot,close,end
    
end % loop through subjects

% ----------------------------------------------------------------
% * change zeros (empty rows) to NaNs
% ----------------------------------------------------------------
for j = 1:length(cl_summary)
    
    for k = 1:length(cl_summary(j).hrf_ind)
            cl_summary(j).hrf_ind{k}(find(sum(cl_summary(j).hrf_ind{k} == 0,2)),:) = NaN;
            cl_summary(j).sel_avg_ind{k}(find(sum(cl_summary(j).sel_avg_ind{k} == 0,2)),:) = NaN;
        end
    
    myfields = {'custombeta2' 'custombeta3' 'custombeta4' 'DXb_ind' 'oDXb_ind'};
    if isfield(cl_summary,'custombeta'), myfields{end+1} = 'custombeta';, end
    
    for k = myfields
        if isfield(cl_summary,k{1})
                eval(['cl_summary(j). ' k{1} '(find(sum(cl_summary(j).' k{1} ' == 0,2)),:) = NaN;'])
                %cl_summary(j).DXb_ind(find(sum(cl_summary(j).DXb_ind == 0,2)),:) = NaN;
        end
    end
end
    
    
% ----------------------------------------------------------------
% * get avg deconv betas and custom betas
% ----------------------------------------------------------------
for j = 1:length(cl_summary)    
	cl_summary(j).DXb_avg = nanmean(cl_summary(j).DXb_ind);
	cl_summary(j).oDXb_avg = nanmean(cl_summary(j).oDXb_ind);
    cl_summary(j).DXb_ste = nanstd(cl_summary(j).DXb_ind) ./ sum(~isnan(cl_summary(j).DXb_ind(:,1)));
    
    if any(O.adjustOfInterest)
        if isfield(cl_summary,'custombeta')
            cl_summary(j).custombeta_avg = nanmean(cl_summary(j).custombeta);
            cl_summary(j).custombeta_ste = nanstd(cl_summary(j).custombeta) ./ sum(~isnan(cl_summary(j).custombeta(:,1)));
        end
        
        cl_summary(j).custombeta_avg2 = nanmean(cl_summary(j).custombeta2);
        cl_summary(j).custombeta_ste2 = nanstd(cl_summary(j).custombeta2) ./ sum(~isnan(cl_summary(j).custombeta2(:,1)));
        
        cl_summary(j).custombeta_avg4 = nanmean(cl_summary(j).custombeta4);
        cl_summary(j).custombeta_ste4 = nanstd(cl_summary(j).custombeta4) ./ sum(~isnan(cl_summary(j).custombeta4(:,1)));
        
        if exist('beta_clusters') == 1
            cl_summary(j).custombeta_avg3 = nanmean(cl_summary(j).custombeta3);
            cl_summary(j).custombeta_ste3 = nanstd(cl_summary(j).custombeta3) ./ sum(~isnan(cl_summary(j).custombeta3(:,1)));  
        end
    end
end

% ----------------------------------------------------------------
% get selective averages
% ----------------------------------------------------------------

for k = 1:length(cl_summary(1).sel_avg_ind)

	for j = 1:length(clusters)

		cl_summary(j).sel_avg_ind{k}( cl_summary(j).sel_avg_ind{k} == 0) = NaN;

		cl_summary(j).sel_avg{k} = nanmean(cl_summary(j).sel_avg_ind{k});
		cl_summary(j).sel_avg_ste{k} = nanstd(cl_summary(j).sel_avg_ind{k}) ./ sum(~isnan(cl_summary(j).sel_avg_ind{k}(:,1)));
		
	end

end


if ~(exist('beta_clusters') == 1)
    beta_clusters = [];
end


return