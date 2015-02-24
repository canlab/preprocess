


	OPT.studydir = '/data/placebo';
		% path to your study main directory, no slash at the end

	OPT.nruns = 6;	
		% number of runs in your study
	

	% anatomy
	% ------------------------------------------------------------------------------------
	OPT.coreg_anat = 0;			% coreg structural object to structural overlay object or functional
	OPT.coreg_2_funct = 0;		% coreg object to a functional scan (1st funct) instead of overlay
							% must use getfunctnames = 1 also for this to work properly.

	OPT.target = 'het1overlay.img';
							% name of structural (target) image for coregistration
							% ignore if coreg to first funct - use getfunctnames option instead.
 

	OPT.segment_object = 0;		% segment the object image
	OPT.unseg_object = 'het1spgr.img';
							% name of unsegmented object image

	OPT.object = 'het1spgr.img';	
							% name of t1 image to coregister (object) 
							% and use to determine normalization parameters
							% to use gray matter only, use (e.g.) 'ht1spgr_seg1.img'

	OPT.canonicalT1 = '/usr/private/spm99/templates/sscalped_avg152T1.img';
							% name of canonical template to normalize to.
							% to use gray matter only, use (e.g.) 'T1_seg1.img'

	OPT.secondcanonicalT1 = '/usr/private/spm99/templates/scalped_avg152T1.img';
							% name of canonical template to normalize to.
							% only for twostagenorm: this is the 2nd (rougher) template

	% functionals
	% ------------------------------------------------------------------------------------
	OPT.createmean = 0;			% create mean functional image - spm's is shifted over for me.
	OPT.detNorm = 0;			% determine normalization parameters from anatomicals
	OPT.applyNorm = 0;			% normalize realigned functional volumes (ravol*img)
	OPT.twostagenorm = 0;		% two-stage normalization: choose smoothed template first, then rougher
							% re-applies normalization to nspgr, writes and smooths nnravols. 
							% detnorm and applynorm refer to the FIRST normalization.
							% twostagenorm runs BOTH det and apply for the SECOND normalization.
							% These two can be used independently.
	OPT.smooth = 0;				% smooth normalized volumes
	OPT.voxsize = [3.75 3.75 5];% voxel size in mm
	OPT.kernel = [9 9 9];		% smoothing kernel in mm
	OPT.check_coreg = 0;		% check coregistration
	OPT.handadjust = 0;			% only works with check_coreg.  Adjust image alignment manually.
	OPT.checknorm = 0;
	OPT.bigmask = 0;			% create bigmask image from individual subject anatomy
	OPT.getfunctnames = 1;		% get functional names; useful for normalizing and checking
	
	moviecheck = 0;				% check normalization or realignment by viewing movies
	plotparams = 0;				% plot and save jpegs of movement parameters
	printparams	= 0;			% print movement param figure to printer.

% ------------------------------------------------------------------------
% END USER INPUT
% ------------------------------------------------------------------------

if OPT.getfunctnames
	eval(['cd ' OPT.studydir])
        EXPT.nsess = OPT.nruns;
        clear scan1dir
        for i = 1:length(EXPT.subjects)        

            for j = 1:EXPT.nsess, scan1dir{i}{j} = [EXPT.subjects{i} filesep 'scan' num2str(j)];,end
            EXPT.im_files{i} = str2mat(tor_list_files(scan1dir{i},'snnra*.img'));
	    EXPT.nravols{i} = str2mat(tor_list_files(scan1dir{i},'nnra*.img'));
            
        end

        clear scan1dir
        for i = 1:length(EXPT.subjects) 
            scan1dir{i} = [EXPT.subjects{i} filesep 'scan1'];
            anatdir{i} = [EXPT.subjects{i} filesep 'anatomy'];
        end
            
        % save stuff in EXPT for later reference
        OPT.firstfun = str2mat(tor_list_files(scan1dir,'snnra*0001.img'));
        OPT.firstnra = str2mat(tor_list_files(scan1dir,'nnra*0001.img'));
        OPT.firstra = str2mat(tor_list_files(scan1dir,'ra*0001.img'));
       OPT.nspgr = str2mat(tor_list_files(anatdir,['n' OPT.object]));
    end

EXPT.SUBJECT = OPT;

if plotparams
	eval(['cd ' OPT.studydir])
	for mysub = EXPT.subjects
	   if isempty(mysub{1})
	   else

			eval(['cd ' mysub{1}])


			% check for number of scans found - and adjust if necessary!
			% ---------------------------------------------------------------------
			D = dir; 	% eval(['D=dir(''' OPT.subjcode ''');'])
			m = strmatch('scan',str2mat(D.name));
			if length(m) ~= OPT.nruns
				warning(['Only ' num2str(length(m)) ' scan directories found.  Adjusting OPT.nruns.'])
				OPT.nruns = length(m);
			end

			try
			um_plot_params

			% guidelines for moving spgr.img
			% -------------------------------------
			cd movement
			load real_params
			a = params(end,:);
			a(1:3) = a(1:3) .* pi ./ 180;
			fprintf(1,'\nMovement End\npitch\troll\tyaw\tx\ty\tz\n')
			fprintf(1,'%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t\n',a(1),a(2),a(3),a(4),a(5),a(6))
			gosave = 1;
			catch
				disp(['Problem with movement params for ' pwd '...skipping.'])
				gosave = 0;
			end

			cd ..

			cd ..
			
			if gosave
				saveas(gcf,[mysub{1} '_movement'],'jpg')
				if printparams
					print -dps2 -Pjjon-print1
				else
					pause(5)
				end
				close
			end
		

	   end
	end
end

	for mysub = EXPT.subjects
		if isempty(mysub{1})
		else
			eval(['cd ' OPT.studydir])
			if OPT.coreg_2_funct, OPT.target = deblank(OPT.firstra(strcmp(mysub{1},EXPT.subjects),:));, end
			OPT.subjcode = mysub{1};
			jlab_preproc(OPT)
		end
	end

if moviecheck
	um_make_movie(EXPT.nravols,EXPT.subjects);
end





