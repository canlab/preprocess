function tor_smooth(DIRS,imgwildcard,kernel)
% function tor_smooth(DIRS,imgwildcard,kernel)
%
% DIRS is cell array of image directories
% imgwildcard is something like ra*vol


% -------------------------------------------------------------------------
% * loop through directories, change to each one, get files, and smooth
% -------------------------------------------------------------------------

for i = 1:length(DIRS)

	go = 1;

	str = ['cd(''' DIRS{i} ''')'];
	disp(['Changing to directory with command: ' str])
	
	try
		eval(str)
		[Files,Dirs] = spm_list_files('.',imgwildcard);
		disp(['Your 1st file is: ' Files(1,:) ', ' num2str(size(Files,1)) 'images in queue.'])
		
	catch
		warning('Can''t find directory:')
		disp([str])
		disp('Skipping directory...')
		go = 0;
	end


	if go == 1,
           try

		% lifted almost directly from spm_smooth_ui
		% implement the convolution
		%---------------------------------------------------------------------------
		SPMid = spm('FnBanner',mfilename,'2.4');
		[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Smooth');
		spm_help('!ContextHelp','spm_smooth_ui.m');

		spm('Pointer','Watch');
		spm('FigName','Batch Smooth: working',Finter,CmdLine);
		spm_progress_bar('Init',size(Files,1),'Script Batch Smoothing','Volumes Complete');
		for i = 1:size(Files,1)
			Q = deblank(Files(i,:));
			[pth,nm,xt,vr] = fileparts(deblank(Q));
			U = fullfile(pth,['s' nm xt vr]);
			spm_smooth(Q,U,kernel);
			spm_progress_bar('Set',i);
		end
		spm_progress_bar('Clear',i);
		spm('FigName','Tor Batch Smooth: done',Finter,CmdLine);
		spm('Pointer');

	   catch

		warning('Error in smoothing...disk space? Permissions? Skipping this directory.')

	   end
	else
	end


	

end
