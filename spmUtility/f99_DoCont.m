function status = f99_DoCont(act,contrast_num)

% ===================================================================
% function 'f99_DoCont'
% ===================================================================
% * this function computes contrast and statistical images
% * there are 3 options for argunment 1:
%     - 'all' : all contrasts will be computed
%     - 'selected' : those contrasts defined in arg.2 will be computed
%                    the exact names of the contrasts should be used!!
% PARAMETERS
% - option (see above)
% - contrasts number (in case of 'singleContrast'), see xCon.mat 
%   or list in contrast manager
% - contrast names (defined in f99_GetCont_...)
% ===================================================================

% SPM batch system: Contrast computation
% FORMAT status = spm_bch_DoCont
%_______________________________________________________________________
%
% This function implements the computations of the contrasts specified
% in the m-file.
%
% See also: spm_getSPM.m
%
%_______________________________________________________________________
% @(#)spm_bch_DoCont.m	2.6 Jean-Baptiste Poline & Stephanie Rouquette 99/10/29
% Derived from spm_getSPM.m

%- initialise status
%-----------------------------------------------------------------------
swd = pwd;
Vbeta = [];

%- Get xCon if exist, if doesnt warns and returns 
%-----------------------------------------------------------------------

if exist(fullfile('.','xCon.mat'),'file'), 
	load('xCon.mat'); 
	lxCon = length(xCon);

	%-Canonicalise SPM99b format xCon (which saved mmapped handles) **
	%---------------------------------------------------------------
	for i=1:length(xCon)
	   if isstruct(xCon(i).Vcon), xCon(i).Vcon=xCon(i).Vcon.fname; end
	   if isstruct(xCon(i).Vspm), xCon(i).Vspm=xCon(i).Vspm.fname; end
	end

else, 
	str = ['cannot open ' fullfile('.','xCon.mat') ...
               ' file in spm_bch_DoCont ' swd];
	warning(str);
	return;
end


if exist(fullfile('.','SPM.mat'),'file'), 
	try 
		load(fullfile('.','SPM.mat'),'xX');	
		load(fullfile('.','SPM.mat'),'Vbeta');	
		load(fullfile('.','SPM.mat'),'XYZ');	
		load(fullfile('.','SPM.mat'),'VResMS');	
	catch 
		str = ['cannot open ' fullfile('.','SPM.mat') ...
				' file in spm_bch_GetCont ' swd];
		warning(str);
		status.str = str;
		status.err = 2;
		return;
	end
else 
	str = ['cannot find ' fullfile('.','SPM.mat') ...
				' file in spm_bch_GetCont ' swd];
	warning(str);
	return;
end

% TOR: SPM2 version
if isempty(Vbeta)
    load SPM
    N = fieldnames(SPM);
    for i = 1:length(N), eval([N{i} ' = SPM.' N{i} ';']);,end
end


if isempty(Vbeta), 
    warning('no beta, no contrasts'), return, 
else
    %- remap the files 
    Vbeta  = spm_vol([repmat([swd,filesep],length(Vbeta),1),char(Vbeta)]);
    VResMS = spm_vol(fullfile(swd,VResMS));
end

M     	= Vbeta(1).mat;
dim	= Vbeta(1).dim(1:3);

%- loop over contrasts...
%-----------------------------------------------------------------------
% select contrasts to compute
switch act
  case 'all'
    conmatrix = 1:length(xCon);
  case 'selected'
    min_Ic=1;
    max_Ic=size(contrast_num,1);
    conmatrix = [];
    for t = 1:size(contrast_num,1)
        l = 1;
        ff = 1;
        while  l<=length(xCon) & ff>0
            if strcmp(deblank(xCon(l).name),deblank(strcat(contrast_num(t,:))))
              conmatrix = [conmatrix l];
              ff = ff*-1;
            end
            l = l+1;
        end
        if ff>0
            fprintf(['\n contrast name ' mat2str(xCon(t).name) 'was not found \n']);
        end
    end
end
conmatrix
for i = conmatrix

    if ~isfield(xCon(i),'eidf') | isempty(xCon(i).eidf)
	[trMV,trMVMV] = spm_SpUtil('trMV',...
				spm_FcUtil('X1o',xCon(i),xX.xKXs),xX.V);
        xCon(i).eidf  = trMV^2/trMVMV;
    else
        trMV = []; trMVMV = [];
    end

    %-Write contrast/ESS images?
    %-------------------------------------------------------------------
    if ~isfield(xCon(i),'Vcon') | isempty(xCon(i).Vcon) | ...
        ~exist(fullfile('.',xCon(i).Vcon),'file')
        	
	switch(xCon(i).STAT)

	case 'T'       %-Implement contrast as sum of scaled beta images
	%---------------------------------------------------------------
            fprintf('\t%-32s: %-10s%20s',sprintf('contrast image %2d',i),...
				      '(spm_add)','...initialising') %-#

	    Q     = find(abs(xCon(i).c) > 0);
	    V     = Vbeta(Q);
	    for j = 1:length(Q)
		V(j).pinfo(1,:) = V(j).pinfo(1,:)*xCon(i).c(Q(j));
	    end
	    
	    %-Prepare handle for contrast image
	    %-----------------------------------------------------------
	    xCon(i).Vcon = struct(...
	        'fname',  sprintf('con_%04d.img',i),...
	        'dim',    [dim,16],...
	        'mat',    M,...
	        'pinfo',  [1,0,0]',...
	        'descrip',sprintf('SPM contrast - %d: %s',i,xCon(i).name));

	    %-Write image
	    %-----------------------------------------------------------
	    fprintf('%s%20s',sprintf('\b')*ones(1,20),'...computing')%-#
	    xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);
	    xCon(i).Vcon.pinfo(1,1) = spm_add(V,xCon(i).Vcon);
	    xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);

	    fprintf('%s%30s\n',sprintf('\b')*ones(1,30),...
			sprintf('...written %s',xCon(i).Vcon.fname)) %-#

	case 'F'  %-Implement ESS as sum of squared weighted beta images
	%---------------------------------------------------------------
	    fprintf('\t%-32s: %30s',sprintf('ESS image %2d',i),...
						     '...computing') %-#

            %-Residual (in parameter space) forming mtx
	    %-----------------------------------------------------------
	    h       = spm_FcUtil('Hsqr',xCon(i),xX.xKXs);

	    %-Prepare handle for ESS image
	    %-----------------------------------------------------------
	    xCon(i).Vcon = struct(...
		'fname',  sprintf('ess_%04d.img',i),...
		'dim',    [dim,16],...
		'mat',    M,...
		'pinfo',  [1,0,0]',...
		'descrip',sprintf('SPM ESS - contrast %d: %s',i,xCon(i).name));

	    %-Write image
	    %-----------------------------------------------------------
	    fprintf('%s',sprintf('\b')*ones(1,30))                   %-#
	    xCon(i).Vcon  = spm_create_image(xCon(i).Vcon);
	    xCon(i).Vcon  = spm_resss(Vbeta,xCon(i).Vcon,h);
	    xCon(i).Vcon  = spm_create_image(xCon(i).Vcon);

	otherwise
	    error(['unknown STAT "',xCon(i).STAT,'"'])

	end % switch(xCon(i).STAT)

    end % if ~isfield(xCon(i)

    %-Write statistic image(s)
    %-------------------------------------------------------------------
    if ~isfield(xCon(i),'Vspm') | isempty(xCon(i).Vspm) | ...
        ~exist(fullfile('.',xCon(i).Vspm),'file')

	fprintf('\t%-32s: %30s',sprintf('spm{%c} image %2d',xCon(i).STAT,i),...
                                                    '...computing')  %-#

	switch(xCon(i).STAT)
	case 'T'                                  %-Compute SPM{t} image
	%---------------------------------------------------------------
	Z   = spm_sample_vol(xCon(i).Vcon, XYZ(1,:),XYZ(2,:),XYZ(3,:),0)./...
		(sqrt(spm_sample_vol(VResMS,  XYZ(1,:),XYZ(2,:),XYZ(3,:),0)*...
					(xCon(i).c'*xX.Bcov*xCon(i).c) ));
	str = sprintf('[%.2g]',xX.erdf);
	
	case 'F'                                  %-Compute SPM{F} image
	%---------------------------------------------------------------
	if isempty(trMV)
	   	trMV = spm_SpUtil('trMV',spm_FcUtil('X1o',xCon(i),xX.xKXs),xX.V);
	end
	Z =(spm_sample_vol(xCon(i).Vcon,XYZ(1,:),XYZ(2,:),XYZ(3,:),0)/trMV)./...
           (spm_sample_vol(VResMS, XYZ(1,:),XYZ(2,:),XYZ(3,:),0));
	
	str = sprintf('[%.2g,%.2g]',xCon(i).eidf,xX.erdf);
	
	otherwise
	%---------------------------------------------------------------
            error(['unknown STAT "',xCon(i).STAT,'"'])

        end %- switch(xCon(i).STAT)

        %-Write full statistic image
        %---------------------------------------------------------------
        fprintf('%s%30s',sprintf('\b')*ones(1,30),'...writing')      %-#

        xCon(i).Vspm = struct(...
	    'fname',  sprintf('spm%c_%04d.img',xCon(i).STAT,i),...
	    'dim',    [dim,16],...
	    'mat',    M,...
	    'pinfo',  [1,0,0]',...
	    'descrip',sprintf('SPM{%c_%s} - contrast %d: %s',...
	                       xCon(i).STAT,str,i,xCon(i).name));

	tmp = zeros(dim);
	tmp(cumprod([1,dim(1:2)])*XYZ - sum(cumprod(dim(1:2)))) = Z;

	xCon(i).Vspm  = spm_write_vol(xCon(i).Vspm,tmp);

	fprintf('%s%30s\n',sprintf('\b')*ones(1,30),...
          sprintf('...written %s',xCon(i).Vspm.fname)) %-#

    end %- if ~isfield(xCon(i),'Vspm')

end %- for i = 1:length(xCon)

%- save xCon 
%-----------------------------------------------------------------------
try
	save('xCon.mat','xCon')
catch
	str = ['Can''t write xCon.mat to the results directory: ' swd];
	warning(str);
end
