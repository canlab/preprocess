function dstats = descriptives(in,varargin)
% function dstats = descriptives(ts/mask,basename [opt],nimages[opt],covariate [opt],threshold [opt],aformat[opt])
% OR, if ts is a timeseries, dstats = descriptives(ts,covariate [opt])
% **********************************
% in is a timeseries of a voxel OR mask OR [] (default 64 x 64 x 28 whole vol.)
% 		pre-loaded into Matlab
%
% covariate: regresses out cov. matrix and temporal derivative before computing autocorrelation.
%			enter just the regressor, with no derivative or intercept.  these added automatically.
%
% example: 
%			s8dstats = descriptives([],'ravol_e3424_11_16_100_',200,X(:,2));
%			this gets the whole-volume stats for 200 images with a covariate.
%
% nimages is either a number or a row vector of specific image numbers.
% other options: threshold, aformat 'long' or 'short' to give just simple stats; default is long.
% 
%
% Tor Wager, March 2001.  Last edit 3/08/01 to add the last feature above.

% set up arguments
if isempty(in) | size(size(in),2) > 2
    funct = 'volume';
elseif size(in,2) == 1
    funct = 'voxel';
elseif size(in,2) > 1 & size(size(in),2)
    funct = 'multi';
else error('enter ts col. vector(s) or volume mask or [] for 1st arg.')
end
if nargin > 3
	basename = varargin{1};,nimages = varargin{2};cov = varargin{3};
elseif nargin > 2
    basename = varargin{1};,nimages = varargin{2};
elseif nargin > 1 & nargin < 3
	cov = varargin{1};
end

if exist('nimages'),
   if size(nimages,2) == 1, nimages = 1:nimages;,end
end

% 600 threshold is appropriate for intext.  tested 3/03/01.
ARsize = 16;                                                        % default order for Yule-Walker autoregressive.
threshold = 600;                                                    % default masking threshold for whole volume
if nargin > 4, if ~isempty(varargin{4}),threshold = varargin{4};,end,end
aformat = 'long';																 % format of program: long = all calculations,
if nargin > 5, if ~isempty(varargin{5}),aformat = varargin{5};,end,end
																					% short = omit most.
switch funct

case 'voxel'
   %disp('input is timeseries vector for one voxel.')
   if size(in,2) < size(in,1), in = in';,end						% make a col vector
   dstats.mean = mean(in);
   dstats.std = std(in);
   dstats.min = min(in);
   dstats.max = max(in);
   if strcmp(aformat,'long')
   		temp(:,1) = abs(fft(in))';									    % compute power
   		dstats.fft(:,1) = temp(1:ceil(size(temp,1)/2),1);	            % take only 1st half of power
   		%dstats.ar = aryule(in,ARsize);
   		if exist('cov') == 1
			if ~isempty(cov)
	   			dstats.xc = getv('get',in',cov);
       	 else dstats.xc = getv('get',in');  
       	 end
   		else
   			dstats.xc = getv('get',in');
      end
   end
      
case 'volume'
   disp('input is basename and nimages for a volume.')
   disp(['output format is ' aformat]);
    sumfft = zeros(size(nimages,2)/2,1);
    %sumar = zeros(1,ARsize + 1);
    sumxc = zeros(1,floor(size(nimages,2)/4)+1);
    countfft = 0;
    if isempty(in),                                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        disp(['loading reference image and thresholding at ' num2str(threshold)])     % if no mask is entered, use whole volume...   %
       		                                      %                                              %
       [ref,hdr] = readim2([basename '0001']);                           %        % load 1st image as reference         %
       in = zeros(hdr.xdim,hdr.ydim,hdr.zdim);  
       disp(['mask dimensions are ' num2str([hdr.xdim hdr.ydim hdr.zdim])])
        in(ref > threshold) = 1;    clear ref;  drawnow;                  %        % threshold at 600                    %
    end                                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for z = 1:size(in,3)
       disp([' loading slice ' num2str(z) '...']); drawnow;
        slice = timeseries('slice',basename,nimages,z); % load slice
        slice = permute(slice,[4 1 2 3]);
        disp(['     computing descriptives...']); drawnow;
        for x = 1:size(in,2)
            fprintf('.')
            for y = 1:size(in,1)
                if in(y,x,z) > 0
                   ts = slice(:,y,x);
						if nargin > 3
							ds = descriptives(ts,[],[],cov,[],aformat);
						else
                     ds = descriptives(ts,[],[],[],[],aformat);
                   end
                 	dstats.mean(y,x,z) = ds.mean;
   			        	dstats.std(y,x,z) = ds.std;
   			        	dstats.min(y,x,z) = ds.min;
                 	dstats.max(y,x,z) = ds.max;
                 	if strcmp(aformat,'long')
                 		sumfft = ds.fft + sumfft;
                    	countfft = countfft + 1;
                    	%sumar = sumar + ds.ar;
                    	sumxc = sumxc + ds.xc;
                 	end   
                 end
            end
        end
        clear slice;
     end 
   
     disp('computing avg power spectrum...')
     if ~(exist('dstats')),dstats = [];,end
     if ~(isfield(dstats,'std')),dstats.std = 0;,end
     
     if strcmp(aformat,'long')
     	dstats.avgfft = sumfft / countfft;     
   		%dstats.avgAR = sumar / countfft; 
   		dstats.avgxc = sumxc / countfft; 
         %dstats.V = sparse(getv('make',dstats.avgxc,size(ts,1)));
     end
   	  dstats.noisevar = mean(mean(mean(dstats.std))) .^ 2;
   
case 'multi'
    disp(['computing stats for multiple voxels (' num2str(size(in,2)) ') entered...'])
	nimages = size(in,1);
    sumfft = zeros(nimages/2,1);
    %sumar = zeros(1,ARsize + 1);
    sumxc = zeros(1,floor(nimages/4)+1);
    countfft = 0;
    for i = 1:size(in,2)
       if nargin > 4
			ds = descriptives(in(:,i),cov);
	    else
            ds = descriptives(in(:,i));
		 end
        dstats.mean(i) = ds.mean;
        dstats.std(i) = ds.std;
   	     dstats.min(i) = ds.min;
        dstats.max(i) = ds.max;
     	 if strcmp(aformat,'long')
        	sumfft = ds.fft + sumfft;
        	countfft = countfft + 1;
        	%sumar = sumar + ds.ar;
           sumxc = sumxc + ds.xc;
        end  
        fprintf('.')
     end 
   disp('computing avg power spectrum...')
   if strcmp(aformat,'long')
   		dstats.avgfft = sumfft / countfft;     
   		%dstats.avgAR = sumar / countfft; 
   		dstats.avgxc = sumxc / countfft; 
   		%dstats.V = sparse(getv('make',dstats.avgxc,nimages));
   end      
   dstats.noisevar = mean(mean(mean(dstats.std))) .^ 2;
   
end		% end switch
return


