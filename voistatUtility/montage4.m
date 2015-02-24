function [T,redT] = montage4(anatomical,orient,range,varargin)
% function [T,redT] = montage(anatomical,orient,range,Tmap1,threshold,color,Tmap2,threshold,color...etc.)
% to PLOT ONLY:
% redT = montage4(T,'mask' OR 'real')
%
% Output: 3d arrays for anatomical, all t-masks specified
%
% Function montage4: overlays t-maps at thresholds you specify over anatomical images
% anatomical and Tmap1, etc. should be names of .img files in single quotes,
% LEAVING OFF the .img extension.
%
% RETURNS: 3-d volume matrix of all t-maps specified
% 			  t-maps (everything after 1st image) are masked by their respective thresholds
%			  1st output is vector of handles of slice axes
%
% orient specifies the display orientation.  'ax', 'sagg', or 'cor'
% 
% range specifies the range of slices you want to plot. e.g., [20 80], or 'all'
%
% example command:
% T = montage4('/data1/intext/sub1/Anatomy/t1','ax','spmT_0003',2,'r.')
% T = montage4([],'sagg',[30 90],[],2,'r.')
% T = montage4('s8t1','sagg','all','ravol_e3424_11_16_100_0008',1000,'rs')
% [T,redT] = montage4('s8t1','cor',[60 80],'ravol_e3424_11_16_100_0008',1000,'rs',[],1100,'ys')
%
% order of transformations to image:
%   flipy
%   rotate
%   scalex, scaley

% USER must specify THESE DEFAULTS in this SCRIPT
% ===========================================================================================
adjustfactor = 35; 		% you set this: more is less cut off of anatomical in case it's too big.
imgformat = 'Neuro';  	% Neuro, Rad, or Doug (for Doug's recon program as of Jan 2001)
%messedupSPMadjust = 0;  % if you (oops) wrote Tmaps in radio, but computed on Neuro data.

% output is in NEUROLOGICAL convention.
% ===========================================================================================

if nargin == 2, funct = 'plot';,elseif nargin > 2, funct = 'compute';,else error('incorrect argument spec.'),end

switch funct
	
case 'compute'
% make the T structure with everything you need to plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing maps.')
format compact
nimgs = 1 + ((nargin-3) / 3);
if nimgs < 1, nimgs = 1;, end
if ~(nimgs == round(nimgs)), error('enter 3 arguments per t-map.'), end
if strcmp(range,'all') | isempty(range), range = 'all';,end
	
% read anatomical
T.nimgs = nimgs;
T.orient = orient;
T.range = range;
T.name{1} = anatomical;
T.thresh(1) = 0;
T.color{1} = 'bone';			% the colormap of the anatomical
switch imgformat
case 'Neuro'
	[T.im{1},T.hdr{1}] = readim2(anatomical,orient);
	T.flipy = 'no';
case 'Rad'
	[T.im{1},T.hdr{1}] = readim2(anatomical,orient); 
	%switch orient
	%case 'ax'
	%	T.im{1} = flipdim(T.im{1},1)
	%case 'sagg'
	%	T.im{1} = flipdim(T.im{1},3)
	%case 'cor'
	%	T.im{1} = flipdim(T.im{1},2)
	%end
	T.flipy = 'no';
case 'Doug'
	[T.im{1},T.hdr{1}] = readim2(anatomical,orient,'flipy');
	T.flipy = 'yes';
end

% read other images
for i = 2:nimgs
	T.name{i} = varargin{3*i-5};
	switch imgformat
	case 'Neuro'
		[T.im{i},T.hdr{i}] = readim2(T.name{i},orient,'t');
	case 'Rad'
		[T.im{i},T.hdr{i}] = readim2(T.name{i},orient,'t');
		%switch orient
		%case 'ax'
		%	T.im{i} = flipdim(T.im{i},1)
		%case 'sagg'
		%	T.im{i} = flipdim(T.im{i},3)
		%case 'cor'
		%	T.im{i} = flipdim(T.im{i},2)
		%end
	case 'Doug'
		[T.im{i},T.hdr{i}] = readim2(T.name{i},orient,'t','flipy');
	end
	%if messedupSPMadjust
	%	switch orient
	%	case 'ax'
	%		T.im{i} = flipdim(T.im{i},1)
	%	case 'sagg'
	%		T.im{i} = flipdim(T.im{i},3)
	%	case 'cor'
	%		T.im{i} = flipdim(T.im{i},2)
	%	end
	%end
		
	T.thresh(i) = varargin{3*i-4};
	T.color{i} = varargin{3*i-3}; 
	T.realt{i} = T.im{i};                                   % original t-image copy
	T.im{i} = gettmask(T.im{i},T.thresh(i));			    % t-mask at specified threshold
    %T.scalex(i) = size(T.im{1},2) / size(T.im{i},2);
    %T.scaley(i) = size(T.im{1},1) / size(T.im{i},1);  
    % ======  get x,y on plot scalefactors ======== 
    switch orient
    case 'sagg'
        T.scalex(i) = T.hdr{i}.ysize / T.hdr{1}.ysize;
        T.scaley(i) = T.hdr{i}.zsize / T.hdr{1}.zsize;
    case 'cor'
        T.scalex(i) = T.hdr{i}.xsize / T.hdr{1}.xsize;
        T.scaley(i) = T.hdr{i}.zsize / T.hdr{1}.zsize;
    otherwise
        T.scalex(i) = T.hdr{i}.ysize / T.hdr{1}.ysize;  % axial: scalex = cols of matrix = y dim. in hdr
        T.scaley(i) = T.hdr{i}.xsize / T.hdr{1}.xsize;  %        because matlab images rows = y, cols = x.
    end % end switch
end

% ====== drop_n_slices : compare fields of view and cut down anatomical if necessary ========
% last argument is adjustment factor, in mm.  positive = cuts less off.
%T.hdr{1}.zsize = T.hdr{1}.zsize+.3;
%[T.im{1},drop_n_slices] = shrinkfov(T.im{1},orient,T.hdr{1},T.hdr{2},adjustfactor);

%redT = montage4(T,'real');
%redT = montage4(T,'mask');


case 'plot'
% plot the T structure as a montage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('plotting images.')
T = anatomical; funct = orient; 
if strcmp(funct,'real'),
    for i = 2:T.nimgs, T.mask{i} = T.im{i}; T.im{i} = T.realt{i};,end
end

for i = 1:T.nimgs
	zsize(i) = T.hdr{i}.zsize; 
	if strcmp(T.orient,'sagg'), zsize(i) = T.hdr{i}.xsize;			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif strcmp(T.orient,'cor'), zsize(i) = T.hdr{i}.ysize; 		% figure out voxel size in mm
	end			
	
	if i == 1														% specify range for 'all'
		if ischar(T.range)
			if T.range == 'all'
				T.range = [0 size(T.im{1},3) * zsize(i)];
			end
		end
	end

    T.originalslice{i} = 1:size(T.im{i},3);
	mm{i} = 1:size(T.im{i},3);
	mm{i} = mm{i} .* zsize(i) - zsize(i)/2;							% mm values of middle of ALL slices
	T.inrange{i} = mm{i} >= T.range(1) & mm{i} <= T.range(2);
    T.fullimdims{i} = size(T.im{i});                                % full image dims, for mask reconstruction
    T.fullmask{i} = T.im{i};
	T.im{i} = T.im{i}(:,:,T.inrange{i} == 1);						% cut all maps down to range
	mm{i} = mm{i}(:,T.inrange{i} == 1);                             % also cut mm and originalslice
	T.mm{i} = mm{i};
    T.originalslice{i} = T.originalslice{i}(:,T.inrange{i} == 1);   % to preserve record of which slice in vol.
end

% ======= plot the anatomical =========
[temp1,temp2,T.hand{1},T.whichslices,T.rows,T.cols,T.figh] = readim2(T.im{1},'p');		% plot, get whichslices and handles on reduced anatomical
colormap(T.color{1})
mm{1} = mm{1}(:,T.whichslices);										% restrict mm values to plotted slices
T.originalslice{1} = T.originalslice{1}(:,T.whichslices);	

% ======= set the data aspect ratio ========
for j = 1:size(T.hand{1},2)
            axes(T.hand{1}(j));
            dasp = get(gca,'DataAspectRatio'); 
            switch T.orient
            case 'sagg'
                daspect([dasp(1)*.78 dasp(2) 1])
            case 'cor'
                daspect([dasp(1)*.78 dasp(2) 1])
            otherwise
                %daspect([dasp(1)*1/.78 dasp(2) 1])
            end
end
          
% ====== label the images with mm values from origin =======
x = T.hdr{1}.origin(1) * T.hdr{1}.xsize;							% origin coordinates in mm
y = T.hdr{1}.origin(2) * T.hdr{1}.ysize;
z = T.hdr{1}.origin(3) * T.hdr{1}.zsize;
if strcmp(T.orient,'sagg'), omm = mm{1} - x; 			
elseif strcmp(T.orient,'cor'), omm = mm{1} - y;		
else omm = mm{1} - z;	
end			
for i = 1:size(T.hand{1},2)-1
	subplot(T.hand{1}(i))
	title([num2str(omm(i)) ' mm'])
end
drawnow																%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T.mm = mm;

% ====== figure out the closest slice for each t-map =======
nslices = size(mm{1},2);	
counter = 1:nslices;												% index slices plotted
for i = 2:T.nimgs													
	for j = 1:size(T.im{i},3)									
	  mmdif = abs(mm{1} - mm{i}(j));					
      	  temp = counter(mmdif == min(mmdif));		
      	  whichplot(j) = temp(1);                       			    % figure out plot index
      	  if whichplot(j) == 0, whichplot(j) = 1;,end					% for each tmap slice included
      	  if whichplot(j) > nslices, 
		  warning('picked slice greater than last slice.')
		  whichplot(j) = nslices;
	  end															%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	end
	T.whichplot{i} = whichplot;
    
    % ====== save a record of the original functional slice for each plot =======
    	T.osliceplot{i} = zeros(size(T.hand{1},1))                     % for each plot handle, the original functional slice
	for j = 1:size(whichplot,2)                                     % zero if no funct slice corresponding to this anatomical
        T.osliceplot{i}(whichplot(1,j)) = T.originalslice{i}(j);    % this indexes the original functional slices for each
	end                                                             % plot handle, to reconstruct z-values later
    

    % ===== plot the mask images or the actual t-maps ========
    if strcmp(funct,'real'),
        % plot the actual maps
        % ==================================================================================
        T = plotrealT(T,i);
    else
        % plot the on/off mask 
        % ==================================================================================
        plotmask(T.im{i},T.hand{1},T.color{i},[T.scalex(i) T.scaley(i)],T.whichplot{i});
        for j = 1:size(T.hand{1},2)
			set(T.hand{1}(j),'CameraViewAngle',max(T.rows,T.cols));
		end
    end
	drawnow
end

end % end switch
return


% ===================================== Sub-functions =====================================

function [anatarray,drop_n_slices] = shrinkfov(anatarray,ori,ahdr,thdr,adjustfactor)
% ====== compare fields of view and cut down anatomical if necessary ========
% this is needed if the images you're superimposing cover different regions of the brain/skull.
% it automatically cuts off the bottom of the anatomical, as needed...additional adjustment may be
% needed (adjustfactor) if the top of the skulls do not line up.
% last argument is adjustment factor, in mm.  positive = cuts less off.
fovdiff = ahdr.zsize * ahdr.zdim - thdr.zsize * thdr.zdim;
drop_n_slices = round((fovdiff - adjustfactor) / (2 * ahdr.zsize));         % *2 because assumes same center, shrinks top and bottom equally
if drop_n_slices > 0
    disp(['dropping ' num2str(drop_n_slices) ' from bottom of anatomical to match field of view with t-maps.'])
    switch ori
    case 'sagg'
        anatarray = anatarray(drop_n_slices+1:end-(drop_n_slices+7),:,:);
    case 'cor'
        anatarray = anatarray(drop_n_slices+1:end,:,:);
    otherwise
        anatarray = anatarray(:,:,drop_n_slices+1:end);
    end
end
return




function T = plotrealT(T,i)
        % ===== darken the anatomical axes =====
        axes(T.hand{1}(1));
        maxclim = get(gca,'CLim'); maxclim = maxclim(2);
        for j = 1:size(T.hand{1},2)
            axes(T.hand{1}(j)); set(gca,'CLim',[0 maxclim*1.7])
            %dasp = get(gca,'DataAspectRatio'); now done above
            %daspect([dasp(1)*.78 dasp(2) 1])
        end

        zoom = 6/max(T.rows,T.cols);
        colormap(jet)
        clim = [min(min(min(T.im{i}))) max(max(max(T.im{i})))]
        for j = 1:size(T.im{i},3)
            a = get(T.hand{1}(T.whichplot{i}(j)),'Position');
            T.hand{i}(j) = axes('Position',a);   
	        
            axis off; hold on; axis square
            
            slice = T.im{i}(:,:,j);
            for k = 1:size(slice,1)
               for m = 1:size(slice,2)
                  if (abs(slice(k,m)) > T.thresh(i)), newslice(k,m) = slice(k,m);,else newslice(k,m) = 0;,end
               end
            end
            
            T.imh{i}(j) = imagesc(newslice,[clim]);
            %set(T.imh{i}(j),'AlphaData',.5)
            set(gca,'YDir','reverse')
            alpha('color')
            switch T.orient
            case 'sagg'
                daspect([size(T.im{i},2)*.78 size(T.im{i},1) 1])
            case 'cor'
                daspect([size(T.im{i},2)*.78 size(T.im{i},1) 1])
            otherwise
                daspect([size(T.im{i},2)*.78 size(T.im{i},1) 1])
            end
            camzoom(zoom)
            camzoom(1.2)
            drawnow
        end   
        axes(T.hand{1}(end)); hold off;
        temp = imagesc(T.im{i}(1),[clim]); axis off;
        colorbar('horiz'); text(0,0,'t-scores');
        delete temp; cla   
return

