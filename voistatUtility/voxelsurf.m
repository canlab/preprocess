function vox = voxelsurf(redT,varargin)
% function vox = voxelsurf(mean_img_name,funct_basename,nimages,voxsize,avgover)
%
% function vox = voxelsurf(redT,vox)

if nargin > 2
	% mean functional volume name (no .img)
	% functional volume basename (no numbers or .img)
	% number of images

	redT2 = getfiles([redT]);
	if ~isempty(redT2),
		redT = redT2;
		redT = redT{1}(1:end-4);
	end
end

% not equipped to handle different voxel sizes on images, etc.
% define orthoh.
% if only T structure is entered, will load image for first use. 

if nargin > 2
    disp('setting up vox structure for first surfer use.')
%[T,redT] = montage4('t1','sagg',[80 100],'spmT_0002',3,'bs');
% im = readim2(redT.name{1});
    vox.im = readim2(redT);
    vox.basename = varargin{1};
    vox.nimages = varargin{2};
    vox.voxsize = varargin{3};
    vox.avgover = varargin{4};
	vox.HChoice = varargin{5};
	vox.TR = varargin{6};
    vox.orthoh = [];
    vox.voifigh = [];
    vox.pointh = [];
else
    % ======== don't set up, just get the timeseries. ==========
    vox = varargin{1};   
    figure(redT.figh)
    
    % get voxel coordinates
    
    vox.coord = getvoxels('single',[],redT);
    
    % image that slice on the original brain
    % figure;imagesc(vox.im(:,:,vox.coord(1,3)));colormap(gray);,hold on,plot(vox.coord(1,1),vox.coord(1,2),'bs')
    
    [vox.orthoh,vox.orthocoord] = ortho(vox.im,vox.coord,vox.orthoh,vox.voxsize);
    drawnow
    
    vox.ts = timeseries('voxel',vox.basename,vox.nimages,vox.coord);
    [vox.voifigh,vox.avgdata] = voistat('voifig',vox.voifigh,vox.ts,[],vox.avgover,vox.HChoice,vox.TR);
    drawnow
end   
return
 