function varargout = slice_overlay(varargin);
% Function to display + manage slice display 
% Slice display works on a global structure SO
% with fields 
%  - img - array of images to display
%        - img structs contain fields
%             vol - vol struct info (see spm_vol)
%             cmap - colormap for this image
%             nancol - color for NaN. If scalar, this is an index into
%                    the image cmap.  If 1x3 vector, it's a colour
%             prop - proportion of intensity for this cmap/img
%                    if not in inclusive range 0-1, gives split cmap
%                    effect where values of this cmap override previous
%                    image cmap values
%             func - function to apply to image before scaling to cmap
%                    (and therefore before min/max thresholding. E.g. a func of
%                    'i1(i1==0)=NaN' would convert zeros to NaNs
%             range - 2x1 vector of values for image to distribute colormap across
%                    the first row of the colormap applies to the first
%                    value in 'range', and the last value to the second
%                    value in 'range'
%             outofrange - behavior for image values to the left and
%                    right of image limits in 'range'.  Left means
%                    colormap values < 1, i.e for image values <
%                    range(1), if (range(1)<range(2)), and image values >
%                    range(1) where (range(1)>range(2)). If missing,
%                    display min (for Left) and max (for Right) value from colormap. 
%                    Otherwise should be a 2 element cell array, where
%                    the first element is the colour value for image values
%                    left of 'range', and the second is for image values
%                    right of 'range'.  Scalar values for
%                    colour index the colormap, 3x1 vectors are colour
%                    values.  Empty arrays default to black for
%                    transparent colour (where SO.prop <= 1), and empty
%                    for split colour  
%            
% - transform - 4x4 transformation to apply to image slice position,
%             relative to mm given by slicedef, before display
% - slicedef - 2x3 array specifying dimensions for slice images in mm
%             where rows are x,and y of slice image, and cols are neg max dim,
%             slice separation and pos max dim
% - slices   - vector of slice positions in mm in z (of transformed image)
% - figure    - figure handle for slice display figure
% - xsliceno  - no of slices to display across figure (defaults to an optimum)
% - hold      - resampling order for images (see spm_sample_vol)
% - refreshf  - flag - if set or empty, refresh axis info for figure
%             else assume this is OK
% - cbar      - if empty, missing, no colourbar.  If an array of integers, then indexes
%             img array, and makes colourbar for each cmap for that img
% - labels - struct can be missing (no labels) or contain fields 
%                  colour - colour for label text 
%                  size - font size in units normalized to slice axes 
%                  format - if = cell array of strings =
%                  labels for each slice in Z.  If is string, specifies
%                  sprintf format string for labelling in distance of the
%                  origin (Xmm=0, Ymm=0) of each slice from plane containing
%                  the AC, in mm, in the space of the transformed image
%  
%  
% V 0.2 3/5/00  
% Way way left of beta - take great care.  Please report problems to  
% Matthew Brett matthew@mrc-cbu.cam.ac.uk

global SO 


if nargin < 1
  action = 'display';
else
  action = lower(varargin{1});
end

switch action
 case 'check'
  checkso;
 case 'display'

% check and/or initialise SO struct
% ------------------------------------------------------------------------------------
checkso;

% get coordinates for plane
% ------------------------------------------------------------------------------------
X=1;Y=2;Z=3;
dims = SO.slicedef;
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);
zmm = SO.slices;
[y x] = meshgrid(ymm,xmm');
vdims = [length(xmm),length(ymm),length(zmm)];

% no of slices, and panels (an extra for colorbars)
% ------------------------------------------------------------------------------------
nslices = vdims(Z);
npanels = nslices;
cbars = 0;
if is_there(SO,'cbar')
  cbars = length(SO.cbar);
%  cbars = size(SO.cbar,1)
  scbars = SO.scbar;
  npanels = npanels+cbars;
end

% get figure data
% ------------------------------------------------------------------------------------
% if written to, the axes may be specified already
figno = figure(SO.figure);
axisd = findobj(gcf, 'Type','axes','Tag', 'slice panel');
if isempty(axisd) | length(axisd) < npanels
  SO.refreshf = 1;
end
t = SO.orient;

% extract descriptions
% ------------------------------------------------------------------------------------
Descrip1 = SO.Descrip1;
Descrip2 = SO.Descrip2;
Descrip3 = SO.Descrip3;

% (re)initialize axes and stuff
% ------------------------------------------------------------------------------------
if SO.refreshf
  % clear figure and axis store
  clf
  axisd = [];

  % figure dimensions, in pixels
  figubu = get(figno, 'Units');
  set(figno, 'Units', 'Pixels');
  figdims = get(figno, 'Position');
  figsz = figdims(3:4);
  set(figno, 'Units', figubu);
  set(figno, 'Color',[1 1 1]);
  

  % prevent print inversion problems
  set(figno,'InvertHardCopy','off');
  
  % by default, make most parsimonious fit to figure
  yxratio = length(ymm)*dims(Y,2)/(length(xmm)*dims(X,2));
  if ~is_there(SO, 'xslices')
    % iteration needed to optimize, surprisingly.  Thanks to Ian NS
    xaxlen=figsz(1):-1:1;
    yaxlen=yxratio*xaxlen;
    xslices = floor(figsz(1)./xaxlen);
%    yslices  = floor(figsz(2)./yaxlen);
    yslices  = floor((figsz(2)-70)./yaxlen); % ruimte boven vrij laten
    estnpanels = xslices.*yslices;
    tmp = find(estnpanels >= npanels);
    if isempty(tmp)
      error('Whoops, cannot fit panels onto figure');
    end
    b = tmp(1); % best fitting scaling
    xslices = xslices(b);
    xaxlen = xaxlen(b);
    xaxlen = xaxlen-2;
  else
    % if xslices is specified, assume X is flush with X figure dimensions
    xslices = SO.xsliceno;
    xaxlen = figsz(X)/xslices;
    xaxlen = xaxlen-2;
  end
  
  % Axis dimensions are in pixels.  This prevents aspect ratio rescaling
  yslices = ceil(npanels/xslices);
  yaxlen = xaxlen*yxratio;
%  cstart = figsz(Y)-yaxlen*yslices;
  cstart = figsz(Y)-(yaxlen*yslices)-70;
  
  % make axes for panels
  r=0;c=1;
  for i = 1:xslices*yslices
    axpos = [(r*xaxlen)+(r*2) (yslices-c)*yaxlen+cstart-(c*2) xaxlen yaxlen];
    axisd(i) = axes(...
	'Parent',figno,...
	'XTick',[],...
	'XTickLabel',[],...
	'YTick',[],...
	'YTickLabel',[],...
	'Box','on',...
	'XLim',[1 vdims(X)],...
	'YLim',[1 vdims(Y)],...
	'Units', 'pixels',...
	'Position',axpos,...
	'Tag','slice panel'...
	);
    r = r+1;
    if r >= xslices
      r = 0;
      c = c+1;
    end
  end
end

% sort out labels
% ------------------------------------------------------------------------------------
if is_there(SO,'labels')
  labels = SO.labels;
  labels.format;
  if iscell(labels.format)
    if length(labels.format)~=vdims(Z)
      error(...
	  sprintf('Oh dear, expecting %d labels, but found %d',...
		  vdims(Z), length(labels.contents)));
    end
  else
    % format string for mm from AC labelling
    fstr1 = labels.format;
    labels.xlabel_sag = '%s';
    fstr2 = labels.xlabel_sag;
    labels.format      = cell(vdims(Z),1);
    labels.xlabel_sag  = cell(vdims(Z),1);
    acpt = SO.transform * [0 0 0 1]';
    for i = 1:vdims(Z)
      if t==1
        labels.format(i) = {sprintf(fstr1,zmm(i)-acpt(Z))};
      elseif t==2
        labels.format(i) = {sprintf(fstr1,-1*(zmm(i)-acpt(Z)))};
%        labels.format(i) = {sprintf(fstr1,(zmm(i)-acpt(Z)))};
      elseif t==3
        labels.format(i) = {sprintf(fstr1,zmm(i)-acpt(Z))};
        if zmm(i)-acpt(Z)>0
          labels.xlabel_sag(i) = {sprintf(fstr2,['R'])};
        elseif zmm(i)-acpt(Z)<0
          labels.xlabel_sag(i) = {sprintf(fstr2,['L'])};
        elseif zmm(i)-acpt(Z)==0
          labels.xlabel_sag(i) = {sprintf(fstr2,[])};
        end
      end
    end
  end
end


% cycle through slices displaying images
% ------------------------------------------------------------------------------------
nimgs = length(SO.img);
nvox = prod(vdims(1:2));
pandims = [vdims([2 1]) 3]; % NB XY transpose for display

zimg = zeros(pandims);
for i = 1:nslices
  ixyzmm = [x(:)';y(:)';ones(1,nvox)*zmm(i);ones(1,nvox)];
  img = zimg;
  for j = 1:nimgs
    thisimg = SO.img(j);
    % to voxel space of image
    vixyz = inv(SO.transform*thisimg.vol.mat)*ixyzmm;
    % raw data 
    i1 = spm_sample_vol(thisimg.vol,vixyz(X,:),vixyz(Y,:),vixyz(Z,:), ...
			 SO.hold);
    if is_there(thisimg, 'func')
      eval(thisimg.func);
    end
    % transpose to reverse X and Y for figure
    i1 = reshape(i1, vdims(1:2))';
    % rescale to colormap
    if ~is_there(thisimg,'outofrange')
      thisimg.outofrange = {[0 0 0],[0 0 0]};
    end
    [csdata cmap lvals rvals nans] = scaletocmap(...
	i1,...
	thisimg.range(1),...
	thisimg.range(2),...
	thisimg.cmap,...
	thisimg.outofrange{1},...
	thisimg.outofrange{2},...
	thisimg.nancol...	
    );
    % take indices from colormap to make true colour image
    iimg = reshape(cmap(csdata(:),:),pandims);
    if thisimg.prop < 1 % truecolor overlay
      img = img + iimg*thisimg.prop;
    else % split colormap effect
	 % deal with nans and values above and below range
	 msk = nans;
	 if isempty(thisimg.outofrange{1})
	   msk = msk | lvals;
	 end
	 if isempty(thisimg.outofrange{2})
	   msk = msk | rvals;
	 end
	 tmp = repmat(~msk,[1 1 3]);
	 img(tmp) = iimg(tmp);
    end
  end
  % threshold out of range values
  img(img>1) = 1;
  image('Parent', axisd(i),...
	'CData',img);
  if is_there(SO,'labels')
    text('Parent',axisd(i),...
         'Position',[3 2],...
	 'Color', labels.colour,...
	 'FontUnits', 'normalized',...
	 'VerticalAlignment','bottom',...
	 'FontSize',labels.size,...
	 'String', labels.format{i});
    if t==1
      text('Parent',axisd(i),'Position',[5 ...
            100],'String','R','FontSize',10,'Color', labels.colour)
      text('Parent',axisd(i),'Position',[80 ...
            100],'String','L','FontSize',10,'Color', labels.colour)
    elseif t==2
      text('Parent',axisd(i),'Position',[4 ...
            130],'String','R','FontSize',10,'Color', labels.colour)
      text('Parent',axisd(i),'Position',[85 ...
            130],'String','L','FontSize',10,'Color', labels.colour)
    elseif t==3
      text('Parent',axisd(i),'Position',[4 ...
            60],'String',labels.xlabel_sag{i},'FontSize',10,'Color', labels.colour)
    end
  end
end
for i = (nslices+1):xslices*yslices
   set(axisd(i),'Color',[0 0 0]);
end

% add colorbar(s) 
% ------------------------------------------------------------------------------------
for i = 1:cbars
  axno = axisd(end-cbars+i);
  cbari = SO.img(SO.cbar(i));
  cml = size(cbari.cmap,1);
  p = get(axno, 'Position'); % position of last axis
  cw = p(3)*0.2;
  ch = p(4)*0.75;
  pc = p(3:4)/2;
  [axlims idxs] = sort(cbari.range);
  a=axes(...
      'Parent',figno,...
      'XTick',[],...
      'XTickLabel',[],...
      'Units', 'pixels',...
      'YLim', axlims,...   
      'FontUnits', 'normalized',...
      'FontSize', 0.075,...
      'YColor',[1 1 1],...
      'Tag', 'cbar',...
      'Box', 'off',...
      'Position',[p(1)+pc(1)-cw/2,p(2)+pc(2)-ch/2+10,cw,ch]...
      );
  [p(1)+pc(1)-cw/2,p(2)+pc(2)-ch/2+10,cw,ch]
  ih = image('Parent', a,...
	'YData', axlims(idxs),...
	'CData', reshape(cbari.cmap,[cml,1,3]));
  if t==1, pos=[47 10];
     elseif t==2, pos=[47 6]; 
     elseif t==3, pos=[55 6];
  end
  text('Parent',axno,...
       'Position',pos,...
       'HorizontalAlignment','center',...
       'String',scbars(i,:),'FontSize',14,'Color', [1 1 1])
end

% Add comments on top
% ------------------------------------------------------------------------------------
figdim   = get(figno,'Position');
figdim   = figdim(3:4);
pos_x1   = figdim(1)*0.012;
pos_x2   = figdim(1)*0.4796;
pos_y    = figdim(2)*0.941;
len_x    = figdim(1)*0.3597;
len_y    = figdim(2)*0.052;
leg = axes(...
      'Parent',figno,...
      'XTick',[],...
      'YTick',[],...
      'XTickLabel',[],...
      'Box','on',...
      'Units', 'pixels',...
      'Position',[pos_x1 pos_y 300 60],...
      'Color',[.8 .8 .8]);

text('Parent',leg,...
     'Position',[0.03 0.8],...
     'Color',[0 0 0],...
     'String',Descrip1,...
     'Fontsize',10,...
     'FontWeight','bold')
 
text('Parent',leg,...
     'Position',[0.03 0.22],...
     'Color',[0 0 0],...
     'String',Descrip2,...
     'Fontsize',10)
 
leg = axes(...
      'Parent',figno,...
      'XTick',[],...
      'YTick',[],...
      'XTickLabel',[],...
      'Box','on',...
      'Units', 'pixels',...
      'Position',[pos_x2 pos_y 400 60],...
      'Color',[.8 .8 .8]);

text('Parent',leg,...
     'Position',[0.03 0.5],...
     'Color',[0 0 0],...
     'String',Descrip3,...
     'Fontsize',10)
 
% end switch action
end

return


% '======================================================================='
function checkso
% '======================================================================='
% checks and fills SO structure
% (sometime)
global SO

return


% '======================================================================='
function tf = is_there(a, fname)
% '======================================================================='
% returns true if field fname is present in struct a, and not empty
tf = isfield(a, fname);
if tf
  tf = ~isempty(getfield(a, fname));
end
return


% '======================================================================='
function [img,cmap,lvals,rvals,nans]=scaletocmap(inpimg,mn,mx,cmap,l,r,n)
% '======================================================================='
cml = size(cmap,1);
scf=(cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
rvals = img>cml;
if isempty(r)
  r = [0 0 0];
end
if prod(size(r))==3
  cmap = [cmap; r(:)'];
  r = size(cmap, 1);
end
img(rvals) = r;
lvals = img<1;
if isempty(l)
  l = [0 0 0];
end
if prod(size(l))==3
  cmap = [cmap; l(:)'];
  l = size(cmap, 1);
end
img(lvals) = l;
nans = isnan(img);
if isempty(n)
  n = [0 0 0];
end
if prod(size(n))==3
  cmap = [cmap; n(:)'];
  n = size(cmap, 1);
end
img(nans) = n;
return

