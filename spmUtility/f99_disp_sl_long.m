function f99_display_slices(dir_struct,struct,dir_activ,act,struct_intens,range_factor,colmap,tmap_range,s_t,col_intens,col_bar,t,lx,hx,dx,name,comment1,comment2,comment3,write_file)

% '==================================================================================='
%                   'Overlay of activations on a structural image'
% '==================================================================================='

% ====================================================================================
%                             --- PARAMETERS ---
% ====================================================================================

% 'where to find the images'
% - dir_struct    : directory where to get the structural image
% - struct        : name and path of T1-anatomical image (relative to fmriDIR)
% - dir_activ     : - directory where to get the activation images
%                   - directory where to write overlay images!!!
% - act           : name and path of activation-image (spmT!!!) 
%                   (relative to fmriDIR)
% 'adjusting the intensity of your structural image'
% - struct_intens : scale factor for intensity of structural image
% - range_factor  : sacle factor to reduce the upper limit of the
%                   intensity range of the structural image
%                   e.g. in case the scull is hyperintense
% 'properties of the activation image'
% - colmap        : matrix with string of names of colormaps for activations
% - tmap_range    : minimum and maximal t-value to rescale the t-values
% - s_t           : Transparent or Split overlay [1 2]
% - col_intes     : intensity of colormap (for transparent overlay)
% 'what panels to display?'
% - col_bar       : matrix of string names for different colorbars 
%                                               (empty if no)
% - t             : oreintation of the slices (1=axial; 2=coronal; 3=sagittal)
% - lx            : position of the lowest slice
% - hx            : position of the highest slice
% - dx            : slice separation in mm
% 'adding comments'
% - name          : name of the overlay to be displayed on top 
%                   also the name of the saved file!!!
% - comment1      : any comment of maximal 53 characters (e.g. subtraction name)
% - comment2      :   idem
% - comment3      :   idem
% 'saving the overlay image?'
% - write_file    : write overlay images to disk? (1=yes / 0=no)
% ------------------------------------------------------------------------------------

% request a whole load of parameters for slice_overlay routine
% SO structure contains all the parameters for display
% See slice_overlay.m for detailed comments

clear global SO
global SO fmriDIR


% change to working directory
% ---------------------------
str = ['cd ' fmriDIR];
eval(str)
CWD = pwd;
spm_input('!SetNextPos', 1);


% load images
% ------------------------------------------------------------------------------------
% structural image
str = ['cd ' dir_struct];
eval(str)
img1   = [pwd filesep struct];
% activation image
str = ['cd ' CWD];
eval(str)
%str = ['cd ' dir_activ];
%eval(str)
for k = 1:size(act,1)
  str = ['cd ' dir_activ(k,:)];
  eval(str)
  img2s(k,:) = [pwd filesep act(k,:)];
  str = ['cd ' CWD];
  eval(str)
end
noaimgs = size(img2s, 1);

str = ['cd ' CWD];
eval(str)


% params for structural image
% ------------------------------------------------------------------------------------
SO.img(1).vol = spm_vol(img1);
img = spm_read_vols(SO.img(1).vol); % for min and max
SO.img(1).range = [min(img(:)) max(img(:))/range_factor];
SO.img(1).cmap = gray; % structural colormap
%if noaimgs
%  SO.img(1).prop = spm_input('Intensity for structural', '+1', 'r', 0.6, 1, [0 ...
%		    2]);
%else
%  SO.img(1).prop = 1;
%end
SO.img(1).prop = struct_intens;
SO.img(1).nancol = [0 0 0]; % NaNs are black


% load stuff for activation images
% ------------------------------------------------------------------------------------
remcol = 1-SO.img(1).prop;
SO.cbar = [];
for i = 1:noaimgs;
  imgno = i+1;
  % colormap
%  actc = [];
%  while isempty(actc)
%    actc = getcmap(spm_input(['Colormap - image ' num2str(i)],...
%				  '+1','s', 'actc'));
%  end
  map = getcmap(deblank(colmap(i,:)));
  SO.img(imgno).cmap = map; % colormap for activation
  
  % params for functional image
  SO.img(imgno).vol = spm_vol(deblank(img2s(i,:)));
  img = spm_read_vols(SO.img(imgno).vol);
  img = img(:);
  img_range = [min(img) max(img)];
%  SO.img(imgno).range = spm_input('Img val range for colormap',...
%				  '+1', 'e', [min(img) max(img)], 2);
  SO.img(imgno).range = tmap_range(i,:);
  
%  s_t = spm_input('Type of overlay', '+1', 'Transparent|Split', [1 2], 1); 
  if s_t == 1
%    SO.img(imgno).prop = spm_input('Colormap intensity', '+1', 'e', ...
%			       remcol/noaimgs,...
%			       1);
    SO.img(imgno).prop = col_intens;
    remcol = remcol - SO.img(imgno).prop;
  else  % split colormap
    if xor(SO.img(imgno).range(1) < SO.img(imgno).range(2), ...
	  SO.img(imgno).range(2) < 0)
      SO.img(imgno).outofrange = {[],size(map,1)};
    else
      SO.img(imgno).outofrange={[1], []};
    end
    SO.img(imgno).prop = Inf;
  end
  
%  if spm_input('Colour bar for this image?', '+1', 'Yes|No', [1 0], 1)
  if ~isempty(deblank(col_bar(i,:)))
    SO.cbar = [SO.cbar imgno]
    SO.scbar  = col_bar;
  end
  
  SO.img(imgno).nancol = [0 0 0];
end


% image orientation
% ------------------------------------------------------------------------------------
%ts = [0 0 0  0    0   0     -1 1 1;...
%      0 0 0 pi/2  0   0     -1 1 1;...
%      0 0 0 pi/2  0  -pi/2  -1 1 1];
ts = [  0  0  0   0    0   0     -1  1  1;...
        0  0  0  pi/2  0   0     -1  1  1;...
      -35  0  0  pi/2  0  -pi/2  -1  1  1];
% ------------------------------------------------------------------------------------

vxsz = 2.73;
%bbx = [-90 91 -126 91; -90 91 -50 85; -126 91 -50 85]; % MNI template
bbx = [-175 175 -175 175];
%t = spm_input('Image orientation', '+1', 'Axial|Coronal|Sagittal',...
%			 [1 2 3], 1);
SO.transform = spm_matrix(ts(t,:));
SO.orient = t;
SO.slicedef = [bbx(t,1) vxsz bbx(t,2); bbx(t,3) vxsz bbx(t,4)];
SO.slicedef

if t==1
  direc = 'TRA';
  lx = [lx:dx:hx];
elseif t==2
  direc = 'COR';
  lx = [-lx:-dx:-hx];
elseif t==3
  direc = 'SAG';
  lx = [hx:-dx:lx];
end


% slices for display
% ------------------------------------------------------------------------------------
%lx = spm_input('Lowest slice (mm)', '+1', 'e');
%if prod(size(lx))==1
%  hx = spm_input('Highest slice (mm)', '+1', 'e');
%  dx = spm_input('Slice separation (mm)', '+1', 'e', 4);

%lx = [hx:-dx:lx]
%end
SO.slices = lx;


% various plausible defaults
% ------------------------------------------------------------------------------------
SO.hold = 1; % resampling is trilinear.  See spm_sample_vol
SO.refreshf = 1; % always refresh figure window
SO.labels.colour = [1 1 1]; % color for slice labels
SO.labels.size   = 0.075; % font size normalized to image axis
SO.labels.format = '%+3.0f'; % format string for slice labels
if isempty(spm_figure('FindWin','fres'))
 spm_figure('Create','fres','fMRI Results','on');
end
SO.figure = spm_figure('FindWin','fres');  % use new SPM window
Fres = SO.figure;
%SO.figure = spm_figure('FindWin'); % use SPM figure window


% make description strings
% ------------------------------------------------------------------------------------
[a,datum] = unix('date');
datum = spm_str_manip(datum,'l30');
%datum = datum(40:69);
Descrip1 = name;
Descrip2 = spm_str_manip(CWD,'a50');
Descrip2 = str2mat(Descrip2,datum);
comment1 = spm_str_manip(comment1,'a53');
comment2 = spm_str_manip(comment2,'a53');
comment3 = spm_str_manip(comment3,'a53');
Descrip3 = str2mat(comment1,comment2,comment3);
SO.Descrip1 = Descrip1;
SO.Descrip2 = Descrip2;
SO.Descrip3 = Descrip3;


% and do the display
% ------------------------------------------------------------------------------------
f99_sl_overlay_long


% writing every thing to file
% ------------------------------------------------------------------------------------
% PRINTS the postscript file
if write_file
    str = ['cd ' dir_activ(1,:)];
    eval(str)
    psfile = [name '_' direc '.ps'];
    global PRINTSTR
    old_PRINTSTR=PRINTSTR;
    PRINTSTR = ['print -dpsc2 ' psfile];
    spm_print(Fres);
    % PRINTS the tif file
    fname = [name '_' direc];
    PRINTSTR = ['print(''-noui'',''-painters'',''-dtiff'',[''' fname '.tif'']);'];
    spm_print(Fres);
    PRINTSTR=old_PRINTSTR;
    str = ['cd ' CWD];
    eval(str)
end

