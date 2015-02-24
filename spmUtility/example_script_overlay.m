% '===================================================================================' 
%                   'Overlay of activations on a structural image' 
% '===================================================================================' 

%*************************************************************************** 
% Obligatory first part of the automation script                           % 
%*************************************************************************** 
%  fmriTEST = 1 to run the simulation mode 
%           = 0 to run the actual post-processing 
%  fmriDIR  = the directory in which all your data is situated 
%--------------------------------------------------------------------------- 
global fmriDIR 
fmriDIR = ['/images_spike1/beatse/sesRFX']; 
  

% function f99_display_slices 
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
f99_display_slices(... 
    ['anato'],... 
    ['nmpr_1.img'],... 
    ['5_dim_scal';... 
     '4_dim_scal'],... 
    ['sf_filtered_masked.img';... 
     'sd_filtered.img       '],... 
    0.75,2,... 
    ['col1  ';... 
     'col2  '],... 
    [3.21 12 ;... 
     3.21  7 ],... 
    2,0.6, 
    ['Hand-Rest';'Feet-Rest'],... 
    1,... 
    -50,86,4,... 
    'SID_FIX',... 
    'sesRFX',... 
    'thresh: 0.001/20',... 
    'non_corr',... 
    1) 
  


