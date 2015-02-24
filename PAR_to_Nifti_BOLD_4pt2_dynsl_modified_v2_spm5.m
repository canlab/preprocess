% This code converts ASL Philips v4.1 files to nifti file format. The data
% in philips files must be listed dynamics first then slices.
% Ajna, Iris, August 2008, Partly based on Eric Zarahan's Philips to AVW
% code.
%
% Modified May 2011 by Yoni Ashar to support command line use
%   and automatically converting all files in a directory
% input: dir with PAR files, and array of file name patterns to exclude from conversion,
%     i.e.  ['SPGR', 'WIP']
%
% will search thru runsDir for .par and .PAR files and convert them
% skips any file that matches an excludePattern
%
% if folder of same name already exists, will NOT convert

function PAR_to_Nifti_BOLD_4pt2_dynsl_modified_v2_spm5(runsDir,excludePatterns)

warnings_state=warning('off','all');

%get listing of all .par and .PAR files
parFiles = dir([deblank(runsDir) '*.PAR']);
parFiles = cat(1, parFiles, dir([runsDir '*.par']));
num_files = size(parFiles,1);

single_flag=0;
for files=1:num_files


%STRING MANIPULATIONS
current_PAR_file= [runsDir  deblank(parFiles(files).name)];

%SKIP excludes
skip = false;
for tmp_cntr =1:size(excludePatterns, 2)
    if size(findstr(excludePatterns{tmp_cntr}, current_PAR_file),1) > 0 
        skip = true;
    end
end
if skip
   continue
end
    
temp=findstr('/',current_PAR_file);
current_path=current_PAR_file(1:temp(length(temp)));

temp=findstr('.',current_PAR_file);
prefix=current_PAR_file(1:temp(length(temp)));
if strcmp('PAR',current_PAR_file(length(current_PAR_file)-2:length(current_PAR_file)))
current_REC_file=[prefix 'REC'];
else
  current_REC_file=[prefix 'rec'];  
end

temp=findstr('/',prefix);
file_prefix=prefix(temp(length(temp))+1:length(prefix)-1);

%if folder w/ same name exists, skip
if exist(file_prefix, 'dir')
    continue
end

temp2=findstr('run',current_PAR_file);
temp=findstr('.',current_PAR_file);

functional_flag=0;


%READING PAR FILE (CONTAINS HEADER INFO)
[patient_name,description,TR,bits_per_pixel,n_images,x_dim,y_dim,z_dim,x_vox_size,y_vox_size,z_vox_size,scale_slope]=PAR_read_Univ_v4pt2(current_PAR_file);

if ~(strcmp(patient_name(2,:),'NO LAST NAME ENTERED'))
patient_name_save=strcat(deblank(patient_name(2,:)),'_',deblank(patient_name(1,:)));
else
patient_name_save=deblank(patient_name(1,:));
end

fprintf(1,'%s \n','');
fprintf(1,'%s \n',['There are ' num2str(z_dim) ' slices in this dataset']);

dim=[x_dim y_dim z_dim-single_flag*(z_dim-1) 4];
pinfo=zeros(3,1);
pinfo(1,1)=1;
affine=zeros(4,4);
affine(1,:)=[x_vox_size         0         0   -1*(x_dim)*x_vox_size/2];
affine(2,:)=[0         y_vox_size         0   -1*(y_dim)*y_vox_size/2];
affine(3,:)=[0         0         z_vox_size   -1*(z_dim+1)*z_vox_size/2];
affine(4,:)=[0         0         0  		  1            ];

AVW_hdr_info=struct('fname',[],'dim',dim,'mat',affine,'pinfo',pinfo,'descrip',strcat(description));

fid=fopen(current_REC_file,'r');

four_d_array_ui=fread(fid,inf,'*uint16');
four_d_array_ui=reshape(four_d_array_ui,x_dim,y_dim,n_images,z_dim);

%MAKING DIRECTORY FOR AVW FILES
temp=['cd ' current_path ];
eval(temp);

temp=['!mkdir ' file_prefix];
save_dir=[current_path file_prefix];

eval(temp);


% Added by Ajna on May 2005, image ordered to vary through dynamics first,
% then slices.

%WRITING AVW VOLUMES
    for rrep=1:n_images
    temp=double(squeeze(four_d_array_ui(:,:,rrep,:)));
    temp=temp./scale_slope(1);

%flipping images in y so that they end up in radiologic orientation
%(exactly match images from Mark Perera's conversion which we assume to 
%yield radiologic images, and which also yielded expected laterality of motor
%cortex activation; EZ 7/02)
temp=flipdim(temp,2);
temp=flipdim(temp,1);  %added by ajna for spm5

%joew        save_name=[save_dir '/' file_prefix '_BOLD_' num2str(rrep) '.nii'];
        save_name=sprintf('%s/%s_BOLD_%.3d.nii', save_dir, file_prefix, rrep);
tempdat=file_array;
tempdat.fname=save_name;
tempdat.dim=[x_dim y_dim z_dim];
tempdat.dtype='float32';
tempdat.offset=0;
tempdat.scl_slope=1;
tempdat.scl_inter=0;
tempnifti=nifti;
tempnifti.dat=tempdat;
tempnifti.mat=affine;
tempnifti.mat0=affine;
tempnifti.dat(1:x_dim,1:y_dim,1:z_dim)=temp(1:x_dim,1:y_dim,1:z_dim);
create(tempnifti)
clear tempdat tempnifti temp
    end
%end for REC files
end

end %main func
