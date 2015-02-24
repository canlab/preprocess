% This code converts Philips v4.2 structural files to nifti file format
% Iris Asllani and Ajna Borogovac, August 2008
% Code is based on Eric Zarahan's spm5_Philips_to_AVW
% Edited Dec 2009 to change the way it reads the voxel size so that it can
% be used for all 3 image orientation.
%
% Modified May 2011 by Yoni Ashar to support command line use
%   and automatically converting all files in a directory
% input: dir with PAR files, and 'y' or 'n' for whether to create PVEc
%     directory as well

function PAR_to_Nifti_Structural_univ_4pt2_spm5(runsDir, structFileIdentifier, PVEc)


warnings_state=warning('off','all');

disp('This code converts Philips v4.1 structural files to nifti file format')
disp('Iris Asllani and Ajna Borogovac, August 2008')
disp('Code is based on Eric Zarahans spm5_Philips_to_AVW')

%get listing of all .par and .PAR files
parFiles = dir([runsDir '*' structFileIdentifier '*.PAR']);
parFiles = cat(1, parFiles, dir([runsDir '*' structFileIdentifier '*.par']));
num_files = size(parFiles,1);


analysistype=PVEc;
for files=1:num_files

%STRING MANIPULATIONS
current_PAR_file= [runsDir  deblank(parFiles(files).name)];

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
anal_flag='_Conventional';
if exist([file_prefix anal_flag], 'dir')
    continue
end

temp2=findstr('run',current_PAR_file);
temp=findstr('.',current_PAR_file);
if ~(isempty(temp2))
current_exam_number=current_PAR_file(temp2(length(temp2))+3:temp(length(temp))-1);
functional_flag=1;
else
functional_flag=0;
end

%READING PAR FILE (CONTAINS HEADER INFO)
[patient_name,description,TR,bits_per_pixel,n_images,x_dim,y_dim,z_dim,x_vox_size,y_vox_size,z_vox_size,scale_slope]=PAR_read_Univ_v4pt2(current_PAR_file);

if ~(strcmp(patient_name(2,:),'NO LAST NAME ENTERED'))
patient_name_save=strcat(deblank(patient_name(2,:)),'_',deblank(patient_name(1,:)));
else
patient_name_save=deblank(patient_name(1,:));
end

dim=[x_dim y_dim z_dim 4];
pinfo=zeros(3,1);
pinfo(1,1)=1;
affine=zeros(4,4);
affine(1,:)=[x_vox_size         0         0   -1*(x_dim)*x_vox_size/2];
affine(2,:)=[0         y_vox_size         0   -1*(y_dim)*y_vox_size/2];
affine(3,:)=[0         0         z_vox_size   -1*(z_dim+1)*z_vox_size/2];
affine(4,:)=[0         0         0  		  1            ];

AVW_hdr_info=struct('fname',[],'dim',dim,'mat',affine,'pinfo',pinfo,'descrip',strcat(description));

four_d_array=zeros(x_dim,y_dim,z_dim,n_images);

fid=fopen(current_REC_file,'r');

temp=fread(fid,inf,'uint16');
temp_new=reshape(temp,x_dim,y_dim,z_dim,n_images);
if n_images==1
four_d_array=temp_new./scale_slope(1);
elseif n_images>=2
    if z_dim==1
        four_d_array(:,:,1:n_images/2,1)=temp_new(:,:,1:n_images/2,1)./scale_slope(1);
        four_d_array(:,:,1+(n_images/2):n_images,1)=temp_new(:,:,1+(n_images/2):n_images,1)./scale_slope(2);
    else
        four_d_array(:,:,:,1:n_images/2)=temp_new(:,:,:,1:n_images/2)./scale_slope(1);
        four_d_array(:,:,:,1+(n_images/2):n_images)=temp_new(:,:,:,1+(n_images/2):n_images)./scale_slope(2); 
    end
end

%flipping images in y so that they end up in radiologic orientation
%(exactly match images from Mark Perera's conversion which we assume to 
%yield radiologic images, and which also yielded expected laterality of motor
%cortex activation; EZ 7/02)
four_d_array=flipdim(four_d_array,2);
four_d_array=flipdim(four_d_array,1);  %also added this for spm5

%MAKING DIRECTORY FOR AVW FILES
temp=['cd ' current_path ];
eval(temp);



temp=['!mkdir ' file_prefix anal_flag];
save_dir=[current_path file_prefix anal_flag];

eval(temp);

%WRITING AVW VOLUMES

for rep=1:n_images
    if z_dim==1
       temp=squeeze(four_d_array(:,:,:,rep));
    else
    temp=squeeze(four_d_array(:,:,:,rep));
    end
save_name=[save_dir '/' file_prefix '_' num2str(rep) '.nii'];
tempdat=file_array;
tempdat.fname=save_name;
tempdat.dim=[x_dim y_dim z_dim];
tempdat.dtype='float64';
tempdat.offset=0;
tempdat.scl_slope=1;
tempdat.scl_inter=0;
tempnifti=nifti;
tempnifti.dat=tempdat;
tempnifti.mat=affine;
tempnifti.mat0=affine;
tempnifti.dat(1:x_dim,1:y_dim,1:z_dim)=temp(1:x_dim,1:y_dim,1:z_dim);
create(tempnifti)
end
if analysistype=='y'
    anal_flag='_PVEc';
    temp=['!mkdir ' file_prefix anal_flag];
save_dir=[current_path file_prefix anal_flag];

eval(temp);

%WRITING AVW VOLUMES

for rep=1:n_images
    if z_dim==1
       temp=squeeze(four_d_array(:,:,:,rep));
    else
    temp=squeeze(four_d_array(:,:,:,rep));
    end
save_name=[save_dir '/' file_prefix '_' num2str(rep) '.nii'];
tempdat=file_array;
tempdat.fname=save_name;
tempdat.dim=[x_dim y_dim z_dim];
tempdat.dtype='float64';
tempdat.offset=0;
tempdat.scl_slope=1;
tempdat.scl_inter=0;
tempnifti=nifti;
tempnifti.dat=tempdat;
tempnifti.mat=affine;
tempnifti.mat0=affine;
tempnifti.dat(1:x_dim,1:y_dim,1:z_dim)=temp(1:x_dim,1:y_dim,1:z_dim);
create(tempnifti)
end
end




%end for REC files
end

end %main function
