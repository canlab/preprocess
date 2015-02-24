
function [patient_name,description,TR,bits_per_pixel,n_images,x_dim,y_dim,z_dim,x_vox_size,y_vox_size,z_vox_size,scale_slope]=PAR_read_Univ_v4pt2(infile)


%Eric Zarahn, July 2002
% November 2005 - Ajna modified to take care of the rescaling issue
%Ajna Borogovac, November 2005
%modified to correct for rescaling factor difference between dynamics (like in Gad_CBV data
% Now we correct for all files, be it SPGR or ROXY, or whatever.
% Ajna, May 2006 - Modified to accept PAR/REC format after the update (version 4).
% 26May2010, Ajna and Iris, changed the code to the 4pt2 for the upgrade.

fid=fopen(infile,'r');

%NAME
checker=1;
s1='';s2='';
while checker==1
s1=fscanf(fid,'%s',1);
if and(strcmp(s2,'Patient'),strcmp(s1,'name'))
checker=0;
fscanf(fid,'%s',1);
patient_name1=fscanf(fid,'%s',1);
patient_name2=fscanf(fid,'%s',1);
if ~(patient_name2=='.')
patient_name=strvcat(patient_name1,patient_name2);
else
patient_name=strvcat(patient_name1,'NO LAST NAME ENTERED');
end
end
s2=s1;
end


%DESCRIPTION
checker=1;
s1='';s2='';
while checker==1
s1=fscanf(fid,'%s',1);
if and(strcmp(s2,'Examination'),strcmp(s1,'name'))
checker=0;
fscanf(fid,'%s',1);
description=fscanf(fid,'%100c',1);
temp=findstr(description,'.');
description=description(1:temp(1)-1);
end
s2=s1;
end


%SCAN DURATION
checker=1;
s1='';s2='';
while checker==1
s1=fscanf(fid,'%s',1);
if and(strcmp(s2,'Scan'),strcmp(s1,'Duration'))
checker=0;
fscanf(fid,'%s',1);
scan_duration=fscanf(fid,'%d',1);
end
s2=s1;
end


%NUMBER OF SLICES
checker=1;
s1='';s2='';
while checker==1
s1=fscanf(fid,'%s',1);
if and(strcmp(s2,'of'),strcmp(s1,'slices/locations'))
checker=0;
fscanf(fid,'%s',1);
z_dim=fscanf(fid,'%d',1);
end
s2=s1;
end

%n_images
checker=1;
s1='';s2='';
while checker==1
s1=fscanf(fid,'%s',1);
if and(strcmp(s2,'of'),strcmp(s1,'dynamics'))
checker=0;
fscanf(fid,'%s',1);
n_images=fscanf(fid,'%d',1);
end
s2=s1;
end

%TR
checker=1;
s1='';s2='';
while checker==1
s1=fscanf(fid,'%s',1);
if and(strcmp(s2,'Repetition'),strcmp(s1,'time'))
checker=0;
fscanf(fid,'%s',2);
TR=fscanf(fid,'%f',1);
%n_images=1000*scan_duration/TR;
end
s2=s1;
end


%Ajna Added on May 2 2005 -- getting scale slope (ss) values
%Ajna and Iris changed to get the voxel size directly from numbers at the
%bottom of the PAR file/

checker=1;
s1='';
s2='';
fseek(fid,0,'bof');
while checker==1
    s1=fscanf(fid,'%s',1);
    if strcmp(s1,'diffusion')
        s2=fscanf(fid,'%s',1);
        if strcmp(s2,'L.ty')
        checker=0;
        [a,count]=fscanf(fid,'%f',inf);
        bits_per_pixel=a(8); %BITS PER PIXEL
        x_dim=a(10); %Reconstruction resolution in x dimension
        y_dim=a(11); %Reconstruction resolution in y dimension
            %VOXEL SIZE
            x_vox_size=a(29);
            y_vox_size=a(30);
  
            %x_vox_size=FOV_y_z_x(3,1)/x_dim;
            %y_vox_size=FOV_y_z_x(1,1)/y_dim;
            %z_dim=FOV_y_z_x(2,1)/z_vox_size;
        slice_thickness=a(23); %Slice thickness
        slice_gap=a(24); %Slice gap
        z_vox_size=slice_thickness+slice_gap;
    end
end
end

fclose(fid);
scale_slope(1)=a(14);

if n_images>=2
scale_slope(2)=a(length(a)-35);
end
