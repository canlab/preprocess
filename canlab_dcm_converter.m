%usage: canlab_dcm_converter(fileheader,struct_wildcard, basedir, Nseries, series_wildcard)
%author: Scott
%date: 5/17/2010
%purpose:   This takes a directory location of .dcm structural images and
%           converts them into .nii images by default.  The mean image is
%           saved in the base directory (default = pwd) and the dicom
%           images are moved to basedir/dicom_slices.
%
%           Nseries and series_wildcard are optional variables that are
%           used to specify series locations in the dicom images.  If
%           unspecified, they default to Nseries = 1 and series_wildcard = '';
%
%           All structurals are identified by
%           fullfile(basedir,[series_wildcard[series num]*,struct_wildcard]), so if the
%           dicom images are located within separate directories inside of
%           basedir, you need to include the path with the proper directory
%           syntax in your wildcard call (e.g. dir/dir for unix/mac or
%           dir\dir for pc)
%
%           series_wildcard is a single char array that describes the
%           beginning of a series before the number that defines it.  It is
%           expanded to a Nseries x columns char matrix where each row is
%           'wildcard1*'
%           'wildcard2*', etc.  If Nseries is 10 or greater, the output
%           looks like 'wildcard01*', etc.
%
%           If you want to convert functional .dcm files, you should simply
%           use somethings like the following code:
%               hdr = spm_dicom_headers(dcm_files);
%               spm_dicom_convert(hdr, 'all', 'flat', 'nii');
%           where dcm_files is a char matrix of .dcm filenames
%


function canlab_dcm_converter(fileheader,struct_wildcard, basedir, Nseries, series_wildcard)

if nargin < 3
    basedir = pwd;
end

if nargin < 4
    Nseries = 1;
end

if nargin < 5
    series_wildcard = '';
end

if Nseries > 1
    if isempty(series_wildcard)
        error('You have specified multiple series without definining them using series_wildcard');
    else
        if Nseries > 9
            swild = [series_wildcard, '01*'];
            swild = repmat(swild,Nseries,1);
            for n = 2:9
                swild(n,:) = [series_wildcard, '0', num2str(n), '*'];
            end
            for n = 10:Nseries
                swild(n,:) = [series_wildcard, num2str(n), '*'];
            end
        else
            swild = [series_wildcard, '1*'];
            swild = repmat(swild,Nseries,1);
            for n = 2:Nseries
                swild(n,:) = [series_wildcard, num2str(n), '*'];
            end
        end
    end
else
    swild = series_wildcard;
end



for j = 1:Nseries
    if isempty(swild)
        dcm_files = filenames(fullfile(basedir, struct_wildcard),'char', 'absolute');
    else
        dcm_files = filenames(fullfile(basedir, [swild(j,:), struct_wildcard]),'char', 'absolute');
    end
    n_in_series = size(dcm_files, 1);
    
    if n_in_series
        fprintf('Found %3.0f slices.\n', n_in_series);
        
        hdr = spm_dicom_headers(dcm_files);
        spm_dicom_convert(hdr, 'all', 'flat', 'nii');
        
        scan_wildcard = '*.nii';
        scan = filenames(scan_wildcard, 'char', 'absolute');
        movefile(scan, [fileheader '.nii']);
        
        slicedirname = fullfile(basedir,'dicom_slices');
        if ~exist(slicedirname, 'dir')
            mkdir(slicedirname);
%         else warning([slicedirname ' already exists!']) %#ok
        end
        for n = 1:n_in_series
            movefile(deblank(dcm_files(n,:)), slicedirname);
        end
    else
        error('No DICOM Series found.')
    end
end


nifti_imgs = filenames('*nii', 'char', 'absolute');
disp(nifti_imgs);
spm_check_registration(nifti_imgs)
