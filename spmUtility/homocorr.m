function vol_homocor2(vol_path,varargin) 
% -----------------------------------------------------------
%       to correct for image nonuniformity
%       first fills in holes, then generates a low frequency 
%       background image and makes correction factor
%       input image is stored in img, output is Iout
%
%       homocor.m  rev 0    4/8/00                Gary Glover
%       ----------------------------------------------------- 
%
%   Modified to deal with SPM99 analize format image files.
%   Produces a 'homocor.img' volume in the working directory.
%
%   Works with volumes acquired in any orientation, The inplane
%   is taken to be the plane where two of the volume dimensions
%   match. If all dimensions match, it will default to axial.
% -----------------------------------------------------------
%   @(#)vol_homocor.m    1.10  Kalina Christoff    2000-05-29
%
%   Modified 9/11/01 by Tor Wager to take input file names.
%   Usage:
%   vol_homocor2(vol_path,outfile)
%      vol_path = path + input volume name or just file name (finds path)
%      outfile = name of output file to save
%           [optional: default is homocor.img]
%

%
% get the path
% ------------
if isempty(vol_path),
  vol_path=spm_get(1, 'img', 'Please select image to correct');
else
  vol_p=which(vol_path);
  if exist(vol_p) == 2
	vol_path = vol_p;
  else % else leave it alone, as input - contains relative path name
  end

end

% load the volume
% ---------------
if ischar(vol_path)
  vol = spm_vol(vol_path);
end

% make a guess as to what orientation are the slices
% eg. vol.dim=[256 256 124 4] probably axial, etc
% --------------------------------------------------

if       vol.dim(1)==vol.dim(2) % probably axial
            inpl_crd=[1 2]; thrpl_crd=3;
            orientation='axial';

 elseif  vol.dim(2)==vol.dim(3) % probably sagittal
            inpl_crd=[2 3]; thrpl_crd=1;
            orientation='sagittal';

 elseif  vol.dim(1)==vol.dim(3) % probably coronal
            inpl_crd=[1 3]; thrpl_crd=2;
            orientation='coronal'; 
end


  fprintf('\nAssuming the planes were acquired in ');
  fprintf(orientation); fprintf( ' orientation.\n\n');
      
foo = vol.dim(inpl_crd);
npix = foo(1);
np2 = npix/2;

%thres = input('gimme threshold % [5] =');
%if(isempty(thres))
  thres = 5;
%end;
thres = thres*.01;

Y = spm_read_vols(vol,0);

% initialize the corrected output volume Yc
Yc=zeros(vol.dim(1:3));

fprintf('Please wait - now working on plane    ');

% begin looping through slices
% ----------------------------
for slice = 1:vol.dim(thrpl_crd);
     
    if slice<10; fprintf('\b%d',slice); 
     elseif slice<100,  fprintf('\b\b%d',slice); 
     elseif slice<1000, fprintf('\b\b\b%d',slice);
    end

    % surprisingly, Y needs to be squeezed to remove the singleton 
    % dimensions for coronal or sagittal... interesting command.
    % ---------------------------------------------------------------
    if thrpl_crd==1; img=squeeze(Y(slice,:,:)); end;  
    if thrpl_crd==2; img=squeeze(Y(:,slice,:)); end;
    if thrpl_crd==3; img=squeeze(Y(:,:,slice)); end; % though no need here

    Imax = max(max(img));
    Ithr = Imax*thres;

    %  fill in holes

    mask = (img>Ithr);
    Imask = img.*mask;
    Iave = sum(sum(Imask))/sum(sum(mask));
    Ifill = Imask + (1-mask).*Iave;

    %  make low freq image

    z = fftshift(fft2(Ifill));

    %ndef = 3;
    ndef = np2/8;
    n = ndef;

    y2 = (1:np2).^2;  % use 15 as default win for guassian
    a2 = 1/n;
    for x=1:np2
      r2 = y2 +x*x;
      win1(x,:) = exp(-a2*r2);
    end

    win(np2+1:npix,np2+1:npix) = win1;
    win(np2+1:npix,np2:-1:1) = win1;
    win(1:np2,:) = win(npix:-1:np2+1,:);
    zl = win.*z;

    Ilow = abs(ifft2(fftshift(zl)));

    %corrected image

    Ilave = sum(sum(Ilow.*mask))/sum(sum(mask));
    Icorf = (Ilave./Ilow).*mask;  % correction factor
    Iout = img.*Icorf;

    % update the volume with the corrected values
    % -------------------------------------------
    if thrpl_crd==1; Yc(slice,:,:) = Iout; end
    if thrpl_crd==2; Yc(:,slice,:) = Iout; end
    if thrpl_crd==3; Yc(:,:,slice) = Iout; end
  
    clear img;
  
% end looping through slices
% ---------------------------
end

    fprintf('\n\nDone.\n\n')

% write corrected image (homocor.img)
% ==================================

Yc_vol = vol;
if nargin > 1
   Yc_vol.fname = varargin{1};
else
   Yc_vol.fname   ='homocor';
end
Yc_vol.descrip ='homocor';

spm_write_vol(Yc_vol,Yc);

