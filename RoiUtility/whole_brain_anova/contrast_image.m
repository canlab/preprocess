function Vo = contrast_image(P,Q,myc,dr)
% function Vo = contrast_image(P,Q,myc,dr)
%
% Tor Wager
%
% Creates a contrast image called Q (do not include path)
% Given a list of img files P (spm format, with path)
% and a contrast vector myc
% In the directory of 1st image in P.
%
% Edit: optional input to specify directory - Scott 10/31/13

V = spm_vol(P);
vols = spm_read_vols(V);

if ~(length(myc) == size(vols,4))
	error('Contrast vector length is not equal to number of image files.')
end

myc2 = zeros(size(vols));

for i = 1:length(myc) 
	myc2(:,:,:,i) = myc(i);
end

cvol = vols .* myc2;
cvol = sum(cvol,4);

% make sure row vec for name saving purposes
if length(myc) ~= size(myc, 2), myc = myc'; end

% -------------------------
% write
% -------------------------
if exist('dr','var')
    Q = fullfile(dr,Q);
else
    [d,f,e] = fileparts(V(1).fname);
    Q = fullfile(d,Q);
end
Vo = V(1);
Vo.fname = Q;
Vo.descrip = ['Contrast ' num2str(myc)];

spm_write_vol(Vo,cvol);
