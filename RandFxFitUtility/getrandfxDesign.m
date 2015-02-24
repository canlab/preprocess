function [model,mydelta] = getrandfxDesign(TR,boldreps,irf,varargin)
%function model = getrandfxDesign(TR,boldreps,irf,varargin)
% :::		:::		:::		:::		:::		:::		:::
% TR: in seconds
% boldreps: number of images in scan
% irf: filename for your subject's IRF, in a column vector,or 'spm' for spm's gamma function HRF
% varargin{1} - x : behData matrices for your scans.
% behData: onset times in col. 1, condition nums in col 2.
%
% getDesign5: sampling rate is every .01 seconds, NOT every ms as in getDesign4

if nargin < 4,error('no behavioral data vectors entered.'),end

%offset = 125 - 10000; 						% offset from onset time to actual stim presentation
model = [];
mydelta = [];

if not(strcmp(irf,'spm') | strcmp(irf,'SPM'))						% use IRF from input
   disp('building HRF...using custom function')
   disp(['loading ' irf])
   irf = load(irf);
   if (size(irf,1) == 400) & (TR == 2)
      disp(['using long irf']) 
      HRF = interp(irf,5);   
   elseif size(irf,1) == 10
      error('use the long irf - av. over 20 trials, and interpolate before saving the timeseries for intext.')
      % HRF = interp(irf,1000*TR);
   else error('IRF should be 400 elements.  In Voxbo, use trial averaging for 20 trials, and interpolate before saving.')
   end
   
   HRF = HRF / max(HRF);

else																				% use SPM's gamma function
   disp('building HRF...using spm_hrf gamma function')
   HRF = spm_hrf(.01);
   HRF = HRF / max(HRF);
end

for i = 1:nargin - 3
   disp(['starting scan ' num2str(i)])
	behData = varargin{i};   
   %behData(:,1) = behData(:,1) + offset;	% adjust to reflect real onset times

   nConds = size(behData,1);
   nms = ceil(max(behData(:,1))/10);					% no of ms in scan / 10: # of samples
   disp('	making delta function...')
	deltaF = zeros(nms,nConds);				% make Delta function in 10 ms res
	for i = 1:nConds
      deltaF(round(behData(i,1)/10),i) = 1;
   end
   if size(deltaF,1) > boldreps*TR*100,deltaF = deltaF(1:boldreps*TR*100,:);,end
   mydelta = [mydelta;deltaF];

	disp('	convolving regressors...')
	for i = 1:nConds
   	disp(['		starting ' num2str(i) ' / ' num2str(nConds)])
   	reg = conv(deltaF(:,i),HRF);
   	reg = resample(reg,1,TR*100);			% resample at TR
   	reg = reg(1:boldreps);					% stop at last boldrep in scan
   	scan(:,i) = reg;
   end
   model = [model;scan];
end

return