function [b,bcen,hr,hrf] = gsr_deconv(hr,scanspersess,DX,tp)
% deconvolve GSR
% b = estimated evoked responses, rows are time (s), cols are conditions
% bcen =
% hr = downsampled hr
% hrf = filtered hr ready for evoked resp estimation
%
% e.g. gsr_deconv(physio.rea15.task_data.hr)
%
% fixed: HPlen, samprate,desired_samprate
%
% example:
% [b,bcen,hr,hrf] =
% gsr_deconv(physio.(ni_code).task_data.hr,scanspersess,DX);

HPlen = 120;
samprate = 1;
desired_samprate = 1;
TR = 2;

if(~exist('tp') || isempty(tp))
    tp = 20;
end


% unnecessary, hr is already in s
%hr = resample(hr,1,samprate);

%len = min(length(hr),size(DX,1));
len = sum(scanspersess * TR); % length in s

% pad hr to be correct length
m = mean(hr);
hr = hr - m;
hr = pad(hr,len - length(hr));
hr = hr + m;
hr = hr(1:len);

time = 0:length(hr)-1;

% inputs: data, samp. rate of data, HP filter in s, spersess
[hrf,I,S] = hpfilter(hr,desired_samprate,HPlen,scanspersess * TR);
hrf = trimts(hrf,3,[]);

PDX = pinv(DX);
b = PDX * hr;

% number of conditions, gives error i not an integer
nconds = ((size(DX,2)-length(scanspersess)) ./ tp);

% columns of b are evoked responses for each condition
% b is tp x nconds
b = reshape(b(1:end-(length(scanspersess)))',tp,nconds);

base = mean(b(1:2,:));
bcen = b - repmat(base,size(b,1),1);

return