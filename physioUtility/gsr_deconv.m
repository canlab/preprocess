function [b,bcen,gsr,gsrf] = gsr_deconv(gsr,scanspersess,DX,tp)
% deconvolve GSR
% b = estimated evoked responses, rows are time (s), cols are conditions
% bcen =
% gsr = downsampled gsr
% gsrf = filtered gsr ready for evoked resp estimation
%
% e.g. gsr_deconv(physio.rea15.task_data.gsr)
%
% fixed: HPlen, samprate,desired_samprate
%
% example:
% [b,bcen,gsr,gsrf] =
% gsr_deconv(physio.(ni_code).task_data.gsr,scanspersess,DX);

HPlen = 120;
samprate = 100;
desired_samprate = 1;
TR = 2;

if(~exist('tp') || isempty(tp))
    tp = 20;
end

gsr = resample(gsr,1,samprate);

%len = min(length(gsr),size(DX,1));
len = sum(scanspersess * TR); % length in s

% pad gsr to be correct length
m = mean(gsr);
gsr = gsr - m;
gsr = pad(gsr,len - length(gsr));
gsr = gsr + m;
gsr = gsr(1:len);

time = 0:length(gsr)-1;

% inputs: data, samp. rate of data, HP filter in s, spersess
[gsrf,I,S] = hpfilter(gsr,desired_samprate,HPlen,scanspersess * TR);
gsrf = trimts(gsrf,3,[]);

PDX = pinv(DX);
b = PDX * gsr;

% number of conditions, gives error i not an integer
nconds = ((size(DX,2)-length(scanspersess)) ./ tp);

% remove intercept cols
% columns of b are evoked responses for each condition
% b is tp x nconds
b = reshape(b(1:end-(length(scanspersess)))',tp,nconds);

base = mean(b(1:2,:));
bcen = b - repmat(base,size(b,1),1);

return