function [DX,delta] = physio_make_DX(ons,scanspersess,tp)
% [DX,delta] = physio_make_DX(ons,scanspersess)
%
% build deconvolution matrix for estimating evoked responses
% DX = deconv mtx
%
% ons = onsets in seconds, cell array, 1 cell per condition (collapse
% across session times)
%
% scanspersess = images per session, in TRs
%
% see script for fixed study-specific parameters!!!

% make deconv matrix in 1 s time bins
% 
samprate = 100; % sampling rate of physio
TR = 2;         % rep time of scanning, sec per image

scanspersess = [192 196 196 184 190 192] .* TR;     % in seconds  % s per session; length is num sessions*TR

% sampling rate, seconds per sample of physio
desired_samprate = 1;     % 1 = leave in seconds; 2, downsample by 2, etc.; delta is in time-bins of desired_samprate sec


len = sum(scanspersess); % length in s

if(~exist('tp') || isempty(tp))
    tp = 20;
end
 
% assumes ons is measured in seconds
[X,delta,delta_hires,hrf] = onsets2delta(ons,desired_samprate,len);
 
[DX,sf] = tor_make_deconv_mtx3(delta,tp,1,0,1,0,scanspersess);

return

