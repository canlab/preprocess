function [DX,delta] = build_fir_model(ons,firpoints,spersess)
% function [DX,delta] = build_fir_model(tr,timepoints,spersess)
%
% Purpose: to take spm-style vectors of onsets, arranged in a cell array
% with one cell per condition nested within session, and build an FIR
% deconvolution model (DX)
%
% Inputs:
% tr, cell array of onset times (column vectors within cells), each cell is
% a condition, do all conditions in order and then replicate vector for
% each session
% 
% timepoints, a number or vector (k conditions) of time points to estimate
% with FIR model.
%
% spersess, a vector of how many images are in each session, e.g., [320
% 320] for 2 sessions of 320 images each
%
% see also: spm2dx (takes SPM.mat file)
% tor_make_deconv_mtx3.m  which does most of the work


[X,delta] = onsets2delta(ons,1,sum(spersess));

DX = tor_make_deconv_mtx3(delta,firpoints,1,0,1,0,spersess);

return

