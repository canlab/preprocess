function bc = beta2conditions(b,DX,varargin)
% bc = beta2conditions(b,DX,varargin)
% DX is structure of info about design matrix - e.g., EXPT.DX
% created with tor_setup_deconvolution_design_ui.m or whatever
% optional NO SCALING parameter is varargin, can be anything

if ~isempty(DX.baseline) & length(varargin) == 0
    meanind = (length(b) - DX.nsess + 1):length(b);
    b = beta_scale(b,DX.baseparams,meanind);
end

for i = 1:length(DX.dxtrialonsets)
    
    bc{i} = b(DX.dxtrialonsets(i):DX.dxtrialonsets(i) + DX.numframes - 1);
    
end

return
