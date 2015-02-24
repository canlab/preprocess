function [hrf_est,b,DX] = tor_deconv(y,O,xX)

% O.tp = timepoints in estimate

myzeros = zeros(O.tp,1);
xlen = size(xX.sF{1},1);

% -------------------------------------------------------------------
% * make deconvolution matrix DX
% -------------------------------------------------------------------

index = 1;
for i = 1:size(xX.sF,2)
    DX(:,index) = xX.sF{i};
    index = index + 1;
    for j = 2:O.tp
        reg = [myzeros(1:j-1); xX.sF{i}];
        reg = reg(1:xlen);
        DX(:,index) = reg;
        index = index + 1;
    end
end

% -------------------------------------------------------------------
% * fit deconvolution matrix DX
% -------------------------------------------------------------------

% b = DX \ y;
b = pinv(DX) * y;

index = 1;
for i = 1:O.tp:length(b)
    hrf_est(:,index) = b(i:i+O.tp-1);
    index = index + 1;
end

return