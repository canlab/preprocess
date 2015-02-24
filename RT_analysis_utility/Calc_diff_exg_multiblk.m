function [sqerr,stop_phat] = Calc_diff_exg_multiblk(phatstop,phat_sigma,ssoffset,rt,truestop,fmode)
%
%[sqerr,stop_phat]
% phatstop is a vector of proposed stop-signal distribution parameters
%       e.g., mu, sigma, tau for ex-gaussian
%       to vary mu only, for example, pass in [phatstop phat(2) phat(3)]
% 
% ssoffset is time from onset of trial to onset of stop signal
%
% rt is a cell array of GO rts for block1, 2, 3, etc.
%       e.g., rt{1} is the vector of rts for block 1. rt{2} for block 2,
%       etc.
% 
% truestop is the actual (obs) prob. of stopping on each block
%   column vector of n elements, where n is # blocks
%
% define function handle as follows:
% f= @(phatstop) Calc_diff_exg_multiblk(phatstop,phat,ssoffset,rt,truestop);


nblocks = length(rt);

for i = 1:nblocks
    
    switch fmode
        case 'exgauss'
            stop_phat(i,1) = Calc_diff_exgauss(phatstop,phat_sigma,0,ssoffset(i),rt{i});
        case 'ssrt'
            stop_phat(i,1) = nansum(rt{i} > phatstop + ssoffset(i)) ./ length(rt{i});
            % phatstop is ssrt
    end
end

% get summed squared errors across all blocks
sqerr = sum((stop_phat - truestop).^2);  

return

