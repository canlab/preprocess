% Count runs of contiguous spikes
% Return index of each observation, whose value is the spike/exclude run
% length
% This is useful for identifying long runs of missing data that should be
% excluded
%
% y, data vector
% spikes, list of elements considered outliers/missing/etc
% spikelen = physio_missing_runlength(y, spikes)
%
% Tor Wager, June 2009

function spikelen = physio_missing_runlength(y, spikes)
    

isspike = zeros(size(y));
isspike(spikes) = 1;

isspike(1) = 0;  %  bug if first is spike

[cnt, tot, lenmat] = cnt_runs(isspike);

lenmat = lenmat(1:cnt);

indx = 1;
updateindx = 0;

spikelen = zeros(size(y));

% if spike, then replace with run length, and update flag = on for when
% run ends
% if no spike, check whether update is on.  If so, update indx to the next
% run, and turn off so we only update once before next run.

for i = 1:length(y)
    if isspike(i), spikelen(i) = lenmat(indx); updateindx = 1;  end
    
    if ~isspike(i) && updateindx, indx = indx + 1; updateindx = 0; end
    
end

end