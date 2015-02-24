


function results = match_som2motor(PMotor,PSOM)

motorVol = spm_read_vols(spm_vol(PMotor));

iMotor = find(motorVol);

results = [];

for iP = 1:size(PSOM,1)
    somVol = spm_read_vols(spm_vol(PSOM(iP,:)));
    results = [results sum(somVol(iMotor))];
end

return

%
% all done.
%