vol1 = create_smoothed_pvals(64,64,16,4,0);

figure; hist(vol1(:),50)
% threshold the volume
vol1(vol1 >= .05) = NaN;
vol1(vol1 < .05) = 1;
vol1(isnan(vol1)) = 0;

totalvox = prod(size(vol1));

fp1 = sum(sum(sum(vol1)));
disp([num2str(fp1) ' (' num2str(100*fp1./totalvox) '%) false positives for vol1.'])

vol2 = create_smoothed_pvals(64,64,16,4,0);
figure; hist(vol2(:),50)
vol2(vol2 <= .95) = NaN;
vol2(vol2 > .95) = 1;
vol2(isnan(vol2)) = 0;

fp2 = sum(sum(sum(vol2)));
disp([num2str(fp2) ' (' num2str(100*fp2./totalvox) '%) false positives for vol2.'])

cvol = vol1 .* vol2;
fpc = sum(sum(sum(cvol)));

disp([num2str(fpc) ' (' num2str(100*fpc./totalvox) '%) false positives for conjunction.'])

expfp = prod(size(vol1)) .* .05 .* .05;
disp([num2str(expfp) ' (' num2str(100*expfp./totalvox) '%) false positives expected.'])

