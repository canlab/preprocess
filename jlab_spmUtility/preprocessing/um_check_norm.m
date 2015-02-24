function um_check_norm(OPT)
%
% OPT is as defined in example scripts
% Tor Wager 2/12/02

% for normalized contrast images

eval(['cd ' OPT.studydir filesep OPT.subjcode])

objectloc = [OPT.studydir filesep OPT.subjcode filesep 'anatomy' filesep OPT.object];
	%set objectloc so that we get to the correct directory where the object T1 is

first_nra = spm_get('Files','scan1','nra*img');
first_nra = first_nra(1,:);

last_nra = spm_get('Files',['scan' num2str(OPT.nruns)],'nra*img');
last_nra = last_nra(end,:);


first_snra = spm_get('Files','scan1','snra*img');
first_snra = first_snra(1,:);

last_snra = spm_get('Files',['scan' num2str(OPT.nruns)],'snra*img');
last_snra = last_snra(end,:);

P = str2mat(OPT.canonicalT1,objectloc,first_nra,last_nra,first_snra,last_snra);



spm_check_registration(P);
spm_orthviews('Interp',0);
disp('Position crosshairs and click on the image:');
gtext([OPT.subjcode ' template/spgr/first_nra/last_nra/first_snra/last_snra'],'FontSize',14,'Color','b')

spm_print

return