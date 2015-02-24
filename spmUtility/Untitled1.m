   P = str2mat(['/data3/biman/spm99/spm99/templates/T1.img']);
   P2 = P;
   
	for JJ = 1:length(subs)
   		snum = num2str(subs(JJ));
   		dnum = num2str(drive(JJ));
   
      
      datapath = [fmriDIR '/RESULTS/model' model '/sub' snum];
      % datapath = [studypath filesep 'sub' num2str(snum) filesep 'Anatomy'];
		P = strvcat(P,str2mat([datapath filesep 'sncon_0002.img']));
		P2 = strvcat(P,str2mat([datapath filesep 'sncon_0003.img']));

	end

spm_check_registration(P);
spm_orthviews('Interp',0);
spm_print
spm_check_registration(P2);
spm_orthviews('Interp',0);
spm_print
