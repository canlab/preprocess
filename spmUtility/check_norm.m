% for normalized contrast images

if ~(exist('subs') == 1),
   subs = 	input('Enter subject numbers [x x x]: ');
end
if ~(exist('model') == 1),
   model = 	input('Enter model number: ','s');
end
if ~(exist('fmriDIR') == 1),
   fmriDIR = 	input('Enter base directory (fmriDIR): ','s');
end

% For Jonides 7
P = str2mat(['/data3/biman/spm99/spm99/templates/T1.img'],['/data3/biman/spm99/spm99/templates/EPI.img']);

% For Picasso
P = str2mat(['/usr/private/spm99/templates/T1.img'],['/usr/private/spm99/templates/EPI.img']);

   P2 = P;
   
	for JJ = 1:length(subs)
   		snum = num2str(subs(JJ));
   		%dnum = num2str(drive(JJ));
   
      
      datapath = [fmriDIR '/RESULTS/model' model '/sub' snum];
      % datapath = [studypath filesep 'sub' num2str(snum) filesep 'Anatomy'];
		P = strvcat(P,str2mat([datapath filesep 'sncon_0002.img']));
		P2 = strvcat(P,str2mat([datapath filesep 'sncon_0003.img']));
      
      
   end
   
spm_check_registration(P);
spm_orthviews('Interp',0);
str = input('Position crosshairs and then input text title for this page: ','s');
gtext(str,'FontSize',14,'Color','b')

spm_print

spm_check_registration(P2);
spm_orthviews('Interp',0);
str = input('Position crosshairs and then input text title for this page: ','s');
gtext(str,'FontSize',14,'Color','b')

spm_print