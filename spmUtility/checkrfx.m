path('/data3/biman/spm99/diag',path) 

load SPM
if ~(exist('ResRMS.img') == 2)
   spm_imcalc_ui(['ResMS'],['ResRMS'],'sqrt(i1)');
end
spm_check_registration(str2mat(VY.fname,'spmT_0002','ResRMS'))
spm_orthviews('addcolorbar',1:11)
spm_orthviews('window',1:9,[-1.2 1.2])




%

%Diag_1

%diag_sptl
%diag_temp
%diag_spatemp
