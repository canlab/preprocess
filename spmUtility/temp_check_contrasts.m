% special contrast checker for Biman4

load xCon.mat
load SPM.mat

e = sum(xX.X) == 0;


disp('Contrast 2')
a = [1 1 1 1 -1 -1 -1 -1 0 0];
b = tor_make_T_contrast_vector(a,Sess);
c = xCon(2).c';
sum(b - c(1:length(b)))
xCon(2).c(e)'

disp('Contrast 3')
a = [1 1 -1 -1 0 0 0 0 0 0];
b = tor_make_T_contrast_vector(a,Sess);
c = xCon(3).c';
sum(b - c(1:length(b)))
xCon(2).c(e)'

disp('Contrast 8')
a = [0 0 0 0 0 0 0 0 1 0]; 
b = tor_make_F_contrast_vector(a,Sess);
c = xCon(8).c';
sum(sum(b - c(:,1:length(b))))
xCon(2).c(e)'

