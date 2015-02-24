sn = input('enter sub num: ')

eval(['cd /data4/intext2/sub' num2str(sn) '/task'])

mvto = ['/data4/intext2/RESULTS/model5/sub' num2str(sn) '/'];

for i = {'RPV*' 'SPM*' 'ResMS*' 'Y*' 'mask*' 'xCon*'}
   
   eval(['mv ' i{1} ' ' mvto])
   
end

   