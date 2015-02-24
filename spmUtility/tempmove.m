%sn = input('enter sub num: ')

for j = 1 % [2 4 6 8 9 10 11 12]
   sn = j

eval(['cd /data4/intext2/sub' num2str(sn) '/task'])

mvto = ['/data4/intext2/RESULTS/model5/sub' num2str(sn) '/'];

for i = {'beta*'} % {'RPV*' 'SPM*' 'ResMS*' 'Y*' 'mask*' 'xCon*'}
   
   eval(['!mv ' i{1} ' ' mvto])
   
end

end   