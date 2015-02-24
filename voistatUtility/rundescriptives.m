subs = [10 11 12];
%drive = [1 2 4 4 4];

for JJ = subs
   load irfx
   snum = num2str(JJ);
   %dnum = num2str(drive(subs == JJ));
   eval(['cd sub' snum])
   !ls *ra*0001*img > file.txt
   fname = textread('file.txt','%s');
   fname = fname{1}(1:end-8)
   varname = ['s' snum 'dstats'];

   str =[varname ' = descriptives([],''' fname ''',200,X(:,2))']
   eval(str)
   
   str = [varname '.fname = ''' fname '''']
   eval(str)
   
   save mydescriptives
   
   str = (['Vi = getv(''make'',' varname '.avgxc)'])
   eval(str);
   save myVi Vi
   
   clear
   
   cd ..
   
end
