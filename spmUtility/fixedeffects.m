function FEstats = fixedeffects(funct,path1,drivelist,path2,sslist,path3,basename)
% FEstats = fixedeffects(funct,path1,drivelist,path2,sslist,path3,basename)
% 1st argument is function: 'copy', 'compute', or 'both'
%
%
% stats = fixedeffects('both','/data',[1 1 1 1 1 2 2 4 4 4],'/intext/sub',[1 2 4 5 6 8 9 10 11 12],'/task/model6','sncon_0003');
% CON_0003 = fixedeffects('both','/data',[1 1 1 1 1 2 2 4 4 4],'/intext/sub',[1 2 4 5 6 8 9 10 11 12],'/task/model2','con_0003');
% Tstats = fixedeffects('/data','/intext/sub','/task/model6',[1 2 4 5 6 8 9 10 11 12;1 1 1 1 1 2 2 4 4 4],'snspmT_0003');
% - use this to compute mean T1 image the same way.
% - to not copy files, use [] for all 3 paths.
% - leave some paths [] if you need to.

if strcmp(funct,'copy') | strcmp(funct,'both')
% copy files over
% =================================================================================================
% build file list
disp(['current directory is ' pwd])
disp(['copying .img and .hdr files for subjects ' num2str(sslist(1,:)) ])

wild = 0;
for i = 1:size(basename,2)
   if strcmp(basename(i),'*'), wild = 1;,end
end

origbasename = basename;

for i = 1:size(sslist,2)
   ss = num2str(sslist(1,i));
   disp(['	...Subject ' ss])
   
   datadrive = num2str(drivelist(1,i));
   
   basename = [path1 datadrive path2 ss path3 filesep origbasename]
   
   if sslist(1,i) < 10, ssout = ['000' ss];,else ssout = ['00' ss];,end
   
	if wild
   		P = getfiles(basename)
	else P = {basename}
	end

   for j = 1:size(P,1)
      myfile = P{j};
      start = 1;
      for k = 1:length(myfile)
         if strcmp(myfile(k),'/'), start = k; end
      end
      justbase = myfile(start+1:end-4);
      extension = myfile(end-3:end);
      
   		str = ['!cp ' P{j}(1:end-4) extension ' ./' justbase '_' ssout extension]; disp([str])
   		eval(str);
      %str = ['!cp ' P{j}(1:end-4) '.hdr ./' justbase '_' ssout '.hdr'];
      %eval(str);
      %str = ['!cp ' P{j}(1:end-4) '.mat ./' justbase '_' ssout '.mat'];
      %eval(str);
	end
% end
end

end % if copy or both


if strcmp(funct,'compute') | strcmp(funct,'both')
% compute fixed effects map
% =================================================================================================  
FEstats.basename = basename;
FEstats.path = [path1 '*' path2 '*' path3];
FEstats = descriptives([],[basename '_'],sslist,[],-5000,'short');

FEstats.ss = sslist;
FEstats.N = size(sslist,2);
FEstats.fixedtmap = FEstats.mean .* sqrt(size(sslist,2));
FEstats.randtmap = FEstats.mean ./ (FEstats.std / sqrt(size(sslist,2)));

end

return
