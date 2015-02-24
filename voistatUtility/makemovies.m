% movie(f2,M{i},1,4);
%addpath /afs/math.lsa.umich.edu/group/psygrad/torw/scripts/voistat
%cd /data2/intext/sub9/task
imgs = 'sra*img';
scandir = 'scan';
loadfiles = input('Load images from files? 1 for yes, 0 if movies are in M{i} in the workspace. ')

if loadfiles
	dir = input('Enter directory name in '' '' or return for current directory: ');
	imgs = input('Enter image wildcard in '' '' or return for sra*img: ');
	if isempty(imgs),imgs = 'sra*img';,end

	eval(['cd ' dir])
end

format compact; clear R
startscan = input('Enter starting scan number: ');
endscan = input('Enter ending scan number: ');
scandir = input('Enter scan subdirectory in '' '' or return for ''scan'': ');
disp(''); interactive = input('Interactive?  1 for yes, 0 for no. ');
if interactive
	f1 = figure; set(gcf,'Position',[631 134 560 420]);
   f2 = figure; set(gcf,'Position',[51 123 560 420]);
   drawnow
end


% loop through scans
% =========================================================

for i = startscan:endscan
   
   if loadfiles
      
   		disp(['starting scan ' num2str(i)]);
         	str = ['!ls ' scandir num2str(i) '/' imgs ' > file.txt'];
         % str = ['!ls ' imgs ' > file.txt'];
   		disp(str);eval(str)
   		filenames = textread('file.txt','%s');
   
   		% pre-pend scan directory onto file if it's missing
   		% =========================================================
   		%if ~strcmp(filenames{1}(1:4),'scan')
      %		for j = 1:size(filenames,1)
      %   		filenames{j} = ['scan' num2str(i) '/' filenames{j}];
      %		end
   		%end
   
   		% make the movie
   		% =========================================================
   		basename = [filenames{1}(1:end-8)]
   		[M{i},slicearray{i},Cube{i}] = imgmovie(filenames,20);
   		Cube{i} = squeeze(Cube{i});
   end
   
   if interactive
   		% plot the images and movie
   		% =========================================================
   		close
   		figure(f1)
   		plot(Cube{i}); drawnow
      
   		figure(f2)
   		imagesc(slicearray{i}(:,:,Cube{i} == min(Cube{i})));
   		title('minimum global timeseries value')
   		input('Press return for the max value');
   		imagesc(slicearray{i}(:,:,Cube{i} == max(Cube{i})));
   		title('maximum global timeseries value')
   		R = 20;
   		while ~(isempty(R))
   			R = input('See movie?  Enter frames/sec or return to quit. ');
  			if ~isempty(R), movie(f2,M{i},1,R);,end
   			%	else input('Press return for the next scan');
   		end
   		clf   
   end
      
end

if loadfiles,save Movies,end
