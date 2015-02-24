function M = um_make_movie(imgs,names)

fprintf(1,'\num_make_movie.m\t')

what2do = input('Enter 1 for norm check, 2 for realign check ');


switch what2do

case 1

newi = imgs{1}(1,:);
for i = 2:length(imgs)
	newi = str2mat(newi,imgs{i}(1,:));
end

V = spm_vol(newi);
vols = spm_read_vols(V);

for sl = 1:size(vols,3)
	
	figure; set(gcf,'Color','k');colormap gray
	for j = 1:size(vols,4)

		imagesc(vols(:,:,sl,j)); axis image;
		
		axis off; title(names{j})

		text(.1,2,names{j},'Color','y','FontSize',18)
		text(.1,.1,['Slice ' num2str(sl)],'Color','y','FontSize',18)

		drawnow
		pause(.5)

	end
	close

end




case 2

sl = input('Input slice number: ');
myskip = input('View every nth slice - enter n: ');

for i = 1:length(imgs)

	%fprintf(1,'mapping...')
	%V = spm_vol(imgs{i});
	
	fprintf(1,'\nloading...')
	%vols = spm_read_vols(V);
	%vols = vols(:,:,sl,:);
	
	O.z = sl;
	O.include{1} = 1:myskip:size(imgs{i},1);
	vols = timeseries2('slice',imgs{i},O);

	figure; set(gcf,'Color','k');colormap gray
	M = moviein(size(vols,3));

	for j = 1:size(vols,3)

		imagesc(vols(:,:,j)); axis image;
		
		axis off; title(names{i})

		text(.1,2,names{i},'Color','y','FontSize',18)
		text(.1,.1,num2str(O.include{1}(j)),'Color','y','FontSize',18)

		drawnow

		M{i}(:,j) = getframe(gcf);
	end
	
	clear vols
	close
	%disp(['Showing movie ' names{i}])
	%movie(M{i})
	%close

end

end	% end switch


return
