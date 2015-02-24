%Outline Mask and Stack
% 6/19/01 by James
% Basic design: use readim2 to get 3-d array of tvalues, then turn it into a mask based
% on t-threshold. Run said mask through outline to remove all interior voxels, then apply the mask to the original 
% array.Batch processes all subjects and stacks output images so they can be checked for overlap.
% added blobfinder functionality. Blobfinder threshold can be set to trim data. 
%determine which subX directories exist

%config
tThresh=2;
minblob=10; %minimum size, in contiguous voxels, for a blob to pass through blobfinder

%begin looping by subject
%for subject=1:9
%   sub=num2str(subject);
%   directory=strcat('sub',sub);
   %check for existence of subject directory, then switch to it
%   if exist(directory,'dir')==7
%      eval(['cd ' directory])
subject=1;      
      
%turn origArray into a imgmask based on a min tvalue
[origArray,hdr,haxes]=readim2('spmT_0003');
origArray=double(origArray);
mask=max(0, origArray-tThresh);
mask=~~mask;
disp('Finding regions...')
%run it through blobfinder
mask=blobfinder(mask,minblob);

%run it through outline
mask=outline(mask);
mask(:,:,3)
%apply the mask

maskedArray=mask.*origArray;

%get size vector
maskSize=size(mask);

%arbitrary decision: make subplots 4 cols by n to plot all slices
sp=ceil(maskSize(3)/4);

subplot(sp,4,1);

%for-loop through each page of 3d and subplot on individual axes.
slice=1;
%format color and style of output such that each subject looks different

switch subject
      case 1
         color='ws';
      case 2
         color='ws';
      case 3
         color='ws';
      case 4
         color='ws';
      case 5
         color='ws';
      case 6
         color='ws';
      case 7
         color='ws';
      case 8
         color='ws';
      case 9
         color='ws';
      otherwise
         color='ws';
end

axis manual;
hold on

for slice=1:maskSize(3)
 	disp('Plotting slice')	
	subplot(sp,4,slice);
	axis ij
	axis([0 64 0 64]);
	for i=1:64
		for j=1:64  
			if mask(i,j,slice)~=0
 				number=num2str(mask(i,j,slice));   	
				text(i,j,number)
			end
		end
	end


	% plot(j,i,color)
end
 %  eval(['cd ..']); 
 %	end %end of fileexist if
%  end  %end of subject loop
