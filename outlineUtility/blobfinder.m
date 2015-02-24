function array=blobfinder(oldarray, sizethresh)

%Blobfinder searches a t-mask for contiguous non-zero entries and zeros out  blobs that do not meet a given 
%threshold size. It also returns a cell matrix that contains the coords for all of the voxels in each blob, 
%where each cell is a blob.
% Usage: [OutputArray]=blobfinder(InputArray, Threshold) where Threshold is a single digit specifying the
% minimum blob size in contiguous voxels in output. 
[rows,cols,dims]=size(oldarray);

%intialize things...
blob=1; %the number used to index the cell array done
todo=[]; % array of voxels to be processed
done{1}=[]; %cell array, each cell contains the array of coordinates for one blob

%loop through all voxels
for loop1=1:rows
	for loop2=1:cols
		for loop3=1:dims
			%first check: must be nonzero, and not already on the done or todo lists
			go=0;
			if (oldarray(loop1,loop2,loop3)~=0)  
			go=1;
			end
			for n=1:blob-1
			 	if (jsearch(done{n}, [loop1, loop2, loop3])==1) 
					go=0;
				end
			end	
					if go==1
					%add it to the done list
					done{blob}=[loop1, loop2, loop3];
					%check around it for additional voxels
					%+/- 1 row
					if oldarray(loop1-1, loop2, loop3)~=0 
						todo=[todo; [loop1-1, loop2, loop3]];
					end
				
					if oldarray(loop1+1, loop2, loop3)~=0 
						todo=[todo; [loop1+1, loop2, loop3]];
					end
					%+/- 1 col
					if oldarray(loop1, loop2-1, loop3)~=0 
						todo=[todo; [loop1, loop2-1, loop3]];
					end

					if oldarray(loop1, loop2+1, loop3)~=0 
						todo=[todo; [loop1, loop2+1, loop3]];
					end
					%+/- 1 dim
					if loop3>1 
						if oldarray(loop1, loop2, loop3-1)~=0 
							todo=[todo; [loop1, loop2, loop3-1]];
						end
					end
		%extra if is to make sure indexes do not exceed dimensions
					if loop3<3 & dims>1
						if oldarray(loop1, loop2, loop3+1)~=0 
							todo=[todo; [loop1, loop2, loop3+1]];
						end
					end
			%keep checking voxels until the queue is empty, also covers singleton voxels 
					while ~isempty(todo)
						
						%add current working voxel to done list
						done{blob}=[done{blob}; todo(1,:)];
						
						%check around working voxel for others	
						if ((oldarray(todo(1,1)-1, todo(1,2),todo(1,3))~=0) & ~(jsearch(todo,[todo(1,1)-1, todo(1,2),todo(1,3)]))& ~(jsearch(done{blob},[todo(1,1)-1, todo(1,2),todo(1,3)]))==1) 
							todo=[todo; [todo(1,1)-1, todo(1,2),todo(1,3)]];
						end
				
						if ((oldarray(todo(1,1)+1, todo(1,2),todo(1,3))~=0) & ~(jsearch(todo,[todo(1,1)+1, todo(1,2),todo(1,3)])) & ~(jsearch(done{blob},[todo(1,1)+1, todo(1,2),todo(1,3)]))==1) 
							todo=[todo; [todo(1,1)+1, todo(1,2),todo(1,3)]];
						end
					
						if ((oldarray(todo(1,1), todo(1,2)-1,todo(1,3))~=0) & ~(jsearch(todo,[todo(1,1), todo(1,2)-1,todo(1,3)])) & ~(jsearch(done{blob},[todo(1,1), todo(1,2)-1,todo(1,3)]))==1)
							todo=[todo; [todo(1,1), todo(1,2)-1,todo(1,3)]];
						end
						if ((oldarray(todo(1,1), todo(1,2)+1,todo(1,3))~=0) & ~(jsearch(todo,[todo(1,1), todo(1,2)+1,todo(1,3)])) & ~(jsearch(done{blob},[todo(1,1), todo(1,2)+1,todo(1,3)]))==1)
						
							todo=[todo; [todo(1,1), todo(1,2)+1,todo(1,3)]];
						end
						
						if todo(1,3)>1 
							if ((oldarray(todo(1,1), todo(1,2),todo(1,3)-1)~=0) & ~(jsearch(todo,[todo(1,1), todo(1,2),todo(1,3)-1])) & ~(jsearch(done{blob},[todo(1,1), todo(1,2),todo(1,3)-1]))==1)
								todo=[todo; [todo(1,1), todo(1,2),todo(1,3)-1]];
							end
						end
						
						if (loop3<3 & dims>1 & todo(1,3)<3)
							
						if ((oldarray(todo(1,1), todo(1,2),todo(1,3)+1)~=0) & ~(jsearch(todo,[todo(1,1), todo(1,2),todo(1,3)+1])) & ~(jsearch(done{blob},[todo(1,1), todo(1,2),todo(1,3)+1]))==1)
						
								todo=[todo; [todo(1,1), todo(1,2),todo(1,3)+1]];
							end
						end	
					%done, remove working voxel from todo list
						todo(1,:)=[];
					end %ends while 
				blob=blob+1;
				end %ends the if checking for nonzero voxs
		end %ends dims
	end %ends cols
end %ends rows

%so I now have "done", which is a cell array of all of the blobs. I want to remove blobs that do not meet 
%the size threshold. 
list=1:blob-1;
list=fliplr(list);
for lcv=list
	[m,n]=size(done{lcv});
	if m<sizethresh
		done(lcv)=[];
	end
end

%Create the output matrix, starting with a zeros array and adding in each blob from the "done" cell array

output=zeros(rows,cols,dims);

numblobs=size(done,2);

for currblob=1:numblobs

	while (~isempty(done{currblob}))
		
		output(done{currblob}(1,1), done{currblob}(1,2), done{currblob}(1,3))=currblob;
		done{currblob}(1,:)=[];
	end
end
array=output;