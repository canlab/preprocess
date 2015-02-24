function AVG = getslice(AVG)
% 
% to get average for a condition across all voxels, cond 1, slice 27:
% mean(AVG(27).avg{1}(~isnan(AVG(27).avg{1}(:,1)),:))
% **********************************

% Tor Wager, March 2001.  Last edit 3/08/01 to add the last feature above.

% set up arguments
wildcard = 'sra*img';
nimages = 1200;
threshold = 600; 
whichslices = 26:27;   

% for averaging
window = [-1 6];
subjevents = AVG.subjevents;
groups = [1 2 1 2 0 0];
binres = .5;

numslices = length(whichslices);
numgroups = max(groups);
[time,flindex,bincenter] = continuousavg_gettimes(1:nimages,window,subjevents,groups,binres);

sliceindex = 0;
for i = whichslices
   t = cputime;
   
   % setup slice
   % ===============================================
   AVG(i).numvox = 0;  
   for m = 1:numgroups
		AVG(i).avg{m} = zeros(1,length(bincenter{m}));
   		AVG(i).ste{m} = zeros(1,length(bincenter{m}));
	end
   sliceindex = sliceindex + 1;
   
   % get the slice ts
   % ===============================================
   t1 = cputime;
   slice = timeseries('slice',wildcard,nimages,i);
   disp(['Getting the timeseries for this slice took ' num2str(cputime -t1)]);
   
   
   %format the slice
   % ===============================================
   slice = squeeze(slice);				% dimensions should now be: x y time
   
   % The following line should perhaps be changed in the interests of speed - Ned
   slice = permute(slice,[3 1 2]);		% make time dim first - easy to work with
   
   eval(['save slice' num2str(i) ' slice'])
   
   index = 0;								% indexes which voxel it is in 1D array
   
   for j = 1:size(slice,2)				% j indexes rows of slice
      for k = 1:size(slice,3)			% k indexes cols of slice
         
      		index = index + 1;
      		timeseries = slice(:,j,k);
            
                        
      		if mean(timeseries) > threshold
               
            	% get the continuous average
             % ===============================================
             t1 = cputime;
             [myavg,myste] = ...
                continuousavg_slave(timeseries,window,groups,binres,time,flindex);
             for m = 1:numgroups
                try
                	AVG(i).avg{m}(index,:) = myavg{m};
                   AVG(i).ste{m}(index,:) = myste{m};
                catch
                   warning(['Wrong # of timepoints! '])
                   myavg, myste
                   AVG(i).avg{m}(index,:) = Inf/Inf;
                	AVG(i).ste{m}(index,:) = Inf/Inf;
					end
             end
            
             AVG(i).numvox = AVG(i).numvox+1;
            
            % do the F-test - store in slice array
            % ===============================================
         
         
         else
            for m = 1:numgroups
                AVG(i).avg{m}(index,:) = Inf/Inf;
                AVG(i).ste{m}(index,:) = Inf/Inf;
             end
               %fprintf('%3.2f ',mean(timeseries))
               % cont avg and F-test are null values.
      		end
      		disp(['Averaging for voxel ' num2str(index) ' took ' num2str(cputime - t1)]);

   		end % k loop
	end % j loop

   % write a mat file for the avgs and a slice file for the F value
   % ===============================================
   
   
   disp(['Slice ' num2str(i) ': ' num2str(sliceindex) ' of ' num2str(numslices) ': ' num2str(cputime - t) ' s. ' num2str(index) ' voxels.'])
   
	clear slice   
   
end % slices loop


% write an analyze img file of the F-map
% ===============================================

% use write_img_data



return


