function bool=jsearch(array,srch)
% jsearch returns 1 if coordinates in search found in array, otherwise 0.
%Jsearch searches "array", an n-by-3 list of RxCxD coordinates for the coordinates given by "search". It returns
%1 if the search coordinates are found, 0 otherwise. This is a blobfinder utility.

[rows, columns]=size(array);
if ndims(srch)~=2, error('Search array must be a vector...'); end
bool=0;
for lcv1=1:rows

	if array(lcv1,1)==srch(1,1)
		
		if array(lcv1,2)==srch(1,2)

			if array(lcv1,3)==srch(1,3)

				bool=1;
			end %dims
		end %cols
	end %rows
end %for