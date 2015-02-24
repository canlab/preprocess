function outputmatrix = outline(inputmatrix)
%OUTLINE outlines an array by returning only values that have at least one 0 around them.
% outline takes one input matrix of up to 3 dimensions that consists of 1's and 0's and
% returns a matrix of the same dimensions that has 1's only for border (at least one 0 adjacent)
% values.It does not consider the third dimension or the diagonal when making this comparison.



%Get size of inputmatrix

[inputrows,inputcols, inputdims]=size(inputmatrix);


%intialize output to input values, then remove appropriate 1's
outputmatrix=inputmatrix;


%cycle through all rows, cols, and dims
for currentdim=1:inputdims
	for lcv1=2:(inputcols-1)
		for lcv2=2:(inputrows-1)
			if inputmatrix(lcv2,lcv1,currentdim)~=0
				if (inputmatrix(lcv2-1,lcv1,currentdim)~=0 & inputmatrix(lcv2+1,lcv1,currentdim)~=0 & inputmatrix(lcv2,lcv1-1,currentdim)~=0 & inputmatrix(lcv2,lcv1+1,currentdim)~=0)
				
					outputmatrix(lcv2,lcv1,currentdim)=0;
				end
			end
		end
	end
end
