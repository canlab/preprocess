function c = tor_make_F_contrast_vector(a,Sess)
% a is the contrast weights for one scan, in a vector
% this function will put them on different lines in the F contrast.
% c is the normalized contrast vector
% Sess is session struct

% useful if some columns in design are empty!
% tor wager, 10/18/01


if any(a < 0),error('No negative values allowed in F contrast'),end
   
c = [];

for j=1:length(Sess),
	% ---------------------------------------------------------------
	% Work on normalizing b vector, copy of input vector a
	% ---------------------------------------------------------------.
	b = a;

	% ---------------------------------------------------------------
	% Define session contrast and column dims in session contrast
	% ---------------------------------------------------------------.
	sessC = [];
	numRows = sum(b > 0);	
	column = zeros(numRows,1);

	if length(b) ~= length(Sess{j}.ind)
		error(['Length of input string for session does not match number of regressors in session!'])
	end

	% ---------------------------------------------------------------
	% Eliminate contrast weights with no values!!
	% ---------------------------------------------------------------
	index = 1;
	for k = 1:length(b)
		if sum(Sess{j}.sf{k}) == 0	% then there's no events!!
			b(k) = 0;				% so zero the contrast weight.
		end
		
		% ---------------------------------------------------------------
		% make session F contrast
		% ---------------------------------------------------------------
		myColumn = column;
		if b(k) > 0,
			myColumn(index,1) = b(k);
			index = index + 1;
		end

		sessC = [sessC myColumn];
	end

	% ---------------------------------------------------------------
	% Add session C to overall contrast set c	
	% ---------------------------------------------------------------
	rowsToAdd = size(sessC,1);
	colsToAdd = size(sessC,2);
	
	if j == 1,
		c = sessC;
	else
		sessC = [zeros(size(c,1),colsToAdd); sessC];
		c = [c; zeros(rowsToAdd,size(c,2))];
		c = [c sessC];
	end

end     
	

% ---------------------------------------------------------------	
% Now normalize contrast weights so it's even and each sums to one.
% ---------------------------------------------------------------	
c(c > 0) = c(c > 0) ./ (sum(abs(c(c > 0))));	


if round(sum(abs(c(c > 0)))*1000) ~= 1000
	disp(['Positive: ' num2str(sum(abs(c(c > 0))))])
	error('Contrast normalization error: positive values are not normalized to 1.')
end

return
