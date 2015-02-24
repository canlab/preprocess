function c = tor_make_T_contrast_vector2(a,Sess,whichSessions)
% a is the contrast weights for one scan
% c is the normalized contrast vector
% Sess is session struct
%
% use this version to make contrasts for only half the sessions!
% use version 1 to make contrasts over all sessions.
%
% useful if some columns in design are empty!
% tor wager, 10/18/01


if ~(sum(a) == 0),error('ERROR IN CONTRAST: DOES NOT SUM TO 0'),end
   
c = [];

mynull = zeros(size(a));

for j=1:length(Sess),

	if any(whichSessions == j)
		b = a;
	
		if length(b) ~= length(Sess{j}.ind)
			error(['Length of input string for session does not match number of regressors in session!'])
		end

		% Eliminate contrast weights with no values!!
	
		for k = 1:length(b)
			if sum(Sess{j}.sf{k}) == 0	% then there's no events!!
				b(k) = 0;				% so zero the contrast weight.
			end
		end
	else
		b = mynull;
	end
	c=[c b];
end     
	

	
% Now normalize contrast weights so it's even and each sums to one.	
c(c > 0) = c(c > 0) ./ (sum(abs(c(c > 0))));
c(c < 0) = c(c < 0) ./ (sum(abs(c(c < 0)))); 	


if round(sum(abs(c(c > 0)))*1000) ~= round(sum(abs(c(c < 0)))*1000)
	disp(['Positive: ' num2str(sum(abs(c(c > 0))))])
	disp(['Negative: ' num2str(sum(abs(c(c < 0))))])
	error('Contrast normalization error: positive and negative values are not balanced.')
end

return
