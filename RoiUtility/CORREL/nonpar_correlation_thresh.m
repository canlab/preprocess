function [crit_r,cor_cdf,x,my_corr] = nonpar_correlation_thresh(my_len,crit_p)
% function [crit_r,cor_cdf,x,my_corr] = nonpar_correlation_thresh(my_len,crit_p)
%
% Tor Wager, 10/29/01

	disp(['Generating nonparametric null distribution of correlation values'])
	disp(['Using length of timeseries from first cluster set'])
	disp(['________________________________________________________'])
	for i = 1:8000
		myc = corrcoef(rand(my_len,1),rand(my_len,1));
		my_corr(i) = myc(1,2);
	end

	% --------------------------------------------------------------------------
	% * generate cumulative density function and get critical r value
	% --------------------------------------------------------------------------
	x = [.1:.01:1];
	index = 1;
	for i = x             
		cor_cdf(index) = sum(my_corr > i) ./ length(my_corr);
		index = index + 1;
	end
	crit_r = min(x(find(cor_cdf < crit_p)));

return