% Created 06/04/01 by Tor Wager


% -----------------------------
% original parameter settings
% -----------------------------

lowerbound = 0;



disp('nlfit.m')
disp('==============================================================')

disp('WARNING: warnings suppressed during this script')
disp(['Constraints: Lower bound is ' num2str(lowerbound)])
warning off

% -----------------------------
% string all trials together
% -----------------------------

eval(['ROI = load(''roi_' a ''');'])
eval(['ROI = ROI.roi_' a ';'])


% group 1
% =================================================================

	% -----------------------------
	% string all trials together
	% -----------------------------

	mygroup = 1;

	eval(['avgd = ROI.ntrials.group' num2str(mygroup) ';'])

	tpoints = size(avgd,2);
	trialind = (1:tpoints)';
	data = [];

	for i = 1:size(avgd,1)
        	data = [data; [trialind (avgd(i,:))']];
	end

	ROI.trialdata{mygroup} = data;

    
	% -----------------------------
	% fit nonlinear hrf
	% -----------------------------

	x = ROI.trialdata{mygroup}(:,1);
	y = ROI.trialdata{mygroup}(:,2);

	disp(['	nlfit:  fitting group ' num2str(mygroup)])

	myoptions = optimset('LevenbergMarquardt','on');
	[beta,resnorm,residual,exitflag,output,lambda,J]= lsqcurvefit('nlhrf',[1 6 1],x,y,lowerbound,Inf,myoptions);

	if exitflag < 0
		disp(['Fit for group ' num2str(mygroup) ' did not converge.  Using NaN for estimates.'])
		beta = beta * Inf / Inf;
		ci = ones(size(beta,2),2) * Inf / Inf;
	else
		disp(['	nlfit:  Converged on solution.  Computing confidence interval for group ' num2str(mygroup)])
		ci = nlparci(beta,residual,J); 
	end

	% if it doesn't fit, return nan values instead
	if sum(beta == [1 6 1]) == 3, beta = beta * Inf / Inf;, end

	ROI.params{mygroup} = beta;
	ROI.paramci{mygroup} = ci;
	ROI.fit{mygroup} = nlhrf(ROI.params{mygroup},(1:60)');

	% ---------------------------------
	% save data in first subtract index.
	% ---------------------------------

	ROI.subtractparams{1} = beta;
	ROI.subtractci{1}  = ci;
	ROI.subtractfit{1} = ROI.fit{1};




% group 2-1, 6-5, 11-10
% =================================================================

% ------------------------------
% loop through every other group
% ------------------------------
subtractindex = 2;
for mygroup = 1:2:6

	disp(['	nlfit: starting subtraction ' num2str(subtractindex)'])

	try
		ROI.params{mygroup}
		disp('above are computed params for first group.')
	catch
	
		% --------------------------------
		% prepare trial data for 1st group
		% --------------------------------

		eval(['avgd = ROI.ntrials.group' num2str(mygroup) ';'])
		tpoints = size(avgd,2);
		trialind = (1:tpoints)';
		data = [];
		for i = 1:size(avgd,1)
        		data = [data; [trialind (avgd(i,:))']];
		end
		ROI.trialdata{mygroup} = data;


                % --------------------------------
		% prepare fixed(convolution) response for 1st group
		% --------------------------------
		convwith = zeros(1,60);
		switch mygroup
		case 1
		   convwith([1]) = 1;
                case 3
                   convwith([1 3 5 7 9]) = 1;
                case 5
		    convwith([1 3 5 7 9 11 13 15 17 19]) = 1;
		otherwise
			error('Group number assignment is wrong!')
                end

		% -----------------------------
		% fit nonlinear hrf to 1st grp
		% -----------------------------

		x = ROI.trialdata{mygroup}(:,1);
		y = ROI.trialdata{mygroup}(:,2);

		disp(['	nlfit:  fitting group ' num2str(mygroup)])
		disp('----------------------------------------------------------------')

		myoptions = optimset('LevenbergMarquardt','on');
		[beta,resnorm,residual,exitflag,output,lambda,J]= lsqcurvefit('nlhrf',[1 6 1],x,y,lowerbound,Inf,myoptions,convwith);
	
		if exitflag < 0
			disp(['Fit for group ' num2str(mygroup) ' did not converge.  Using NaN for estimates.'])
			beta = beta * Inf / Inf;
			ci = ones(size(beta,2),2) * Inf / Inf;
		else
			disp(['	nlfit:  Converged on solution.  Computing confidence interval for group ' num2str(mygroup)])
			ci = nlparci(beta,residual,J); 
		end
		
		% if it doesn't fit, return nan values instead
		if sum(beta == [1 6 1]) == 3, beta = beta * Inf / Inf;, end

		ROI.params{mygroup} = beta;
		ROI.paramci{mygroup} = ci;
		ROI.fit{mygroup} = nlhrf(ROI.params{mygroup},(1:60)',convwith);

	end % catch ROI.params...

	% --------------------------------
	% prepare trial data for 2nd group
	% --------------------------------
	
	eval(['nextgroupavg = ROI.avgdata.group' num2str(mygroup+1) ';'])
	ROI.subtract{subtractindex} = nextgroupavg' - ROI.fit{mygroup};


	% --------------------------------
	% prepare delayed (convolution) response for subtracted waveform
	% --------------------------------
	convwith = zeros(1,60);
		switch mygroup+1
		case 2
		   convwith([3]) = 1;
                case 4
                   convwith([11]) = 1;
                case 6
		    convwith([21]) = 1;
		otherwise
			error('Group number for n+1 group is wrong!')
                end


	% --------------------------------
	% fit difference with nonlin hrf
	% --------------------------------

	x = (1:60)';
	y = ROI.subtract{subtractindex};

	[beta,resnorm,residual,exitflag,output,lambda,J]= lsqcurvefit('nlhrf',[1 6 1],x,y,lowerbound,Inf,myoptions,convwith);
	disp(['	nlfit:  confidence interval for group ' num2str(mygroup)])
	ci = nlparci(beta,residual,J); 

	% if it doesn't fit, return nan values instead
	if sum(beta == [1 6 1]) == 3, beta = beta * Inf / Inf;, end

	ROI.subtractparams{subtractindex} = beta;
	ROI.subtractci{subtractindex}  = ci;
	ROI.subtractfit{subtractindex} = nlhrf(ROI.subtractparams{subtractindex},(1:60)');

	subtractindex = subtractindex + 1;

end

warning on


