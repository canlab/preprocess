function [wiSet,btwSet,crit_r,mynames,rx,r_cdf,allts] = correl_cluster_ui(varargin)
% function [wiSet,btwSet,crit_r,mynames,rx,r_cdf,allts] = correl_cluster_ui(varargin)
%
% Optional inputs: 
% 	1st: user-entered predictors for correlation
%	2nd: names of user-entered predictors (not required)
%       this is a text string containing tab-separated names for fprintf output, e.g.,
%       'user1\tuser2\tuser3\t'
% Tor Wager, 10/30/01

% -----------------------------------------------------------------------------------------
% * Get filenames of clusters and user-specified predictors
% -----------------------------------------------------------------------------------------
	myOut = input('Enter output style (full, summary, or anything else for none): ','s');
	if strcmp(myOut,'full') | strcmp(myOut,'summary'),diary corr_results.txt,end

	P = spm_get([1 10],'*.mat','Select mat files containing cluster variables',pwd,0)	

	if strcmp(myOut,'summary'), diary off, end

	index = 1;

    usercnames = [];
	%userc = input('Are there user-entered behavioral predictors? (y or n)','s');
	%if strcmp(userc,'y')
    if length(varargin) > 0
		try
			userc = varargin{1};
			disp(['Found ' num2str(size(userc,2)) ' user covariates of length ' num2str(size(userc,1))]) 
		catch
			disp('Must enter user covariates as optional input argument. Skipping user-entered...')
		end

		if nargin > 1
			usercnames = varargin{2};
            usercnames_str = [];
            for i = 1:size(usercnames,1)
                usercnames_str = [usercnames_str usercnames{i} '\t'];
            end
                
        else
            usercnames = [];
            for i = 1:size(userc,2)
                usercnames{i} = ['user' num2str(i)];
                usercnames_str = [usercnames_str '\tuser' num2str(i)];
            end
		end
	else
		userc = [];
	end		

% -----------------------------------------------------------------------------------------
% * Within cluster set correlations for all cluster sets
% -----------------------------------------------------------------------------------------
    allnms = [];

	for i = 1:size(P,1)

		ts = [];
		load(deblank(P(i,:)),'clusters')

         % cluster names
         % -------------------------------------------------------------------------------
         a = '_'; a = repmat(a,length(clusters),1);
         cnames{i} = [num2str((1:length(clusters))') a num2str(cat(1,clusters.numVox))];
         nms = [];
         for j = 1:length(clusters)
             nms = [nms 's' num2str(i) '_' cnames{i}(j,:) '\t'];
         end
         allnms = [allnms nms]; % add to all names for btwn clusters
         nms = [nms usercnames_str];
         nms = [nms '\n'];
         

		disp(['Cluster set ' num2str(i) ', found ' num2str(length(clusters)) ' regions.'])
		disp(['________________________________________________________'])
        fprintf(1,nms)

		for j = 1:length(clusters)

			ts(:,j) = clusters(j).timeseries;

			mynames{index} = ['Set ' num2str(i) ' cluster ' num2str(j)];
			index = index + 1;
			
		end
		
        % save cluster ts for btwn-clusters
        % -------------------------------------------------------------------------------
        if i == 1 
			allts = ts;
		elseif (size(ts,1) == size(allts,1))
			allts = [allts ts];
		else 
			disp('Not all cluster timeseries are same length.  No correlations between sets.')
		end
        
        if ~isempty(userc)
            ts(:,end+1:end+size(userc,2)) = userc;
        end
        
		wiSet{i} = corrcoef(ts);
		wiSet{i} = tril(wiSet{i});
		wiSet{i}



	end

% -----------------------------------------------------------------------------------------
% * Correlate across clusters and within, overall
% -----------------------------------------------------------------------------------------

	disp(['All cluster sets'])
	disp(['________________________________________________________'])

	% user predictor stuff
	% --------------------------------
	if ~isempty(userc) & size(userc,1) == size(allts,1) 
		allts = [allts userc];
        allnms = [allnms usercnames_str];
		for i = 1:size(userc,2)
			if exist('usercnames') == 1
				mynames{end+1} = usercnames{i};
            else
				mynames{end+1} = ['User_pred ' num2str(i)];
            end
			%%%%index = index + 1; don't need
		end
	elseif ~isempty(userc)
		warning('User predictors different length that cluster timeseries!  ...user preds not added.')
		disp(['User length: ' num2str(size(userc,1)) ' timeseries length: ' num2str(size(allts,1))])
	end

    
	btwSet = corrcoef(allts);
	btwSet = tril(btwSet);
	btwSet(btwSet == 1) = 0;		% zero diagonals
	
    fprintf(1,allnms)
    btwSet
    
    means = mean(allts)
    stdevs = std(allts)

if strcmp(myOut,'summary'),diary corr_results.txt,end

% -----------------------------------------------------------------------------------------
% * Significance testing
% -----------------------------------------------------------------------------------------

	doSig = input('Estimate significance levels with nonparametric r distribution? (y or n)','s');
	if strcmp(doSig,'y')

		crit_p = input('Enter significance threshold (e.g., .05): ');

		% Bonferroni correction
		% ---------------------------------------------------------------------------------
		doBonf = input('Do bonferroni correction? (y or n)','s');
		if strcmp(doBonf,'y')
			total_corr = sum(sum(btwSet ~= 0));
			num_user = size(userc,2);
			num_region = size(btwSet,2) - num_user;
			user_corr = num_user * num_region; 
			user_user_corr = ((num_user .^2) - num_user) ./2 ;
			region_region_corr = ((num_region .^2) - num_region) ./2;
			disp(['Total unique correlations calculated: ' num2str(total_corr) ' on ' num2str(num_region) ' regions and ' num2str(num_user) ' user variables'])
			disp([num2str(region_region_corr) ' between regions, ' num2str(user_user_corr) ' between user vars, and ' num2str(user_corr) ' between user vars and regions.'])
			
			my_bonf = input('Enter num comparisons for bonferroni correction: ');
			crit_p = crit_p ./ my_bonf;
			disp(['Critical p value is ' num2str(crit_p) ' for ' num2str(my_bonf) ' comparisons.'])
		end

		% nonparametric estimation
		% ---------------------------------------------------------------------------------
		[crit_r,r_cdf,rx] = nonpar_correlation_thresh(size(ts,1),crit_p);
		disp(['Critical r value is ' num2str(crit_r)])

		[a,b] = find(abs(btwSet) > crit_r);
		c = btwSet(abs(btwSet) > crit_r);
        
		
		% -----------------------------------------------------------------------------------------
		% * Output of significant correlations
		% -----------------------------------------------------------------------------------------
		try
			for i = 1:size(a,1);
				myn1 = mynames{a(i)};	
				myn2 = mynames{b(i)};
				myp = min(r_cdf(round(rx*100) == round(abs(c(i)) * 100)));
				disp([myn1 ' with ' myn2 ': r = ' num2str(c(i)) ' p = ' num2str(myp)])
			end
		catch
			disp('Trouble printing results table')
			mynames
			whos a, whos b, whos i
			a(i)
		end
	
	end

	diary off
return
