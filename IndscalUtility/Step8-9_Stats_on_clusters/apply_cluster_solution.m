function out = apply_cluster_solution(c,xc,varargin)
    % function out = apply_cluster_solution(c,xc,varargin)
    %
    % Purpose:   Test correlations within and between classes of objects
    %
    % Inputs:
    % c          class assignments for objects (e.g., regions) from clustering
    %            solution (DATA.indscal.classes)
    % xc         3-D correlation matrix, REGIONS x REGIONS x subjects, k x k x n
    %
    %            * Test either AVG of correlations across task conditons
    %            to test overall diffs between and within
    %            or CONTRAST (difference) across task conditions, to test for diffs
    %            between conditions in network coherence
    %            * To test a contrast, enter a 'states' vector that identifies
    %            the state membership of each 1...n in xc and a 'contrast'
    %            value.
    %
    % variable arguments are in 'fieldname',value pairs
    % 'doz',[1] / 0       covert to z-scores using Fisher's r-to-z, [1] default
    % 'names',{'n1' 'n2'} cell array of names for network classes (groups)
    % 'bcov',[]           n x 1 vector of between-subjects covariates
    %
    % 'states'            STATES IS NOT USED YET -- USES ONLY CONTRAST INPUT
    %                     AND ASSUMES CONTRAST-LENGTH STATES NESTED WITHIN SUBJECTS
    %                     vector of which 3rd dim elements in xc go with which
    %                     task states.  Use for applying a contrast across task
    %
    %                     states.  This function assumes states are WITHIN
    %                     subjects, and will compute contrasts within adjacent
    %                     matrices in xc -- e.g., if there are two states, then
    %                     the matrices along the 3rd dim of xc should be
    %                     subject1_state1 subject1_state2 sub2_state1
    %                     sub2_state2, etc., in that order.  States are nested
    %                     within subjects on the 3rd dim of xc.
    %
    % 'contrast',[1 -1]   contrast across task states (must enter vector coding
    %                     for task states as well), e.g., [1 -1] is the
    %                     difference between states 1 and 2.  All output
    %                     reflects the contrast value.  If no contrast or state
    %                     values is entered, the default contrast is 1 or [1
    %                     1 1...] across all states.
    %
    % for each subject, apply cluster indicator vector c to each correlation
    % matrix, save within and between correlation values, covert to z-scores,
    % and test across subjects.
    %
    % a within-cluster correlation is a measure of the functional coherence of
    % a network of regions.  a between-cluster correlation is a measure of the
    % functional relationship between groups of regions.
    %
    % if a between-subjects vector is entered, performs multiple regression on
    % the between- and within- cluster correlations, looking for the best
    % predictors of the between-subjects vector (e.g., behavioral scores)
    %
    % by tor wager, 12/5/04
    %
    % Functions called:
    % contrast3d
    % ttest3d
    % correlation_to_text
    % stepwisefit (stats toolbox)
    %
    % Examples:
    % -------------------------
    % Basic solution, no between-subjects predictors, one state variable:
    % out = apply_cluster_solution(DATA.oGs,DATA.xc);
    %
    % To apply an A - B contrast across tasks within subject, assuming that 3rd
    % dim of xc has tasks within subject (As does output of mvroi_xc DATA.xc)
    % out =
    % apply_cluster_solution(DATA.oGs,DATA.xc,'state',DATA.STATS.indscal.taskvector,'contrast',[1 -1]);


    % ------------------------------------------------------------
    % Defaults
    % ------------------------------------------------------------

    doprint = 0;        % print output
    doz = 1;            % convert to z-scores
    bcov = [];          % between-subjects covariate
    states = [];        % state IDs for each matrix on 3rd dim
    contrast = [];      % contrast across states
    dointeractive = 1;      % interactive questions/output
    n = [];             % number of subjects
    
    % default names: N1-x
    %names for networks (groups of regions) - or add your own as input
    for i = 1:max(c)
        names{i} = ['N' num2str(i)];
    end

    % variable input arguments
    % these should be pairs in the form 'field name', value
    % goes through all var args and enters field names as vars in the workspace
    for i = 1:2:length(varargin)
        N = varargin{i};    % get name of variable
        eval([N ' = varargin{i+1};']);
    end

    % if you enter empty names, prompts to name them
    % this really just works with mvroi tool
    if isempty(names) && exist('DATA', 'var')
        
        name_networks_plugin;
    end

    % fix names if too short | missing
    if isempty(names) || length(names) < max(c)
        for i = 1:max(c)
            names{i} = ['N' num2str(i)];
        end
    end

    % between-group names
    bnames = {};
    for i = 1:max(c)
        for j = i+1:max(c)
            maxlen1 = length(names{i});
            maxlen2 = length(names{j});
            bnames{end+1} = sprintf('%s-%s',names{i}(1:maxlen1),names{j}(1:maxlen2));
        end
    end


    % ------------------------------------------------------------
    % If states and contrast are entered, compute contrast from xc
    % xc MUST BE states nested w/i subjects
    % ------------------------------------------------------------
    if ~isempty(contrast)

        %if max(state) ~= len, error('Length of contrast must equal Number of states'), end
        % apply a contrast to the 3rd dimension of xc
        % replace data with contrast data
        % create states vector as well

        [xc,states] = contrast3d(xc,contrast);  % Conditions are nested within subjects

        % string of contrast values for printing and figures
        constr = ['[' sprintf(repmat('%2.0f ',1,length(contrast)),contrast) ']'];
        out.contrast = contrast;
        out.constr = constr;
    end


    % ------------------------------------------------------------
    % extract within- and between-class scores from xc
    % ------------------------------------------------------------

    wh_el = (ones(length(c),length(c)) - eye(length(c),length(c)));     % mask to include only non-self-correlations
    % multiply this by individual elements
    z_el = zeros(size(wh_el));

    nclust=max(c);

    for i = 1:nclust        % classes are columns, subjects are rows
        for j = 1:nclust
            wh1 = find(c == i); wh2 = find(c == j);
            whmat = z_el;       % which values to choose, in matrix form
            whmat(wh1,wh2) = 1;
            whmat = whmat .* wh_el;     % mask out self-correlations
            wh = find(whmat);   %wh contains the entries for within/between
            % for i == j, we have 2 correlations for each, but this is OK if we're
            % just averaging, which we are

            for s = 1:size(xc,3)
                xcs = xc(:,:,s);        % select subject
                xcs = xcs(wh);          % correls for this subject
                if isempty(xcs)         % one-region cluster
                    classcor(i,j,s) = NaN;
                else
                    classcor(i,j,s) = mean(xcs);  % matrix of correlations within and between classes
                end
            end

        end
    end

    % classcor is nclust*nclust*subs
    cmean = mean(classcor,3);       % mean within (diagonals) and between (off-diagonals)

    % ------------------------------------------------------------
    % Convert to z-scores, if requested, assuming correlation values are input
    % ------------------------------------------------------------

    if doz,
        tmp = classcor;
        zcor = .5 .* log( (1+tmp) ./ (1-tmp) );     % Fisher's r-to-z transform
    end

    % ------------------------------------------------------------
    % Print average correlation matrix and Statistics
    % ------------------------------------------------------------

    [cmeanz,t,sig,out.group_stats] = ttest3d(zcor);

    % for 2-D matrix (between-subjects corrs entered ONLY), get p-values
    if size(xc,3) == 1
        out.group_stats = [];
        cmean(isnan(cmean)) = 0;
        if isempty(bcov) && isempty(n)
            n = input('Enter number of subjects: ');
        elseif isempty(n), n = length(bcov);
        end

        % make vector form of correlations
        % this code does the same as the uncommented code below
        %rtmp = cmean-triu(cmean); rtmp(find(eye(size(rtmp)))) = 0;
        %rtmp = rtmp(:); rtmp(abs(rtmp)<eps) = [];

        if length(cmean(:)) == 1 % for one class only
            rtmp = cmean;
        else
            rtmp = cmean .* (1-eye(size(cmean)));
            rtmp = squareform(rtmp)';
        end

        [rci,sigu,z,p,rcrit] = r2z(rtmp,n,.05);
        out.group_stats.p = squareform(p);
        sigu = squareform(sigu);  % uncorrected
        [rci,sig,z,pc,rcrit] = r2z(rtmp,n,.05./length(rtmp)); sig = squareform(sig);  % corrected

        % diagonals
        [rci,sigu2,z,p,rcritu] = r2z(diag(cmean),n,.05);
        sigu = sigu + diag(sigu2);

        [rci,sig2,z,p,rcrit] = r2z(diag(cmean),n,.05./length(rtmp));
        sig = sig + diag(sig2);

        out.group_stats.sigu = sigu;
        out.group_stats.sig = sig;
    end


    if doprint
        if isempty(contrast), fprintf(1,['Avg. r within (diag) and btwn (off-diag) clusters']);,
        else fprintf(1,['Correlations among clusters, Contrast ' constr]);
        end
        print_correlation(cmean,sig,names);   % bonferroni corrected based on upper tri
    end

    % ------------------------------------------------------------
    % Within and between in 2-D matrix form
    % ------------------------------------------------------------

    wh_wi = find(eye(max(c),max(c)));
    wh_btwn = find(tril(ones(max(c),max(c)) - eye(max(c),max(c)))); %take only lower triangle to avoid repetitions

    bt = [];
    btz = [];
    
    for s = 1:size(xc,3)

        xcs = classcor(:,:,s);        %xcs is nclust by nclust
        wi(s,:) = xcs(wh_wi)';      % correls within

        if ~isempty(wh_btwn)
            bt(s,:) = xcs(wh_btwn)';    % correls between
        end

        xcs = zcor(:,:,s);
        wiz(s,:) = xcs(wh_wi)';      % fisher's z within (do stats on z-transformed)
        if ~isempty(wh_btwn)
            btz(s,:) = xcs(wh_btwn)';    % fisher's z between
        end

    end
    % code: 1 if sig uncorrected, 2 if sig corrected
    wis=sig(wh_wi)'+out.group_stats.sigu(wh_wi)';       % load values in from sig
    
    bts = [];
    if  ~isempty(wh_btwn)
        bts=sig(wh_btwn)'+out.group_stats.sigu(wh_btwn)';     % load values in from sig
    end
    
    % find NaN columns
    wh = find(all(isnan(wi),1));
    wi(:,wh) = 0;
    wiz(:,wh) = 0;

    out.class_corrs = classcor;
    out.mean_corrs = cmean;
    out.fisherz = zcor;
    out.within = wi;
    out.between = bt;
    out.withinz = wiz;
    out.betweenz = btz;
    out.class_names = names;
    out.between_names = bnames;
    out.between_covariate = bcov;



    % ------------------------------------------------------------
    % image the matrix and significance
    % ------------------------------------------------------------
    cmean(isnan(cmean)) = 0;

    % Check for existing figure, and create if necessary
    f1 = create_figure('class_corr_image');
    
% %     scnsize = get(0,'ScreenSize');
% %     myposition = [50 50 scnsize(3)-100 scnsize(4)/2];
% %     myposition(3) = min(myposition(3), 1200);
% %     myposition(4) = min(myposition(4), 500);
% % 
% %     figure('position',myposition,'color','white');
    
    fs=12;
    subplot(1,3,1)
    m = max(abs(cmean(:)));
    imagesc(cmean,[-m m]); colorbar; xlabel('Classes','FontSize',14),ylabel('Classes','FontSize',14),
    set(gca,'XTickLabel',names,'YTickLabel',names,'FontSize',fs,'XTick',1:length(cmean),'YTick',1:length(cmean));
    title('Correlations within (diagonals) and between (off-diagonals) clusters','FontSize',fs);
    if ~isempty(contrast), title(['Correlations among clusters, Contrast ' constr],'FontSize',fs); end

    subplot(1,3,2)
    imagesc(sig,[-.7 .7]), colormap gray, title('Sig. of avg correl (corrected)','FontSize',14)
    xlabel('Classes','FontSize',14),ylabel('Classes','FontSize',14),
    set(gca,'XTickLabel',names,'YTickLabel',names,'FontSize',fs,'XTick',1:length(cmean),'YTick',1:length(cmean));

    subplot(1,3,3)
    imagesc(out.group_stats.sigu,[-.7 .7]), colormap gray, title('Sig. of avg correl (uncorrected)','FontSize',14)
    xlabel('Classes','FontSize',14),ylabel('Classes','FontSize',14),
    set(gca,'XTickLabel',names,'YTickLabel',names,'FontSize',fs,'XTick',1:length(cmean),'YTick',1:length(cmean));

    % set colormap
    z = zeros(30,3)+1; b = z; b(:,[1 2]) = repmat((1:30)'./30,1,2); r = z; r(:,[2 3]) = repmat((30:-1:1)'./30,1,2);
    cm = [b; r];
    colormap(cm)

    drawnow

    % ------------------------------------------------------------
    % Univariate stats and descriptives on extracted data
    % plots
    % Univariate correlations with behavioral scores
    % ------------------------------------------------------------

    %disp('Correlations within networks')
    if length(size(xc))>2
        figure('position',[50 50 scnsize(3)-100 scnsize(4)/1.5],'color','white');
        subplot(2,1,1);
        barplot_columns(wi,'Correlations within networks',bcov,'nofig','noind');
        set(get(gca,'Title'),'FontSize',16)
        xlabel('Network','FontSize',14), ylabel('Mean correlation value','FontSize',14);

        %mark significant ones with an asterix
        for c=1:size(wi,2);
            yheight = abs(mean(wi(:,c)))+(ste(wi(:,c)));
            yheight = yheight + mean(abs(wi(:))) .* .10;
            yheight = yheight .* sign(mean(wi(:,c)));
            if yheight < 0, yheight = yheight - mean(abs(wi(:))) .* .3;,end

            xval = c - 0.1;

            if abs(wis(c))==1        % uncorrected
                text(xval,yheight,'*','color','k','fontsize',40);
            elseif abs(wis(c))==2;   % corrected
                text(xval,yheight,'**','color','k','fontsize',40);
            end

        end
        set(gca,'XTick',1:size(wi,2),'XTickLabel',names)

        %disp('Correlations between networks')
        subplot(2,1,2);
        barplot_columns(bt,'Correlations between networks',bcov,'nofig','noind');
        set(get(gca,'Title'),'FontSize',16)
        xlabel('Networks','FontSize',14), ylabel('Mean correlation value','FontSize',14)
        %if doz, ylabel('Correlation z-score'),else,ylabel('Correlation'),end
        set(gca,'XTick',1:size(bt,2),'XTickLabel',bnames)

        %mark significant ones with an asterix
        for c=1:size(bt,2);
            yheight = abs(mean(bt(:,c)))+(ste(bt(:,c)));
            yheight = yheight + mean(abs(bt(:))) .* .10;
            yheight = yheight .* sign(mean(bt(:,c)));
            if yheight < 0, yheight = yheight - mean(abs(bt(:))) .* .3;,end

            xval = c - 0.15;

            if abs(bts(c))==1        % uncorrected
                text(xval,yheight,'*','color','k','fontsize',40);
            elseif abs(bts(c))==2;   % corrected
                text(xval,yheight,'**','color','k','fontsize',40);
            end

        end

        % equalize_axes(ah,1)
        drawnow
    else
        % 2-D, print text
        fprintf(1,'%3.2f\t',wi);
        fprintf(1,'\n')
        disp('Correlations between networks (based on averaging individual region r values)')
        correlation_to_text(squareform(bt),rcritu,names);
    end

    % ------------------------------------------------------------
    % Correlations with behavioral covariate(s): stepwise regression
    % ------------------------------------------------------------
    if length(size(xc))>2   % only for 3-D (corrs within subjects)

        % only one covariate (y) right now.
        if ~isempty(bcov)

            fprintf(1,'\nCorrelations between behavior and within/between correlations\n')
            fprintf(1,'%s\t',names{:})
            fprintf(1,'%s\t',bnames{:})
            fprintf(1,'\n')
            tmp = corrcoef([bcov wi bt]);
            fprintf(1,'%3.2f\t',tmp(1,:));
            fprintf(1,'\n\n')

            warning off
            fprintf(1,'\n-----------------------------------------------------------------\n');
            fprintf(1,'Stepwise regression: Predictions of behavior with class average scores');
            fprintf(1,'\n-----------------------------------------------------------------\n');


            x = [wi bt];
            out.STEP = stepwise_tor(x,bcov,[names bnames]);
            %[b,se,pval,inmodel,stats] = stepwisefit(x,bcov,'penter',.10,'display','on');
            warning on


        end

    end


    return