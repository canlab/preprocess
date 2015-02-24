function [best_params,fit,beff,in] = tor_ga(gensize,numgen,inputs,ofun,varargin)
    %[best_params,fit,beff,in] = tor_ga(gensize,numgen,inputs,ofun,[fixed inputs])
    %
    % ----------------------------------------------------------------------
    % inputs    a cell array describing the inputs to the optimization function
    % (parameters to be optimized).
    % ----------------------------------------------------------------------
    % Each cell of inputs is a p x q matrix of parameters.
    % Internally, a set of 'organisms' is created that is
    % params x params x organisms (3-D).
    % This matrix is subject to crossover across orgs. separately for each
    % cell.
    % If inputs is a p x q matrix, as it will be placed in a single cell.
    %
    % There can be more than one set of
    % parameters that are combined in some way by ofun to produce a fitness
    % value.  if there is more than one set of input parameters,
    % inputs should be entered as a cell array, one cell per input.
    % inputs should be in ORDER of inputs entered to ofun!
    %
    % ----------------------------------------------------------------------
    % ofun      the objective function that combines the inputs.
    % ----------------------------------------------------------------------
    % There are two options for passing this in:
    % 1) enter the name of the function as a string.  the program creates a handle for the
    % function, and evaluates it using inputs specified in the inputs variable.
    % In this case, pass in fixed inputs after ofun, in the varargin fields
    % fixed inputs  optional, fixed inputs that do not change!  same structure
    % as inputs.
    %
    % 2) You can also enter ofun as a function handle durectly, with fixed inputs already
    % embedded before running the program.
    % The function should take as input a param list, and return fitness.
    % e.g., objhan = @(params) my_function_name(params,fixed_inputs1,fixed_inputs1,fixed_inputs1);
    % e.g., objhan = @(wh) prospect_organism(ceil(wh),pop,truep,iter);
    % Pass in objhan as the 'ofun' input argument
    %
    % by Tor Wager, Last updated: Jan 2007
    %
    % Examples:
    % ----------------------------------------------------------------------
    % Example for fitting indscal model:
    % inputs{1} = X; fixin{1} = sp; fixin{2} = B1; fixin{3} = B2;
    % tor_ga(30,10,inputs,'indscalf',fixin);
    %
    % W = rand(size(W)); W(1,:) = [10 10];, W(2,:) = [-10 -10];
    %inputs{2} = W;
    %
    % ---------------------------------------------------------------------
    % Example: Optimize gambles for prospect theory model
    % See prospect_optimize_design.m for definition of population of
    % gambles from which to draw (pop), truep, iter (all fixed inputs)
    %
    % objhan = @(wh) prospect_organism(ceil(wh),pop,truep,iter);
    % genfun = @() randsample(gindx,ntrials,'true')';
    % [best_params,fit,beff,in] = tor_ga(5,3,wh,objhan,genfun);
    % ---------------------------------------------------------------------

    t0 = clock;

    % --------------------------------------------------------------------
    % * set up inputs
    % --------------------------------------------------------------------

    % objective function
    % two modes: 1)  given string, construct string to evaluate with fixed inputs
    %            2)  given function handle with fixed inputs embedded, evaluate
    %            directly
    if isstr(ofun)
        evalmode = 'string';
        eval(['fun = @' ofun ';']), disp(['Objective function is ' ofun ])
    else
        evalmode = 'handle';
        fun = ofun;

    end

    % format inputs correctly
    if ~iscell(inputs), tmp=inputs; clear inputs; inputs{1} = tmp; end

    switch evalmode
        case 'string'
            % inputs - create string that tells feval what arguments to put in


            estr = 'f = feval(fun,in{1}(:,:,j)';
            for i = 2:length(inputs)
                estr = [estr ',in{' num2str(i) '}(:,:,j)'];
            end
        case 'handle'
    end

    % optional inputs - fixed, non-optimized inputs
    if length(varargin) > 0

        if strcmp(class(varargin{1}),'function_handle')
            % this is the org. generation function
            genfun = varargin{1};

        else
            % this is a fixed input

            fixin = varargin{1};
            if ~iscell(fixin), tmp=fixin; clear fixin; fixin{1} = tmp;,end
            for i = 1:length(fixin),
                estr = [estr ',fixin{' num2str(i) '}'];
            end

        end
    else
    end

    fprintf(1,'___________________________________________________________\n')
    switch evalmode
        case 'string'
            estr = [estr ');'];
            disp(['Evaluation string is ' estr])
        case 'handle'
            disp(['Will evaluate this objective function: '])
            disp(fun)
    end
    fprintf(1,'___________________________________________________________\n')


    % --------------------------------------------------------------------
    % * create start state
    % --------------------------------------------------------------------
    % create start state for fixed inputs

    for i = 1:length(inputs)

        if exist('genfun','var')
            % we have a custom organism-generation function handle
            for j = 1:gensize
                in{i}(:,:,j) = genfun();
            end
        else

            % determine range
            r = [min(inputs{i}(:)) max(inputs{i}(:))];

            % create in variable = cells are inputs, columns are param sets, rows
            % params within sets.
            % first one is always the input you put in!
            in{i}(:,:,1) = inputs{i}; % + randn(size(inputs{i})) .* std(inputs{i}(:));

            % rest of the population
            % 80% is random within 2*range of inputs
            for j = 2:round(gensize./.8)
                tmp = rand(size(inputs{i})); tmp = tmp.*r(2).*2 + r(1);
                in{i}(:,:,j) = tmp;
            end

            % 20% is input + noise
            for j = round(gensize./.8):gensize
                tmp = randn(size(inputs{i})) .* std(inputs{i}(:)) + inputs{i};
                in{i}(:,:,j) = tmp;
            end

        end


    end

    % --------------------------------------------------------------------
    % * iterate
    % --------------------------------------------------------------------

    f1 = tor_fig; set(gcf,'Color','w'); hold on; title('Max fitness'),xlabel('Generation')
    % % f2 = nmdsfig(inputs{1},ones(size(inputs{1},1),1),[]); set(gcf,'Position',[1955         701         560         420]);
    % % if length(inputs) > 1
    % %     f3 = nmdsfig(inputs{2},ones(size(inputs{2},1),1),[]); set(gcf,'Position',[2521         699         560         420]);
    % % end

    meff = [];
    beff = NaN .* zeros(1,numgen);
    gentime = NaN .* zeros(1,numgen);

    for i = 1:numgen

        t1 = clock;

        str = sprintf('Generation: %3.0f  ',i); fprintf(1,str);
        str = sprintf('Eval. fitness '); fprintf(1,str);

        for j = 1:gensize

            % --------------------------------------------------------------------
            % * make models
            % --------------------------------------------------------------------

            % --------------------------------------------------------------------
            % * test models
            % --------------------------------------------------------------------
            switch evalmode
                case 'string'
                    % if using function eval string
                    eval(estr)          % e.g., f = fun(in{1}(:,:,j),fixin{1},fixin{2},fixin{3});
                    fit(j,i) = f;       % fitness of each model in each generation
                case 'handle'
                    fit(j,i) = fun(in{1}(:,:,j));
            end

            % test program
            %estr2 = '[fi2,X,W,XL,XR] = feval(fun,in{1}(:,:,j),fixin{1},fixin{2},fixin{3});'
            %eval(estr2)
            %sp = fixin{1};
            %for ii = 1:size(W,1), Brecon(:,:,ii) = XL * diag(W(ii,:)) * XR';, end
            %fit2 = (sp - Brecon).^2; fit2 = sqrt(sum(fit2(:))./prod(size(Brecon)));
            %fit2 = 1 ./ (1 + fit2);
            %fprintf(1,'RMS logit in orig ga = %3.2f, in f2 = %3.2f, recalc %3.2f\n',fit(j,i),fi2,fit2)

        end

        erase_string(str);
        str = sprintf('Crossover '); fprintf(1,str);

        % --------------------------------------------------------------------
        % * save best and crossover
        % --------------------------------------------------------------------
        eff = fit(:,i)';

        % % %if i == numgen, break, end
        
        if ~(any(eff > median(eff))) | (length(meff) > 100 & length(unique(meff(end - min(100,length(meff)-1):end))) == 1), disp(['System converged at generation ' num2str(i)]), break, end
        % two criteria: no efficiency above median efficiency, and no change in the last 100 generations

        % run crossover for each separate input matrix to be recombined
        % best_params contains best parameters of this generation
        for j = 1:length(inputs)
            % replace in with crossover; b = index of best before xover
            % best_params is best of this generation
            % index of best after xover is 1 (hard-coded)
            [in{j},b,best_params{j}] = xover(in{j},eff);
        end

        %     % run crossover for each separate input matrix to be recombined
        %     % best_params contains best parameters of this generation
        %     hh = findobj('Color','b'); delete(hh);
        %     for j = 1:length(inputs),
        %         [in2{j},b,best_params{j}] = xover(in{j},eff);
        %
        % % % %         tmp = best_params{j} - inputs{j};
        %         %if any(tmp(:)) > eps
        %             %fprintf('Updated input %3.0f\t',j),
        %
        % % % %             % update graphics figure; commented out now (specific to
        % % % %             % indscal)
        % % % %             if j == 1, figure(f2);, elseif j==2, figure(f3);,end
        % % % %             hold on; %plot(best_params{j}(:,1),best_params{j}(:,2),'ko','MarkerFaceColor','r');
        % % % %             wh = find(any(tmp,2));
        %
        %             %plot(best_params{j}(wh,1),best_params{j}(wh,2),'bo','MarkerFaceColor','b');
        % % % %             plot(best_params{j}(:,1),best_params{j}(:,2),'bo','MarkerFaceColor','b');
        %             %end
        %
        %     end
        %fprintf(1,'\n')


        erase_string(str);
        gentime(i) = etime(clock,t1);

        str = sprintf('Time: %3.0f ',etime(clock,t1)); fprintf(1,str);

        % shouldn't need this stuff...
        % %         % print fitness of best of generation
        % %         j = 1;  % evaluates the best one, which is saved by xover in the first slot
        % %
        % %         switch evalmode
        % %                 case 'string'
        % %                     % if using function eval string; f is output
        % %                     eval(estr)          % e.g., f = fun(in{1}(:,:,j),fixin{1},fixin{2},fixin{3});
        % %
        % %                 case 'handle'
        % %                     f = fun(in{1}(:,:,j));
        % %         end

        % use original eff vector to get best, in case fitness is stochastic
        meff = eff(eff == max(eff)); meff = meff(1);
        fprintf(1,'Best: # %3.0f, Fitness: %3.2f\n',b,meff)

        beff(i) = meff;
        figure(f1); plot(beff,'k','Linewidth',2); drawnow

    end


    % final report
    fprintf(1,'\nGA Finished\n')
    fprintf(1,'___________________________________________________________\n')
    fprintf(1,'Generations: %3.0f \nOrganisms per generation: %3.0f\n',numgen,gensize);
    fprintf(1,'\nInitial fitness: %3.4f \nFinal fitness: %3.4f\n',beff(1),beff(end));

    fprintf(1,'\nAverage time per iteration: %3.0f\n',mean(gentime));
    fprintf(1,'Average time per organism: %3.0f\n',mean(gentime)./gensize);
    fprintf(1,'Total time: %3.0f\n',sum(gentime));
    fprintf(1,'___________________________________________________________\n')

    return




    % --------------------------------------------------------------------
    % * get best half, crossover to fill in lists
    % --------------------------------------------------------------------
function [newvec,b,best_params] = xover(paramvec,eff)
    %
    % paramvec is 3-D matrix, 3rd dim is realization (organism)
    % columns are sets of params to be crossed over.
    %
    % b is index of best param vec - may have random noise, etc. added in
    % newvec
    % best_params is original best param vec, stored in newvec(:,:,1)

    % best one before xover
    b = find(eff == max(eff)); b = b(1); best_params = paramvec(:,:,b);


    % add a little random noise to efficiency - 1% of std of efficiency
    eff = eff + randn(1,length(eff)) .* .01 * std(eff);

    w = find(eff > median(eff));

    % we can only do crossover if not all the designs are the same
    % ---------------------------------------------------
    if isempty(w)
        warning('Extremely homogenous sample!')
        % add a little random noise to efficiency - 1% of var of efficiency
        fitness = fitness + randn(1,length(fitness)) .* .01 * mean(fitness);
        w = 1:round(size(fitness,2)./2);

    end

    % save best half
    % ---------------------------------------------------
    newvec = paramvec(:,:,w);
    [nparams,dummy,last] = size(newvec);  % save size for adding random variation to existing


    %  %  CODE commented out for adding random noise right now.
    % % n = length(newvec(:));
    % % n = round(n/10);       % 10% - to add random noise

    % number of crossover points
    nxover = max(1,round(nparams ./ 50)); % 5;

    % fill in 2nd half with crossovers of first half
    % ---------------------------------------------------
    n_to_fill = size(paramvec,3) - last;

    for v = 1:n_to_fill

        % choose two random integers within 1st half
        w = ceil(rand(1,2) * last);

        babyv = [];

        for i = 1:size(newvec,2)        % for each set of parameters (columns)

            %babyv(:,i) = rcomb(newvec(:,i,w(1)),newvec(:,i,w(2)));

            babyv(:,i) = rcomb_multi(newvec(:,i,w(1)),newvec(:,i,w(2)),nparams,nxover);
        end

        % add to others
        newvec(:,:,last+v) = babyv;


    end

    % add random variation to existing best half
    % % % wh = round(rand(1,n) .* n);
    % % % wh(wh < 1) = 1; wh(wh > n) = n;
    % % %
    % % % nv = newvec(:,:,1:last);
    % % % nv(wh) = nv(wh) + randn(size(wh));
    % % % newvec = cat(3,nv,newvec(:,:,last+1:end));

    newvec(:,:,1) = best_params;      % re-insert best one


    return


    % --------------------------------------------------------------------
    % * crossover
    % --------------------------------------------------------------------
function c = rcomb(a,b)
    % combines 2 vectors of numbers at a random crossover point

    w = ceil(rand * (length(a) - 1));
    c = a(1:w);
    c = [c; b(w+1:end)];

    return

function  c = rcomb_multi(a,b,n,nxover)

% % n = length(a);
% % nxover = max(1,round(n ./ 50)); % 5;     % number of crossover points
st = randperm(n); 

% divide the vectors up into chunks at random points
st = [1 sort(st(1:nxover)) n]; 
en = [(st(2:end))-1]; st = st(1:end-1); en(end) = n;

c = a;
for i = 1:2:length(st)
    c(st(i):en(i)) = b(st(i):en(i));
end

% example:
% a = (1:100)';
% b = 1000+a;
% nxover = 10;
% ...run code...
% figure;plot(c)

return


function erase_string(str1)
    fprintf(1,repmat('\b',1,length(str1))); % erase string
    return



