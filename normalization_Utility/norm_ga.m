function [prm, normdat] = norm_ga(prm, norg, numgen, varargin)
    % [prm, normdat] = norm_ga(prm, norg, numgen, ['coarse' 'fine'], ['name', imagename], ['template', templatename], ['mask', maskInfo])
    %
    % Use a genetic algorithm to find nonlinear warping parameters that
    % maximize the monotonic relationship between a template or target image and an
    % object image to normalize.
    %
    % prm is a parameter set from spm_normalize
    % This is stored in the .mat file associated with the normalization
    % Start with an image already normalized by SPM.
    %
    % Optional inputs:
    % 'coarse': parameters sampled from starting distribution
    % 'fine' : parameters that vary between -.5 and .5 added to starting estimates
    % 'mask' : followed by maskInfo volume info structure with list of voxels
    %           in mask (see iimg_read_img, extended output)
    % 'nbf'  : followed by [x y z] num. of basis functions
    %
    % example:
    % prm = load('T1_sn.mat');
    % prm = norm_ga(prm, 70, 20);
    % Vout = prm.VG;
    % Vout.fname = 'normalized_t1.img';
    % spm_write_vol(Vout, normdat);
    % spm_image('init', Vout.fname);
    %
    % tor wager, v1.2, nov. 06
    %
    %
    % IMPROVEMENTS/ISSUES:
    % 1) Add multi-resolution, coarse search then finer search
    % 2) Work on epi images directly?
    % 3) speed speed speed
    % 4) Simulated annealing?
    % 5) Better control over random perturbation?
    % 6) flexible input/output names, path handing (needed?)
    %
    % Example: re-do a subject's GA with lower # basis fcns
    %------------------------------------------------------
    % template = 'mean.img';
    % obj = 'mean1011.img';
    % matname = 'norm_subj010_sn3d_new.mat';
    %
    % flags = [];
    % flags.smosrc = 8;
    % flags.cutoff = 50;
    % flags.reg = .2;
    %
    % VG = spm_vol(template);
    % VF = spm_vol(obj);
    % prm = spm_normalise(VG, VF, matname, [], [], flags);
    %
    % ndat = norm_apply_params(prm, 'data', spm_read_vols(V), 'mat', V.mat);
    % outname = 'w_test_subj0011.img';
    % V.fname = outname;
    % spm_write_vol(V, ndat);
    % spm_image('init', outname);

    % maskInfo = iimg_read_img('mask.img', 2);
    % %%don't do this if using SPM norm above. % prm.Tr = prm.Tr(1:nbf(1), 1:nbf(2), 1:nbf(3),:);
    % prm = norm_ga(prm, 5, 3, 'nbf', nbf, 'name', outname, 'template', 'mean.img', 'mask', maskInfo);
    % %%if using SPM norm above. %
    % prm = norm_ga(prm, 5, 3, 'name', outname, 'template', 'mean.img', 'mask', maskInfo);
    % V = spm_vol(outname);
    % ndat = norm_apply_params(prm, 'data', spm_read_vols(V), 'mat', V.mat);
    % V.fname = outname;
    % spm_write_vol(V, ndat);
    % spm_image('init', outname);
    %
    % Evaluate fitness of outname:
    % [imagedata, volInfo] = iimg_mask('mask.img', outname);
    % [meandata, volInfo] = iimg_mask('mask.img', 'mean.img');
    % fitness = norm_eval_fitness(meandata(volInfo.wh_inmask), imagedata(volInfo.wh_inmask), 1)

    % -----------------------------------------------------
    % Process input arguments
    % -----------------------------------------------------
    popType = 'fine';       % type of population parameters to start with
    dosavegeneration = 0;   % save mat file for each generation in case of crash
    nbf = [];               % # basis fcns: [] = use defaults based on existing prm

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'fine', popType = 'fine';
                case 'coarse', popType = 'coarse';
                case 'write', outname = varargin{i+1}; varargin{i+1} = [];

                case 'name', name = varargin{i+1}; varargin{i+1} = [];
                case 'template', template = varargin{i+1}; varargin{i+1} = [];

                case 'mask', maskInfo = varargin{i+1};

                case 'nbf', nbf = varargin{i+1};

                otherwise
                    error('Unknown option');
            end
        end
    end



    %  load data
    % --------------------------------------------
    fprintf('Loading data.\n');

    if ~exist('name', 'var'), name = 'T1.img'; end
    if ~exist('template', 'var'), template = which('avg152T1.nii'); end

    fprintf('Image to normalize: %s\n Template: %s\n', name, template);

    % %     name = which(name);
    % %     if isempty(name), error('Cannot find image.'); end
    % %
    % %     template = which(template);
    % %     if isempty(template), error('Cannot find template.'); end

    V = spm_vol(name);
    dat = spm_read_vols(V);
    tempdat = spm_read_vols(spm_vol(template));

    tempdat(isnan(tempdat)) = 0;

    %  set sizes, etc.
    % --------------------------------------------
    Trvec = prm.Tr(:);      % vectorize so we can test
    npar = size(Trvec, 1);   % number of parameters in optimization

    [i, j, m, n] = size(prm.Tr);
    sz = [i j m n];

    startpop = zeros(npar, norg);

    % build population
    % --------------------------------------------
    startpop = build_population;


    % initialize GA output and vars
    % --------------------------------------------
    fit_by_gen = zeros(1, numgen);

    f1 = tor_fig(1, 2); set(gcf, 'Color', 'w');
    hold on; title('Max fitness'), xlabel('Generation')
    subplot(1, 2, 2);
    title('Fitness of all vectors');
    ylabel('Fitness'); xlabel('Organism')

    % Create basis set for use in all iterations
    % --------------------------------------------
    [z, X, Y, BX, BY, BZ, def_flags, d] = norm_create_basis(prm, nbf);


    fprintf('Setup done: Starting iterations\n');
    if dosavegeneration
        fprintf('Saving each generation in case of crash:  tmp_norm_ga_bestsolution.mat\n');
    end

    % --------------------------------------------------------------------
    % * iterate over generations
    % --------------------------------------------------------------------
    for i = 1:numgen

        t1 = clock;

        % test population fitness
        % --------------------------------------------
        fitness = test_population;

        % save best and crossover
        % --------------------------------------------

        [startpop, best_vec, max_fitness] = xover(startpop, fitness);

        fit_by_gen(i) = max_fitness;

        % best fitness of each generation to date
        myfit = fit_by_gen(fit_by_gen ~= 0);

        % plot
        % --------------------------------------------
        plot_fitness_lines();

        if exist('genstr', 'var')
            erase_string(genstr);
        end

        genstr = sprintf('Generation %3.0f: Best fitness = %3.4f, elapsed = %3.0f s\n', i, max_fitness, etime(clock, t1));
        fprintf(genstr);

        if dosavegeneration
            save tmp_norm_ga_bestsolution best_vec j max_fitness myfit
        end

        % converge
        % --------------------------------------------
        % two criteria for stopping: no efficiency above median efficiency, 
        % and no change in the last 100 generations
        if i == numgen, break, end
        if ~(any(fitness > nanmedian(fitness))) || ...
                (length(myfit) > 100 && ...
                length(unique(myfit(end - min(100, length(myfit)-1):end))) == 1)

            disp(['System converged at generation ' num2str(i)])
            break
        end

    end

    prm.oldTr = prm.Tr;
    prm.Tr = reshape(best_vec, sz);
    prm.GA.generations = numgen;
    prm.GA.norg = norg;


    % apply parameters
    normdat = norm_apply_params(prm, 'data', dat);

    % write output image
    if exist('outname', 'var')
        [dd, ff, ee] = fileparts(name);
        outname = ['wGA_' ff ee];
        fprintf('Writing image: %s\n', outname);

        Vout = prm.VG;
        Vout.fname = outname;
        spm_write_vol(Vout, normdat);
    end



    % --------------------------------------------------------------------
    %
    %%% End of main function %%%
    %
    % --------------------------------------------------------------------


    %%% Nested functions %%%
    function startpop = build_population
        startpop(:,1) = Trvec;  % original norm result
        for k = 2:norg

            switch popType
                case 'coarse'
                    % parameters sampled from starting distribution
                    startpop(:,k) = randsample(Trvec, npar, true); % with replacement
                case 'fine'
                    % parameters that vary between -.5 and .5 * multfact added to
                    % starting estimates
                    multfact = 4;
                    startpop(:,k) = Trvec +  multfact .* rand(length(Trvec), 1) - .5;
                    %fprintf('Perturbations range from %3.2f to %3.2f'
                otherwise
                    error('Unknown popType. Check function help.');
            end
        end
    end




    function fitness = test_population
        % test population fitness
        % --------------------------------------------
        fitness = zeros(1, norg);

        domask = 0;
        if exist('maskInfo', 'var') && ~isempty(maskInfo)
            domask = 1;
            wh = maskInfo.image_indx;
        end

        for myorg = 1:norg

            fprintf('%04d', myorg);

            mytr = reshape(startpop(:,myorg), sz);
            prm.Tr = mytr;
            %normdat = spm_just_applynorm(prm, 'data', dat);
            normdat = norm_apply_params(prm, 'data', dat, 'mat', V.mat, 'basis', z, X, Y, BX, BY, BZ, def_flags, d);
            normdat(isnan(normdat)) = 0;

            if domask
                fitness(myorg) = norm_eval_fitness(tempdat(wh), normdat(wh), 0);
            else
                fitness(myorg) = norm_eval_fitness(tempdat(:), normdat(:), 0);
            end
            %[fitness, tempvec, normvec, yhat_sorted] = norm_eval_fitness(tempdat, normdat, 1);

            fprintf('\b\b\b\b');
        end
    end


    function plot_fitness_lines
        figure(f1); subplot(1, 2, 1);

        if exist('myline', 'var') && ishandle(myline), delete(myline); end
        myline = plot(1:length(myfit), myfit, 'bo-', 'MarkerFaceColor', 'b');

        subplot(1, 2, 2);
        if exist('myline2', 'var') && ishandle(myline2), delete(myline2); end
        myline2 = plot(sort(fitness), 'b-');
        if i == 1
            plot(sort(fitness), 'r-');
        end
        drawnow
    end
end  % main function








% --------------------------------------------------------------------
% * get best half, crossover to fill in lists
% --------------------------------------------------------------------

function [newvec, best_vec, max_fitness] = xover(paramvec, fitness)
    %
    % paramvec is 2-D matrix, 2nd dim is realization (organism)
    % columns are sets of params to be crossed over.
    %
    % best_indx is index of best param vec - may have random noise, etc. added in
    % newvec
    % best_vec is original best param vec, stored in newvec(:,1)

    % save best one - now done outside xover
    best_indx = find(fitness == max(fitness)); best_indx = best_indx(1);

    best_vec = paramvec(:,best_indx);
    max_fitness = fitness(best_indx);

    % add a little random noise to efficiency - 1% of std of efficiency
    fitness = fitness + randn(1, length(fitness)) .* .01 * std(fitness);

    w = find(fitness > median(fitness));

    % we can only do crossover if not all the designs are the same
    % ---------------------------------------------------
    if isempty(w)
        warning('Extremely homogenous sample!')
        % add a little random noise to efficiency - 1% of var of efficiency
        fitness = fitness + randn(1, length(fitness)) .* .01 * mean(fitness);
        w = 1:round(size(fitness, 2)./2);

    end

    % save best half
    % ---------------------------------------------------
    newvec = paramvec(:,w);
    [nparams, last] = size(newvec);  % save size for adding random variation to existing

    % number of crossover points
    nxover = max(1, round(nparams ./ 50)); % 5;


    % fill in 2nd half with crossovers of first half
    % ---------------------------------------------------
    n_to_fill = size(paramvec, 2) - size(newvec, 2);

    for i = 1:n_to_fill

        % choose two random integers within 1st half
        w = ceil(rand(1, 2) * last);

        % combine into new vector
        %babyv = rcomb(newvec(:,w(1)), newvec(:,w(2)));
        babyv = rcomb_multi(newvec(:,w(1)), newvec(:,w(2)), nparams, nxover);

        % add to others
        newvec(:,last+i) = babyv;

    end

    % add random variation to existing best half and noise to all (?)
    % ---------------------------------------------------
    n = length(newvec(:));
    n = round(n/10);       % 10% - to add random noise
    wh = round(rand(1, n) .* n);
    wh(wh < 1) = 1; wh(wh > n) = n;

    nv = newvec(:,1:last);
    nv(wh) = nv(wh) + randn(size(wh, 1), 1);


    newvec = cat(2, nv, newvec(:,last+1:end));

    % save best in first slot, without noise
    % ---------------------------------------------------
    newvec(:,1) = best_vec;      % re-insert best one
end


% --------------------------------------------------------------------
% * crossover
% --------------------------------------------------------------------
function c = rcomb(a, b)
    % combines 2 vectors of numbers at a random crossover point

    w = ceil(rand * (length(a) - 1));
    c = a(1:w);
    c = [c; b(w+1:end)];
end



function  c = rcomb_multi(a, b, n, nxover)
    % % n = length(a);
    % % nxover = max(1, round(n ./ 50)); % 5;     % number of crossover points
    st = randperm(n);

    % divide the vectors up into chunks at random points
    st = [1 sort(st(1:nxover)) n];
    en = [(st(2:end))-1]; 
    st = st(1:end-1); 
    en(end) = n;

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
end



function erase_string(str1)
    fprintf(repmat('\b', 1, length(str1))); % erase string
end
