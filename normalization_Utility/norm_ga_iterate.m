% norm_ga_iterate(imagenames, prmnames, norg, numgen, iter)
%
% Iterative GA normalization to the mean image in a set
% Start with UNnormalized images and normalization parameters from SPM
% Then run this on them to improve warping.
%
% imagenames = filenames('*/T1.img', 'absolute', 'char');
% prmnames = filenames('*/*_sn.mat', 'absolute', 'char');
% norg = 10;
% numgen = 10;
% iter = 10;
% norm_ga_iterate(imagenames, prmnames, norg, numgen, iter)
%
% Mask is interpolated to the space defined by the first image
% and saved as mask.img
%
% norm_ga_iterate(imagenames, prmnames, norg, numgen, iter, 'mask', maskname)
%
% Tor Wager
% Nov. 06, version 11/26
%
% Create a mask image from previous mean:
% img = 'mean.img'
% V = spm_vol(img);
% v = spm_read_vols(V);
% v = v > 10;
% VO = V; VO.fname = 'startmask.img';
% spm_write_vol(VO, v);
%
% Example: re-do a subject's GA with lower # basis fcns
% maskInfo = iimg_read_img('mask.img', 2);
% prm = norm_ga(prm, 5, 3, 'nbf', [3 4 3], 'name', 'w010_mean1011.img', 'template', 'mean.img', 'mask', maskInfo);
function norm_ga_iterate(imagenames, prmnames, norg, numgen, iter, varargin)
    diary Norm_Ga_Log.txt

    imagenames = char(imagenames);
    prmnames = char(prmnames);
    n = size(imagenames, 1);

    % -----------------------------------------------------
    % SETUP mat files and mask image
    % -----------------------------------------------------
    % read param files, make sure we can load them, 
    % save copies in norm_subj???_sn3d.mat
    % replace prmnames with these names

    setup_param_filenames;
    setup_norm;

    % Create initial mean
    [meandata, fitness, prm] = mean_warped_image(imagenames, prmnames, 'mask', 'mask.img');
    !mv mean.img mean_original.img
    !mv mean.hdr mean_original.hdr
    spm_check_registration(char('mean_original.img', 'mean_original.img', 'mask.img')); 
    spm_orthviews_name_axis('Original', 1); 
    spm_orthviews_name_axis('Latest', 2); 
    spm_orthviews_name_axis('Mask', 3); 

    % -----------------------------------------------------
    % Iterate steps of mean image creation and warping
    % -----------------------------------------------------
    fitness_by_iteration = zeros(1, n);

    for iteration = 1:iter

        % Create mean image
        % -----------------------------------------------------
        % write mean.img, which is the new template
        %
        % also load param files for all subjects and return warping info for all
        % subjects in cells of prm.

        [meandata, fitness, prm] = mean_warped_image(imagenames, prmnames, 'mask', 'mask.img');

        fitness_by_iteration(iteration,:) = fitness;


        banner;

        % Warp each image to mean
        % -----------------------------------------------------
        % run Genetic Algorithm on each image, with mean.img as template

        for i = 1:n

            % setup output name
            [VG, save_subj_str] = get_savename;

            % get a new set of warps for one subject using GA
            myperm = norm_ga(prm{i}, norg, numgen, 'name', deblank(imagenames(i,:)), 'template', 'mean.img', ...
                'mask', maskInfo, popType);

            % save warping parameters for this subjects in the sn3d .mat file
            Tr = myperm.Tr;
            disp(save_subj_str), disp(' ');
            eval(save_subj_str)

            % save figure
            saveas(gcf, sprintf('norm_ga_subj%03d_iter%02d', i, iter), 'png');
            close  
    
        end

        save Norm_Ga_Results fitness_by_iteration prm

        show_averages();

    end  % iterations


    % write final mean image and individual normalized images

    [meandata, fitness, prm] = mean_warped_image(imagenames, prmnames, 'write');


    
    fitness_by_iteration(end+1,:) = fitness;

    diary off

    %%% end main function %%%






    %%% Inline functions %%%

    % ___________________________________________________________________
    %
    % Inline: Print banner for start of one iteration
    % ___________________________________________________________________

    function banner

        fprintf(1, '_____________________________________________________\n\n')
        fprintf(1, 'Iteration : %3.0f\n', iteration)

        fprintf(1, 'Fitness by iteration: mean\n')
        fprintf(1, '%3.4f ', mean(fitness_by_iteration, 2)');
        fprintf(1, '\n')

        fprintf(1, 'Fitness by iteration: min\n')
        fprintf(1, '%3.4f ', min(fitness_by_iteration, [], 2)');
        fprintf(1, '\n')

        fprintf(1, '\n_____________________________________________________\n')

    end

    % ___________________________________________________________________
    %
    % Inline: Setup param file names
    % ___________________________________________________________________

    function setup_param_filenames
        % read norm param files and save in ga_sn3d.mat files
        for i = 1:n

            % read norm params
            prm{i} = load(deblank(prmnames(i,:)));

        end

        % read, make struct, and save
        newprmnames = [];
        docopy = 1;

        for i = 1:n

            % read norm params
            prm{i} = load(deblank(prmnames(i,:)));

            % save params in specifically named files
            outname = sprintf('norm_subj%03d_sn3d.mat', i);

            if i == 1
                newprmnames = outname;

                if exist(outname, 'file')
                    docopy = input('Warning! Overwrite old sn3d.mat files? (1/0) ');
                    %if ~(cont == 1), error('Quitting.'), end
                end

            else
                newprmnames = str2mat(newprmnames, outname);
            end

            if docopy, copyfile(deblank(prmnames(i,:)), outname); end
        end

        prmnames = newprmnames;
        if docopy, disp('Created:');  else, disp('Using existing:'); end
        disp(prmnames);
    end


    % ___________________________________________________________________
    %
    % Inline: Setup mask, population type, template for Normalization
    % ___________________________________________________________________

    function setup_norm

        popType = 'fine';
        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}
                    case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
                    case {'coarse', 'fine'}, popType = varargin{i};

                    otherwise
                        error('Unknown string option');
                end
            end
        end

        % template image
        [dd, ff, ee] = fileparts(prm{1}.VG(1).fname);
        template = which([ff ee]);
        fprintf(1, 'Found template image: \n%s\n\n', template);
        fprintf(1, 'This will be used to define the space of the mask image.\n');
        
        % mask image
        % make sure this is in the space defined by the first image
        disp('Creating mask.img in same space as template');
        if exist('mask', 'var')
            fprintf(1, 'Found input image: %s\n', mask);
            maskdat = iimg_reslice(template, mask, 'write', 'outname', 'mask.img');
        else
            fprintf(1, 'Mask is all ones\n');
            V = spm_vol(template);
            maskdat = ones(V.dim(1:3));
            VO = V;
            VO.fname = 'mask.img';
            spm_write_vol(VO, maskdat);
        end

        maskInfo = iimg_read_img('mask.img', 2); % get info struct w/list of in-mask vox
        
        save norm_ga_SETUP prm template imagenames
    end

    % ___________________________________________________________________
    %
    % Inline: Get name of output param file and save eval string for 1
    % subject
    % ___________________________________________________________________
    function [VG, save_subj_str] = get_savename
        VG = spm_vol('mean.img');
        [dd, ff, ee] = fileparts(deblank(prmnames(i,:)));
        save_subj_str = ['save ' ff ' -append VG Tr'];
        disp('Will save params for this subject to:');
        disp(save_subj_str)
    end
    
end


function show_averages
    spm_check_registration(char('mean_original.img', 'mean.img', 'mean.img'));
    [volInfo, dat] = iimg_read_img(char('mean_original.img','mean.img'), 2);
    datd = diff(dat')';

    dattmp = datd;
    dattmp(abs(dattmp) < eps*100) = [];

    % % lb = prctile(dattmp, 10);
    % % ub = prctile(dattmp, 90);
    % % datd(datd < ub & datd > lb) = 0;


    bb = prctile(abs(dattmp), 99);
    thr = prctile(abs(dattmp), 90);
    datd(datd > -thr & datd < thr) = 0;

    cl = iimg_indx2clusters(datd, volInfo);
    cluster_orthviews(cl, 'add', 'handle', 3);

    cm = spm_orthviews_hotcool_colormap([-bb bb], thr/4);

    dattmp = datd;
    dattmp(abs(dattmp) < eps*100) = [];
    create_figure('Hist'); hist(dattmp, 100);
    plot_vertical_line([-thr thr]);
    spm_orthviews_name_axis('Original', 1);
    spm_orthviews_name_axis('Latest', 2);
    spm_orthviews_name_axis('Difference', 3);

end
