function EXPT = make_htw_contrast_images(EXPT)
    % EXPT = make_htw_contrast_images(EXPT)
    %
    % This prompts for a set of contrasts across height/delay/wid images, and
    % calculates them, saving in the same dir as original images.
    %
    % Uses images stored in EXPT.NLCON  --> from get_htw_image_names.m
    % Saves names stored in EXPT.SNPM   --> for random effects
    %
    % Use with whole brain filter / whole_brain_fir
    % or SPM2 FIR basis set after creating height, time to peak, and width
    % images


    % number of conditions
    num_contrasts = length(EXPT.NLCON.height);
    num_subjs = size(EXPT.NLCON.height{1},1);

    if num_subjs ~= length(EXPT.subjects)
        warning('NO OF SUBJECTS DOES NOT MATCH NO OF NLCON IMAGES: ignoring for now, but behavioral regs, etc. could be off!');
    end

    % define contrasts

    go = 1;

    if isfield(EXPT.NLCON,'contrasts') && isfield(EXPT.NLCON,'connames')
        use_existing_contrasts = input('Existing NLCON.contrasts/connames already exists.  Use it (1/0)? ');
        if use_existing_contrasts
            go = 0;
        else
            EXPT.NLCON.contrasts = [];
            EXPT.NLCON.connames = {};
        end
    else
        EXPT.NLCON.contrasts = [];
        EXPT.NLCON.connames = {};
    end


    while go
        new_contrast = input(['Enter contrast across ' num2str(num_contrasts) ' event types in [], return to quit: ']);

        if isempty(new_contrast)
            go = 0;
        else
            EXPT.NLCON.contrasts(end+1,:) = new_contrast;
            EXPT.NLCON.connames{end+1} = input('Enter short name for this contrast, no spaces or special chars: ','s');
        end
    end


    for i = 1:num_subjs
        for contrast = 1:size(EXPT.NLCON.contrasts,1)  % contrasts are rows
            [h,d,w] = get_contrast(EXPT,i,contrast,num_contrasts);
            disp(sprintf('Done: %s %s %s',h,d,w));

            % save output in SNPM for rfx
            if i == 1
                EXPT.SNPM.heightP{contrast} = h;
                EXPT.SNPM.delayP{contrast} = d;
                EXPT.SNPM.widthP{contrast} = w;
            else
                EXPT.SNPM.heightP{contrast} = str2mat(EXPT.SNPM.heightP{contrast},h);
                EXPT.SNPM.delayP{contrast} = str2mat(EXPT.SNPM.delayP{contrast},d);
                EXPT.SNPM.widthP{contrast} = str2mat(EXPT.SNPM.widthP{contrast},w);
            end
        end
    end

    % extra names for convenience in EXPT.SNPM

    EXPT.SNPM.P = [EXPT.SNPM.heightP EXPT.SNPM.delayP EXPT.SNPM.widthP];
    EXPT.SNPM.connums = 1:length(EXPT.SNPM.P);
    EXPT.SNPM.connames = [];
    for i = 1:length(EXPT.NLCON.connames)
        EXPT.SNPM.connames = str2mat(EXPT.SNPM.connames,[EXPT.NLCON.connames{i} '_Height']);
    end
    for i = 1:length(EXPT.NLCON.connames)
        EXPT.SNPM.connames = str2mat(EXPT.SNPM.connames,[EXPT.NLCON.connames{i} '_Delay']);
    end
    for i = 1:length(EXPT.NLCON.connames)
        EXPT.SNPM.connames = str2mat(EXPT.SNPM.connames,[EXPT.NLCON.connames{i} '_Width']);
    end
    EXPT.SNPM.connames(1,:) = [];

end




function [h,d,w] = get_contrast(EXPT,i,contrast,num_contrasts)
    % i is subject
    % contrast is contrast number
    % myc is contrast values

    myc = EXPT.NLCON.contrasts(contrast,:);

    % get list of images
    P = EXPT.NLCON.height{1}(i,:);
    for j = 2:num_contrasts
        P = str2mat(P,EXPT.NLCON.height{j}(i,:));       % should have path
    end

    Q = [EXPT.NLCON.connames{contrast} '_height.img'];  % no path

    Vo = contrast_image(P,Q,myc);
    h = Vo.fname;



    myc = EXPT.NLCON.contrasts(contrast,:);

    % get list of images
    P = EXPT.NLCON.delay{1}(i,:);
    for j = 2:num_contrasts
        P = str2mat(P,EXPT.NLCON.delay{j}(i,:));       % should have path
    end

    Q = [EXPT.NLCON.connames{contrast} '_delay.img'];  % no path

    Vo = contrast_image(P,Q,myc);
    d = Vo.fname;



    myc = EXPT.NLCON.contrasts(contrast,:);

    % get list of images
    P = EXPT.NLCON.width{1}(i,:);
    for j = 2:num_contrasts
        P = str2mat(P,EXPT.NLCON.width{j}(i,:));       % should have path
    end

    Q = [EXPT.NLCON.connames{contrast} '_width.img'];  % no path

    Vo = contrast_image(P,Q,myc);
    w = Vo.fname;
end