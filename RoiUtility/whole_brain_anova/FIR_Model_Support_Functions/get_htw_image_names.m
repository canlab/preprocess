% EXPT = get_htw_image_names(EXPT,varargin)
% var arg in: get first n images only
%
% tor wager
% start in directory above individual results
% extracts the following images:
%   height*img     height of estimated HRF
%   delay*img      time to peak of est. HRF
%   width*img      FWHM est. of HRF
%
% stores images in EXPT.NLCON.*
%
% see group_fit_model, group_brain_nlfit, nl_contrasts

function EXPT = get_htw_image_names(EXPT, varargin)
    EXPT.NLCON.height = {};
    EXPT.NLCON.delay = {};
    EXPT.NLCON.width = {};

    if length(varargin) > 0 % get first n images
        n = varargin{1};
    end

    for i = 1:length(EXPT.subjects)

        cd(EXPT.subjects{i})

        d = dir('height*.img');

        if isempty(d)
            disp(['NO IMAGES FOR ' EXPT.subjects{i} ': SKIPPING.'])
        else    % we have images
            if length(varargin) > 0, d = d(1:n); ,end  % get first n images
            % store these in NLCON structure for random effects or contrast image
            % calculation
            for j = 1:length(d)

                if i == 1
                    EXPT.NLCON.height{j} = [pwd filesep d(j).name];
                else
                    EXPT.NLCON.height{j} = str2mat(EXPT.NLCON.height{j},[pwd filesep d(j).name]);
                end

            end

            % store these in NLCON structure for random effects
            d = dir('delay*.img');
            if length(varargin) > 0, d = d(1:n); ,end  % get first n images

            for j = 1:length(d)

                if i == 1
                    EXPT.NLCON.delay{j} = [pwd filesep d(j).name];
                else
                    EXPT.NLCON.delay{j} = str2mat(EXPT.NLCON.delay{j},[pwd filesep d(j).name]);
                end

            end

            % store these in NLCON structure for random effects
            d = dir('width*.img');
            if length(varargin) > 0, d = d(1:n); ,end  % get first n images

            for j = 1:length(d)

                if i == 1
                    EXPT.NLCON.width{j} = [pwd filesep d(j).name];
                else
                    EXPT.NLCON.width{j} = str2mat(EXPT.NLCON.width{j},[pwd filesep d(j).name]);
                end

            end

        end     % if no images

        cd ..
    end
end