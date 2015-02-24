% function EXPT = read_ind_results(EXPT)
% Tor Wager
%
% Gets image, SPM.mat, and other info from individual subject directories

function EXPT = read_ind_results(EXPT)
    if ~isfield(EXPT,'SNPM'), EXPT.SNPM.P{1} = [];, end
    if iscell(EXPT.SNPM.P)
        if ~iscell(EXPT.SNPM.P{1}), EXPT.SNPM.P = [];, end
    end

    EXPT.snums = [];
    subindex = 0;
    getfirst = 1;

    D = dir;

    do_r = input('Press 0 to get regular con*imgs or 1 to get Rcon*imgs (residuals from covariate): ');

    for i = 3:length(D)

        if D(i).isdir & ~(strcmp(D(i).name(1:min(3,length(D(i).name))),'rfx') | strcmp(D(i).name(1:min(4,length(D(i).name))),'snpm') | strcmp(D(i).name(1:min(3,length(D(i).name))),'fix'))

            eval(['cd ' D(i).name])
            clear whch

            if exist('SPM.mat') == 2

                load SPM.mat

                % -------------------------------------------------------
                % * get stuff that is the same for all subjects (usually)
                % -------------------------------------------------------
                if getfirst
                    EXPT.nsess = str2num(xsDes.Number_of_sessions);
                    EXPT.regnames = Sess{1}.name;
                    EXPT.HP = xX.K{1}.HParam;
                    EXPT.nscans = size(xX.X,1);
                    EXPT.eres = size(Sess{1}.sf{1},1) ./ length(Sess{1}.pst{1});
                    EXPT.TR = xX.RT;
                    a = diag(xM.VM.mat);
                    EXPT.voxSize = a(1:3);
                end




                % -------------------------------------------------------
                % * get stuff that we need from each subject
                % -------------------------------------------------------

                mysub = subindex + 1;
                subindex = subindex + 1;
                EXPT.snums = [EXPT.snums mysub];
                EXPT.subjects{subindex} = D(i).name;
                EXPT.im_files{mysub} = cat(1,VY.fname);
                EXPT.an_files{mysub} = [pwd filesep 'SPM.mat'];
                [Files,Dirs] = spm_list_files(pwd,'beta*img');
                EXPT.betas{mysub} = Files;

                b = [];
                for j = 1:EXPT.nsess
                    a = cell2mat(cat(1,Sess{j}.sf));
                    b = [b;a];
                end
                EXPT.evts_sf{mysub} = b;

                % -------------------------------------------------------
                % contrast info and filenames for random effects
                % -------------------------------------------------------

                D2 = dir;
                a = str2mat(1,D2(cat(1,D2.isdir) == 0).name);
                for j = 1:size(a,1),
                    a2 = deblank(a(j,:));,
                    % con*imgs

                    if do_r
                        whch(j) = (strcmp(a2(1:min(length(a2),5)),'Rcon_') & strcmp(a2(max(1,end-2):end),'img'));,
                    else
                        whch(j) = (strcmp(a2(1:min(length(a2),4)),'con_') & strcmp(a2(max(1,end-2):end),'img'));,
                    end
                    % height, delay, intercept - if present - WRONG.
                    % nlcon*img defined by nl_contrasts.m
                    %whchnl_h(j) = (strcmp(a2(1:min(length(a2),5)),'cond_') & strcmp(a2(max(1,end-5):end),'ht.img'));,
                    %whchnl_d(j) = (strcmp(a2(1:min(length(a2),5)),'cond_') & strcmp(a2(max(1,end-5):end),'ay.img'));,
                    %whchnl_i(j) = (strcmp(a2(1:min(length(a2),5)),'cond_') & strcmp(a2(max(1,end-5):end),'pt.img'));,
                end
                nl = a;
                a = a(whch,:);
                %clear b % all this stuff is wrong.
                %b{1} = nl(whchnl_h,:);
                %b{2} = nl(whchnl_d,:);
                %b{3} = nl(whchnl_i,:);

                if getfirst
                    if do_r
                        EXPT.SNPM.connums = str2num(a(:,6:9))'; % same for everyone
                    else
                        EXPT.SNPM.connums = str2num(a(:,5:8))'; % same for everyone
                    end

                    load xCon.mat
                    a2 = str2mat(xCon.name);
                    EXPT.SNPM.connames = a2(EXPT.SNPM.connums,:);
                    getfirst = 0;
                end

                for j = 1:size(a,1)
                    EXPT.SNPM.P{j}{subindex} = [pwd filesep deblank(a(j,:))];
                end

                % nonlinear fit filenames
                %for k = 1:length(b)
                %    for j = 1:size(b{1},1)
                %        EXPT.DX.nlconP{k}{j}{subindex} = [pwd filesep deblank(b{k}(j,:))];
                %    end
                %end

            end % load stuff

            cd ..

        end % if is a directory

    end % loop through entries in Dir

    % format contrast image names
    for i = 1:length(EXPT.SNPM.P)
        EXPT.SNPM.P{i} = str2mat(EXPT.SNPM.P{i});
    end

    %for i = 1:length(EXPT.DX.nlconP)
    %    for j = 1:length(EXPT.DX.nlconP{1})
    %        EXPT.DX.nlconP{i}{j} = str2mat(EXPT.DX.nlconP{i}{j});
    %    end
    %end

    EXPT.FILT.HP = EXPT.HP;
    EXPT.FILT.TR = EXPT.TR;
    EXPT.FILT.filtertype = 'spm';
    EXPT.FILT.S = use_spm_filter(EXPT.TR,size(EXPT.im_files{1},1),'none','specify',EXPT.HP);

end
