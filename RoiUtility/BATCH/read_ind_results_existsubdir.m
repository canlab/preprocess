function EXPT = read_ind_results(EXPT)
% function EXPT = read_ind_results(EXPT)
% Tor Wager
%
% Gets image, SPM.mat, and other info from individual subject directories
% Start in dir above ind subject dirs

    if ~isfield(EXPT,'SNPM'), EXPT.SNPM.P{1} = [];, end
    if iscell(EXPT.SNPM.P)
        if ~iscell(EXPT.SNPM.P{1}), EXPT.SNPM.P = [];, end
    end
    
EXPT.snums = [];
subindex = 0;
getfirst = 1;

do_r = input('Press 0 to get regular con*imgs or 1 to get Rcon*imgs (residuals from covariate): ');
do_t = input('Get T-images as well as con*imgs? (1/0) ');
if do_t
    if ~isfield(EXPT.SNPM,'T'), EXPT.SNPM.T = {};,end
end

for i = EXPT.subjects 	% 3:length(D)
    
        
        try, eval(['cd ' i{1}]), catch, disp(['No directory called '  i{1}]),end 
        clear whch
        
        if exist('SPM.mat') == 2
            
            load SPM.mat
            
            % -------------------------------------------------------
            % * get stuff that is the same for all subjects (usually)
            % -------------------------------------------------------
            if getfirst
                
                % SPM2 fix
                if ~(exist('xsDes') == 1), xsDes = SPM.xsDes;, end
                if ~(exist('xX') == 1), xX = SPM.xX;, end
                if ~(exist('Sess') == 1), Sess = SPM.Sess;, end
                if ~iscell(Sess), tmp = Sess; Sess = {}; Sess{1} = tmp;,end
                if ~isfield(Sess,'name'), Sess{1}.name = SPM.condnames;,end
                
                try
                    % spm99
                    EXPT.HP = xX.K{1}.HParam;
                catch
                    % spm2
                    EXPT.HP = xX.K(1).HParam;
                end
                EXPT.nsess = str2num(xsDes.Number_of_sessions);
                EXPT.regnames = Sess{1}.name;

                EXPT.nscans = size(xX.X,1);
                
                try
                    EXPT.TR = xX.RT;
                catch
                    % spm2
                    EXPT.TR = SPM.xY.RT;
                end
                
                try
                    EXPT.eres = size(Sess{1}.sf{1},1) ./ length(Sess{1}.pst{1});
                catch
                    %spm2, i think
                    %size(Sess{1}.U(1).u,1) ./ length(Sess{1}.U(1).pst);
                    EXPT.eres = EXPT.TR ./ Sess{1}.U(1).dt;
                end

                try
			a = diag(xM.VM.mat);
            EXPT.voxSize = a(1:3);
		catch
			warning('no VM field in xM? Don''t know voxel sizes')
		end
            end
            


            
            % -------------------------------------------------------
            % * get stuff that we need from each subject
            % -------------------------------------------------------
            
            mysub = subindex + 1;
            subindex = subindex + 1;
            EXPT.snums = [EXPT.snums mysub];
            EXPT.subjects_check{subindex} = i{1};
            
            % SPM2 fix
            if ~(exist('VY') == 1), VY = SPM.xY.VY;, end
            
            EXPT.im_files{mysub} = cat(1,VY.fname);
            EXPT.an_files{mysub} = [pwd filesep 'SPM.mat'];
            try
                [Files,Dirs] = spm_list_files(pwd,'beta*img');
            catch
                Ptmp = dir('beta*img');
                Files = str2mat(Ptmp.name);
            end
            EXPT.betas{mysub} = Files;
            
            b = [];
            for j = 1:EXPT.nsess
                try
                    a = cell2mat(cat(1,Sess{j}.sf));
                catch
                    %spm2
                    a = cat(1,Sess{1}.U(:).u);
                end
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
                    
                    try, load xCon.mat
                    catch
                        % spm2
                        xCon = SPM.xCon;
                    end
                    
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
                
                if do_t
                    % -------------------------------------------------------
                    % contrast info and filenames for random effects
                    % -------------------------------------------------------
                
                    D2 = dir;
                    a = str2mat(1,D2(cat(1,D2.isdir) == 0).name);
                    for j = 1:size(a,1), 
                        a2 = deblank(a(j,:));, 
                        whch(j) = (strcmp(a2(1:min(length(a2),5)),'spmT_') & strcmp(a2(max(1,end-2):end),'img'));,

                    end
                    nl = a;
                    a = a(whch,:);

                    for j = 1:size(a,1)
                        EXPT.SNPM.T{j}{subindex} = [pwd filesep deblank(a(j,:))];
                    end 
                end
                
                
                
                
 	else, disp(['Warning: no SPM.mat file for ' i{1}])     ,pwd 
	end % load stuff
        
        cd ..
        
end % loop through entries in EXPT.subjects

% format contrast image names
for i = 1:length(EXPT.SNPM.P)
    try,EXPT.SNPM.P{i} = str2mat(EXPT.SNPM.P{i});
	keepit(i) = 1;
	if isempty(EXPT.SNPM.P{i}), keepit(i) = 0;,end
   catch,
	warning(['Problem with contrast ' num2str(i)])
	EXPT.SNPM.P{i} = [];
	keepit(i) = 0;
   end
end

EXPT.SNPM.P = EXPT.SNPM.P(find(keepit));


%for i = 1:length(EXPT.DX.nlconP)
%    for j = 1:length(EXPT.DX.nlconP{1})
%        EXPT.DX.nlconP{i}{j} = str2mat(EXPT.DX.nlconP{i}{j});
%    end
%end

if do_t
    % format contrast image names
    for i = 1:length(EXPT.SNPM.T)
        try,EXPT.SNPM.T{i} = str2mat(EXPT.SNPM.T{i});
        keepit(i) = 1;
        if isempty(EXPT.SNPM.T{i}), keepit(i) = 0;,end
        catch,
            warning(['Problem with contrast ' num2str(i)])
            EXPT.SNPM.T{i} = [];
            keepit(i) = 0;
        end
    end
    
    EXPT.SNPM.T = EXPT.SNPM.T(find(keepit));
end






 EXPT.FILT.HP = EXPT.HP;
 EXPT.FILT.TR = EXPT.TR;
 EXPT.FILT.filtertype = 'spm';
 EXPT.FILT.S = use_spm_filter(EXPT.TR,size(EXPT.im_files{1},1),'none','specify',EXPT.HP);
 
 return
 