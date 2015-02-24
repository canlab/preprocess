function EXPT = getfunctnames(EXPT)
%EXPT = getfunctnames(EXPT)
%
% Get functional image names 
%

OPT.object = 'het1spgr.img';	

if ~isfield(EXPT,'wildcard'), EXPT.wildcard = 'swMC*.img';,end

wild1 = EXPT.wildcard(1:find(EXPT.wildcard == '*')-1);
OPT.wildcard = EXPT.wildcard;

eval(['cd ' EXPT.studydir])
        if ~isfield(EXPT,'nsess'),EXPT.nsess = input('How many runs?');,end

        for i = 1:length(EXPT.subjects)        
            clear scan1dir
            for j = 1:EXPT.nsess, scan1dir{i}{j} = [EXPT.subjects{i} filesep 'scan' num2str(j)];,
                EXPT.im_files{i} = str2mat(tor_list_files(scan1dir{i},EXPT.wildcard));
	            EXPT.nravols{i} = str2mat(tor_list_files(scan1dir{i},EXPT.wildcard(2:end)));
            
            end

            clear scan1dir
            for i = 1:length(EXPT.subjects) 
                scan1dir{i} = [EXPT.subjects{i} filesep 'scan1'];
                anatdir{i} = [EXPT.subjects{i} filesep 'anatomy'];
            end
            
            % save stuff in EXPT for later reference
            OPT.firstfun = str2mat(tor_list_files(scan1dir,[wild1 '*0001.img']));
            OPT.firstnra = str2mat(tor_list_files(scan1dir,[wild1(2:end) '*0001.img']));
            OPT.firstra = str2mat(tor_list_files(scan1dir,[wild1(3:end) '*0001.img']));
            OPT.nspgr = str2mat(tor_list_files(anatdir,['n' OPT.object]));
        end

EXPT.SUBJECT = OPT;
EXPT.OPT = OPT;

return

% simpler, faster way
%for s = 1:length(EXPT.subjects)

%wcard = 'w*img';

%dd = [];   	p = pwd;

%for i = 1:4 
%	d = dir([EXPT.subjects{s} '/scan' num2str(i) '/' wcard]); d = str2mat(d.name); 
	
%	dp = fullfile(p,EXPT.subjects{s},['scan' num2str(i)],filesep);
%	dp = repmat(dp,size(d,1),1);
%	d  = [dp d];

%	dd = [dd; d];
%end

%EXPT.nravols{s} = dd;

end

