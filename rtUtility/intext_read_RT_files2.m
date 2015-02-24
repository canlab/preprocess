% to continue loading new subjects:
% load old DATA.mat file
% make sure subject numbers are the same for new and old, no duplicates
% start with the index number of the new subjects
% run.
%
% type:
% dd = dir('*-1.txt'); substart = 181; subend = length(dd); str2mat(dd.name)
% intext_read_rt_files2
%
% for subject 850
% substart = 202; subend = 202;
%
% try:
% dd = dir('*-1.txt'); [num2str((1:length(dd))') str2mat(dd.name)]
%
% REMEMBER: YOU MUST BE RUNNING MATLAB 6.1 OR HIGHER!!!!
% also try:
% N = fieldnames(RT); for i = 4:length(N),for j = 1:size(RT.j2_ext,2), eval(['a(i,j) = nanmean(RT.' N{i} '{j});']),end,end


cd C:\Tor_Documents\CurrentExperiments\Intext\from_server\eprime_data_behavioral

fname = '850-1.txt';
numc = 201;

dd = dir('*-1.txt');

for subj = substart:subend %length(dd)   % 143 starts new batch Dec. 2002
    
    if strcmp(dd(subj).name,'114-1.txt') | strcmp(dd(subj).name,'120-1.txt'), numc = 166;
    elseif strcmp(dd(subj).name,'231-1.txt'), numc = 61;
    elseif strcmp(dd(subj).name,'300-1.txt'), numc = 196;
    else numc = 201;
    end
    
    if exist('RT') == 1, endloc = length(RT.subjnum) + 1;, else, endloc = 1;, end
    
            RT.subjnum(endloc) = str2num(dd(subj).name(1:3));
            ACC.subjnum(endloc) = RT.subjnum(endloc);
            
    % d and other vars are reserved in read_edat_output
    fname = dd(subj).name;
    fprintf(1,' %s',fname)
    
    clear judge* Attr* trial* Stay* stim* stay* An* color* prac* ext* Sess*
    try
        read_edat_output
        success = 1;
    catch
        disp([' Problem reading ' dd(subj).name])
        success = 0;
        keyboard
    end
    
    if success
        try
            

    % -------------------------------------------------------------------------------------------------------------
    % Accuracy
    % -------------------------------------------------------------------------------------------------------------
    
    ACC.ext1(endloc,1) = sum(judge3_ACCTrial == 1) ./ (sum(judge3_ACCTrial == 0) + sum(judge3_ACCTrial == 1));
    ACC.ext2(endloc,1) = sum(judge4_ACCTrial == 1) ./ (sum(judge4_ACCTrial == 0) + sum(judge4_ACCTrial == 1));
    ACC.int1(endloc,1) = sum(judge1_ACCTrial == 1) ./ (sum(judge1_ACCTrial == 0) + sum(judge1_ACCTrial == 1));
    ACC.int2(endloc,1) = sum(judge2_ACCTrial == 1) ./ (sum(judge2_ACCTrial == 0) + sum(judge2_ACCTrial == 1));
    
    extacc = judge3_ACCTrial == 1 & judge4_ACCTrial == 1; extacc(isnan(judge3_ACCTrial) | isnan(judge4_ACCTrial)) = NaN;
    intacc = judge1_ACCTrial == 1 & judge2_ACCTrial == 1; intacc(isnan(judge1_ACCTrial) | isnan(judge2_ACCTrial)) = NaN;
    
    ACC.ext_trial(endloc,1) = sum(extacc == 1) ./ (sum(extacc == 0) + sum(extacc == 1));
    ACC.int_trial(endloc,1) = sum(intacc == 1) ./ (sum(intacc == 0) + sum(intacc == 1));
    
    ACC.ext_nosw(endloc,1) = sum(extacc == 1 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Stay.bmp')); 
    ACC.ext_nosw(endloc,1) = ACC.ext_nosw(endloc,1) ./ (ACC.ext_nosw(endloc,1) + sum(extacc == 0 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Stay.bmp'))); 
    ACC.ext_attrsw(endloc,1) = sum(extacc == 1 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Stay.bmp')); 
    ACC.ext_attrsw(endloc,1) = ACC.ext_attrsw(endloc,1) ./ (ACC.ext_attrsw(endloc,1) + sum(extacc == 0 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Stay.bmp')));
    ACC.ext_obsw(endloc,1) = sum(extacc == 1 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Switch.bmp')); 
    ACC.ext_obsw(endloc,1) = ACC.ext_obsw(endloc,1) ./ (ACC.ext_obsw(endloc,1) + sum(extacc == 0 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Switch.bmp')));
    ACC.ext_bothsw(endloc,1) = sum(extacc == 1 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Switch.bmp')); 
    ACC.ext_bothsw(endloc,1) = ACC.ext_bothsw(endloc,1) ./ (ACC.ext_bothsw(endloc,1) + sum(extacc == 0 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Switch.bmp')));
    
    ACC.int_nosw(endloc,1) = sum(intacc == 1 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Stay.bmp')); 
    ACC.int_nosw(endloc,1) = ACC.int_nosw(endloc,1) ./ (ACC.int_nosw(endloc,1) + sum(intacc == 0 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Stay.bmp'))); 
    ACC.int_attrsw(endloc,1) = sum(intacc == 1 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Stay.bmp')); 
    ACC.int_attrsw(endloc,1) = ACC.int_attrsw(endloc,1) ./ (ACC.int_attrsw(endloc,1) + sum(intacc == 0 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Stay.bmp')));
    ACC.int_obsw(endloc,1) = sum(intacc == 1 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Switch.bmp')); 
    ACC.int_obsw(endloc,1) = ACC.int_obsw(endloc,1) ./ (ACC.int_obsw(endloc,1) + sum(intacc == 0 & strcmp(AttrSwitch,'Same') & strcmp(Attribute3Trial,'Switch.bmp')));
    ACC.int_bothsw(endloc,1) = sum(intacc == 1 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Switch.bmp')); 
    ACC.int_bothsw(endloc,1) = ACC.int_bothsw(endloc,1) ./ (ACC.int_bothsw(endloc,1) + sum(intacc == 0 & strcmp(AttrSwitch,'Different') & strcmp(Attribute3Trial,'Switch.bmp')));
      
    % -------------------------------------------------------------------------------------------------------------
    % RT collection
    % -------------------------------------------------------------------------------------------------------------  
    
    RT.ext_trials(endloc) = sum(~isnan(extacc));
    RT.int_trials(endloc) = sum(~isnan(intacc));
    
    RT.stim_ext{endloc} = stim_RTTrial; RT.stim_ext{endloc}(~(extacc==1)) = NaN; 
    RT.stim_int{endloc} = stim_RTTrial; RT.stim_int{endloc}(~(intacc==1)) = NaN; 
    RT.color_ext{endloc} = color1_RTTrial; RT.color_ext{endloc}(~(extacc==1)) = NaN; 
    RT.color_int{endloc} = color_RTTrial; RT.color_int{endloc}(~(intacc==1)) = NaN; 
    
    RT.j1_ext{endloc} = judge3_RTTrial; RT.j1_ext{endloc}(extacc == 0) = NaN;
    RT.j1_int{endloc} = judge1_RTTrial; RT.j1_int{endloc}(intacc == 0) = NaN;
    RT.swcue_ext{endloc} = stayswit1_RTTrial; RT.swcue_ext{endloc}(extacc == 0) = NaN;
    RT.swcue_int{endloc} = stayswit_RTTrial; RT.swcue_int{endloc}(intacc == 0) = NaN;
    RT.j2_ext{endloc} = judge4_RTTrial; RT.j2_ext{endloc}(extacc == 0) = NaN;
    RT.j2_int{endloc} = judge2_RTTrial; RT.j2_int{endloc}(intacc == 0) = NaN;
    
    RT.j2_ext_nosw{endloc} = RT.j2_ext{endloc};  RT.j2_ext_nosw{endloc}(~(strcmp(AttrSwitch,'Same')) | ~(strcmp(Attribute3Trial,'Stay.bmp'))) = NaN;
    RT.j2_ext_attrsw{endloc} = RT.j2_ext{endloc};  RT.j2_ext_attrsw{endloc}(~(strcmp(AttrSwitch,'Different')) | ~(strcmp(Attribute3Trial,'Stay.bmp'))) = NaN;
    RT.j2_ext_obsw{endloc} = RT.j2_ext{endloc};  RT.j2_ext_obsw{endloc}(~(strcmp(AttrSwitch,'Same')) | ~(strcmp(Attribute3Trial,'Switch.bmp'))) = NaN;
    RT.j2_ext_bothsw{endloc} = RT.j2_ext{endloc};  RT.j2_ext_bothsw{endloc}(~(strcmp(AttrSwitch,'Different')) | ~(strcmp(Attribute3Trial,'Switch.bmp'))) = NaN;
    
    RT.swcue_ext_nosw{endloc} = RT.swcue_ext{endloc};  RT.swcue_ext_nosw{endloc}(~(strcmp(Attribute3Trial,'Stay.bmp'))) = NaN;
    RT.swcue_ext_obsw{endloc} = RT.swcue_ext{endloc};  RT.swcue_ext_obsw{endloc}(~(strcmp(Attribute3Trial,'Switch.bmp'))) = NaN;
    
    RT.j2_int_nosw{endloc} = RT.j2_int{endloc};  RT.j2_int_nosw{endloc}(~(strcmp(AttrSwitch,'Same')) | ~(strcmp(Attribute3Trial,'Stay.bmp'))) = NaN;
    RT.j2_int_attrsw{endloc} = RT.j2_int{endloc};  RT.j2_int_attrsw{endloc}(~(strcmp(AttrSwitch,'Different')) | ~(strcmp(Attribute3Trial,'Stay.bmp'))) = NaN;
    RT.j2_int_obsw{endloc} = RT.j2_int{endloc};  RT.j2_int_obsw{endloc}(~(strcmp(AttrSwitch,'Same')) | ~(strcmp(Attribute3Trial,'Switch.bmp'))) = NaN;
    RT.j2_int_bothsw{endloc} = RT.j2_int{endloc};  RT.j2_int_bothsw{endloc}(~(strcmp(AttrSwitch,'Different')) | ~(strcmp(Attribute3Trial,'Switch.bmp'))) = NaN;
    
    RT.swcue_int_nosw{endloc} = RT.swcue_int{endloc};  RT.swcue_int_nosw{endloc}(~(strcmp(Attribute3Trial,'Stay.bmp'))) = NaN;
    RT.swcue_int_obsw{endloc} = RT.swcue_int{endloc};  RT.swcue_int_obsw{endloc}(~(strcmp(Attribute3Trial,'Switch.bmp'))) = NaN;
    
    catch
        disp(['trying alternate loading'])
        %try
            alternate_read_RT2_loading
            disp('Success!')
            %catch
            %disp([' Problem compiling ' dd(subj).name])
            %end
    end

    end % if success

    % edat output for trial dependencies
    
    EXP.subjnum(endloc,1) = str2num(dd(subj).name(1:3));
     
    Obsw = Attribute3Trial;Obsw(strcmp(Obsw,'?')) = {'NaN'};
    Attsw = AttrSwitch; Attsw(strcmp(Attsw,'?')) = {'NaN'};
    EXP.Obsw{endloc} = Obsw;
    EXP.Attsw{endloc} = Attsw;
    
    if ~exist('judge4_RTTrial') == 1, judge4_RTTrial = judge4_RT;,end
    alltrials = max([judge4_RTTrial judge2_RTTrial]')';
    [gls,desc] = rt_analyzer(alltrials,Attsw,Obsw,'nback',1);
    EXP.bs(subj,:) = gls.betas;
    try,EXP.means(subj,:) = desc.means;,catch,EXP.means(subj,:) = NaN;, end
    try,EXP.counts(subj,:) = desc.count;,catch, EXP.counts(subj,:) = NaN;, end
    
    [gls,desc] = rt_analyzer(judge4_RTTrial,Attsw,Obsw,'nback',1);
    EXP.Gext_bs(subj,:) = gls.betas;
    try,EXP.Gext_means(subj,:) = desc.means;,catch,EXP.Gext_means(subj,:) = NaN;, end
    try,EXP.Gext_counts(subj,:) = desc.count;,catch, EXP.Gext_counts(subj,:) = NaN;, end
        
    [gls2,desc2] = rt_analyzer(judge2_RTTrial,Attsw,Obsw,'nback',1);
    EXP.Gint_bs(subj,:) = gls2.betas;
    try,EXP.Gint_means(subj,:) = desc2.means;,catch,EXP.Gint_means(subj,:) = NaN;, end
    try,EXP.Gint_counts(subj,:) = desc2.count;,catch,EXP.Gint_counts(subj,:) = NaN;, end
    
end % loop

fprintf(1,'\n')

% remove duplicates!
disp('Removing duplicates, but NOT from EXP!')
[tmp,w] = unique(RT.subjnum); ww = zeros(size(RT.subjnum)); ww(w) = 1; w = ~ww;  w=logical(w);
N = fieldnames(RT);
for i = 1:length(N)
    fprintf(1,'%s ',N{i})
    eval(['tmp = RT.' N{i} ';'])
    eval(['RT.' N{i} '(w) = [];'])
end

[tmp,w] = unique(ACC.subjnum); ww = zeros(size(ACC.subjnum)); ww(w) = 1; w = ~ww;  w=logical(w);
N = fieldnames(ACC);
for i = 1:length(N)
    fprintf(1,'%s ',N{i})
    eval(['tmp = ACC.' N{i} ';'])
    eval(['ACC.' N{i} '(w) = [];'])
end

EXP.gls = gls;
EXP.desc = desc;

% get Accuracy, and do the processing only if accurate
% -------------------------------------------------------------------
clear acc
N = fieldnames(ACC); eval(['acc = [ACC.' N{4} '];']);
for i = 5:length(N)-2, eval(['acc = [acc ACC.' N{i} '];']);,end
acc1 = min(acc')';
acc2 = mean(acc')';
ACC.minacc = acc1;   % something messed up with ACC.ext2 in a few subjects, so we won't use this - use individual acc values
ACC.meanacc = acc2;
ACC.include = acc1>=.75;
ACC.include80 = acc2>.8;

fprintf(1,'Included %3.0f, excluded %3.0f\n',sum(ACC.include),sum(ACC.include==0))



    