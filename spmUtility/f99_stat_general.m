function f99_stat_general(RT,nscan,rep,n_cond,cond_name,LULRES1)

% RT          - repitition time
% nscan       - number of scans per session
% rep         - are sessions replicated exactly (1=yes / 0=no)
% n_cond      - number of conditions or trials
% cond_name   - matrix with names of conditions/trials
% RES         - directory where to write results


global fmriTEST kulRES fmriDIR
kulRES = LULRES1;
% Change to results directory
%----------------------------
if f99_exist(fmriDIR,'RESULTS') == 0
  str = ['!mkdir ' fmriDIR filesep 'RESULTS'];
  eval(str);
end
str = [fmriDIR filesep 'RESULTS' filesep];
if f99_exist(str,kulRES) == 0
  str = ['!mkdir ' str kulRES];
  eval (str);
end

str = ['cd ' fmriDIR filesep 'RESULTS' filesep kulRES];
str = deblank(str);
eval (str)


logfile=[fmriDIR filesep 'LOG' filesep 'design.log'];




st_rep = '';
if rep == 0
    st_rep = 'not ';
end
disp('');
disp('*****************************************************************');
f99_log(logfile,['OPTIONS FOR DESIGN MATRIX']);
disp('-----------------------------------------------------------------');
f99_log(logfile,['GENERAL']);
f99_log(logfile,['  Design consists of ' mat2str(size(nscan,2)) ' sessions']);
f99_log(logfile,['    with number of scans ' mat2str(nscan)]);
f99_log(logfile,['    and sessions are ' st_rep 'replicated exactly']);
f99_log(logfile,['  Repetition Time is ' mat2str(RT)]);
f99_log(logfile,['  The ' mat2str(n_cond) ' trial names are']);
for i = 1:n_cond
    f99_log(logfile,['        ' cond_name(i,:)]);
end

st_rep = '';
if rep == 0
    st_rep = 'not ';
end
f99_log(logfile,['  Results are written in the directory : ' kulRES]);

filename = [fmriDIR filesep 'RESULTS' filesep kulRES filesep 'KUL_stat_gen'];
str = ['save ' filename ' RT nscan rep n_cond cond_name kulRES'];
eval(str)
