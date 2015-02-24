function um_run_check_norm(EXPT)

eval(['cd ' EXPT.SUBJECT.studydir])

for i = 1:length(EXPT.subjects)
    
    EXPT.SUBJECT.subjcode = EXPT.subjects{1};
    
    um_check_norm(EXPT.SUBJECT)
    
end