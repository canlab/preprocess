
for i = 3:length(EXPT.SUBJECT.contrasts)
    
    gvol = group_outline(EXPT,.05,8,i,[3 5 7],0);
        cd ..
        
    % 1st 2 contrasts are all effects; start with 3rd
    % names lag one behind actual con numbers
    
    title(EXPT.SUBJECT.connames{i-1})
    saveas(gcf,['group_' EXPT.SUBJECT.connames{i-1}],'jpg')
    close
    

end