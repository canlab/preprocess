slist = {'sub1' 'sub2' 'sub4' 'sub6' 'sub8' 'sub9' 'sub10' 'sub11'};
index = 1;

for sub = {'s1' 's2' 's4' 's6' 's8' 's9' 's10' 's11'}
    
    eval(['lt = ' sub{1} 'lattimes;'])
    eval(['st = ' sub{1} 'times(1:size(lt,1),:);'])
    
    s = getscans_intext(st,'test',lt);
    
    eval([slist{index} ' = s;'])
    
    index = index + 1;
    
end
    