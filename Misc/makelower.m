function makelower

d = dir;

for i = 1:length(d)
    
    if strcmp(d(i).name,upper(d(i).name)) & ~d(i).isdir
        
        str = ['!mv ' d(i).name ' ' lower(d(i).name)];
        
        disp(str)
        eval(str)
    end
end

return
