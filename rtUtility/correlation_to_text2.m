function str = correlation_to_text(a,cm,varargin)
%str = correlation_to_text(matrix,value to mark with *,[cell array of names],[robust flag])

warning('Function is deprecated; use correlation_to_text.m');

if length(varargin) > 0, nms = varargin{1};,else, nms = {[]};,end

% names

str=sprintf('Crit=%3.2f\t',cm);
for j=1:size(a,2),
    
    if length(nms) < j, nms{j} = ['V' num2str(j)];,end
    
    str=[str sprintf('%s\t',nms{j})];,
end
str=[str sprintf('\n')];,
    

% table

for i = 1:size(a,1),
    
    str=[str sprintf('%s\t',nms{i})];,
    
    for j=1:i,
        
        if abs(a(i,j))>cm & a(i,j) ~= 1,t='*';,
            
        else,t='';,
            
        end,
        
        str=[str sprintf('%3.3f%s\t',a(i,j),t)];,
    end,
    str=[str sprintf('\n')];,
end

disp(str)

return
