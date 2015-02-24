function [str,a] = print_correlation(a,cm,varargin)
%[str,corrmtx] = print_correlation(matrix,value to mark with *,[cell array of names],[pairwise delete])
%
% Input options: 
% 1) matrix is correlation matrix, value is critical value for *'s
% input correlation matrix and critical value (assumes correlations are entered)
% 
% 2) raw data and empty value (runs correlations on matrix within program)
% (determines threshold based on n)
%
% see also print_correlation.m

warning('Function is deprecated; use correlation_to_text.m');

cm2 = Inf; pairwise = 0;

if length(varargin) > 0, nms = varargin{1};,else, nms = {[]};,end
if length(varargin) > 1, pairwise = varargin{2};, end 
    
% NaNs
wh = find(any(isnan(a),2));
if wh,
    if pairwise
        disp('Warning: NaNs! Removing casewise.'),

        a(wh,:) = [];
    else
        % not done yet.
        disp('Warning: NaNs! Removing casewise.'),

        a(wh,:) = [];

    end
end


if isempty(cm)
    % choose default value
    [rci,sig,z,p,cm] = r2z(.5,size(a,1),.05);
    [rci,sig,z,p,cm2] = r2z(.5,size(a,1),.1);
    a = corrcoef(a);
end


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
    
    for j=1:size(a,2),
        
        if abs(a(i,j))>cm & a(i,j) ~= 1,t='*';
        elseif abs(a(i,j))>cm2 & a(i,j) ~= 1,t='+';
            
        else t='';
            
        end,
        
        str=[str sprintf('%3.3f%s\t',a(i,j),t)];
    end,
    str=[str sprintf('\n')];
end

disp(str)

return
