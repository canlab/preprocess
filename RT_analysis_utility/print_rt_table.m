function print_rt_table(dat,varargin)

% input: a i x j cell array,where each cell contains numbers or text for n trials
% optional: cell array of column names
% output: a text table in tab delimited format
% tor wager

% Print names
if length(varargin) > 0, 
    for i = 1:length(varargin{1})
        fprintf(1,'%s\t',varargin{1}{i}), 
    end
    fprintf(1,'\n')
end


% print table
for i = 1:size(dat,1)
    % for each row of cells
    
    print_cell_row(dat(i,:))
    
end

return



% sub-functions


function print_cell_row(dat)

ntrials = size(dat{1,1},1);


for rowind = 1:ntrials
    
    % print a row of data
    for j = 1:size(dat,2)
    
        tmp = dat{1,j};
        
        if ischar(tmp(rowind))
            % text
            fprintf(1,'%s\t',tmp(rowind));
            
            
        else
            % a number
            fprintf(1,'%3.4f\t',tmp(rowind));
            
        end
        
    end
    
    fprintf(1,'\n')
    
    
end % loop thru rowind