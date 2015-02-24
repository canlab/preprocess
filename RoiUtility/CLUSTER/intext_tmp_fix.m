function newcl = intext_tmp_fix(c2m,subcl)
%subclusters = intext_tmp_fix(clusters_to_match,subclusters_to_change)

for i = 1:length(subcl)
    
    for j = 1:length(subcl{i})
        
            % make all field names match.

            newcl{i}(j) = c2m(1);
            N = fieldnames(newcl{i}(j));
            for NN = 1:length(N)
                eval(['newcl{i}(j).' N{NN} ' = [];'])
            end
            
            %N = fieldnames(clusters(i));
            for NN = 1:length(N)
                eval(['newcl{i}(j).' N{NN} ' = subcl{i}(j).' N{NN} ';'])
            end
            
        end
        
        
    end
    
    return
    