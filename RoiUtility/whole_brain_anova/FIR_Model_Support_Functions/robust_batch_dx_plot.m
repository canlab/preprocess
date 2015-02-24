function EXPT = robust_batch_dx_plot(EXPT,actionstr,varargin)
% EXPT = robust_batch_dx_plot(EXPT,actionstr,[startat directory])
%
% actionstr = 'save' 'plot' 'both'

% set up necessary information
EXPT = CreateExpt('extractdx');


if length(varargin) > 0, startat = varargin{1};, else, startat = 1;, end

d = dir('robust*');
d = d(cat(1,d.isdir));

if strcmp(actionstr,'save') | strcmp(actionstr,'both')

for i = startat:length(d)
    
    fprintf(1,'\n*---------\n%s\n',d(i).name)
    
    cd(d(i).name)
    
    d2 = dir('robust*mat');
    
    for j = 1:length(d2)
        
        fprintf(1,'\t%s\n',d2(j).name)
        clear cl
        load(d2(j).name)
        
        if exist('cl','var') && ~isempty(cl)
            [cl,EXPT] = extract_dxbeta_data(EXPT,cl,0);
        end
        
        eval(['save ' d2(j).name ' cl'])
        
    end
    
    cd ..
    
end

end     % save part




if strcmp(actionstr,'plot') | strcmp(actionstr,'both')

    name = 'rob*';
    %name = input('Type name of mat file containing cl variable (clusters), no extension, e.g., robust00*_intercept_pos: ','s');
    name = [name '.mat'];
 
    if ~isfield(EXPT,'FIR'),EXPT.FIR = [];,end
     
    if isfield(EXPT.FIR,'indiv')
        doindiv = EXPT.FIR.indiv;
    else
    
        doindiv = input('Plot whole group (0) or individual diffs high/low split (1): ');
    end

    if isfield(EXPT.FIR,'smoothlen')
        smoothlen = EXPT.FIR.smoothlen;
    else
        smoothlen = input('Enter # tp to weight for exponential smoothing (0 is none, recommend 3 for TR=2): ');
    end
    
    for i = startat:length(d)
    
        fprintf(1,'\n*---------\n%s -- SAVING IN TIMECOURSE_PLOTS SUBDIR\n',d(i).name)
    
        cd(d(i).name)
    
        d2 = dir(name);
        for rr = 1:length(d2)
            fprintf(1,'\t%s\n',d2(rr).name)
            clear cl
            load(d2(rr).name)
    
            if ~isempty(cl)
                plot_dx_hrfs(EXPT,cl,0,1,smoothlen,doindiv);
            end
        
            
        end
        cd ..
    end
    
    
end % do plot part




return