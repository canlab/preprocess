function EXPT = nl_contrasts(EXPT)
% function EXPT = nl_contrasts(EXPT)
% Tor Wager
%
% Defines contrasts (if not already entered) in DX.contrasts
% and creates nlcon*img contrast images for each subject.
%
% loop through
%	1 - subject
%	2 - parameter (height, delay, intercept)
%	3 - contrast

% - - - - - - - - - - - - - - - - - - - - - - - - - -
% > set up contrasts and descriptions if necessary
% - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~isfield(EXPT.DX,'contrasts') | ~isfield(EXPT.DX,'connames')
    EXPT.DX = nl_define_contrasts(EXPT.DX);
elseif isempty(EXPT.DX.contrasts) | isempty(EXPT.DX.connames)
    EXPT.DX = nl_define_contrasts(EXPT.DX);
end

paramlist = {'height' 'delay' 'intercept'};

% - - - - - - - - - - - - - - - - - - - - - - - - - -
% > build list of cond img files P{i} for each param
% - - - - - - - - - - - - - - - - - - - - - - - - - -
        
for p = 1:length(paramlist)
    % build list of cond_img files 
    for c = 1:length(EXPT.DX.regsofinterest)
        if c < 10, myz = '000';, elseif c < 100, myz = '00';, elseif c < 1000, myz = '0';, else myz = [];, end
        if c == 1
            P{p}(1,:) = ['cond_' myz num2str(c) '_' paramlist{p} '.img'];     
        else
            P{p} = str2mat(P{p},['cond_' myz num2str(c) '_' paramlist{p} '.img']);
        end
    end
end

EXPT.DX.condP = P;
EXPT.DX.params = paramlist;

% - - - - - - - - - - - - - - - - - - - - - - - - - -
% > output to workspace
% - - - - - - - - - - - - - - - - - - - - - - - - - -

disp(['nl_contrasts.m   contrasts over nonlinear fits'])
EXPT.DX.contrasts

disp('Cond input files (saved in EXPT.DX.condP)')
for p = 1:length(paramlist)
    P{p}
end

% - - - - - - - - - - - - - - - - - - - - - - - - - -
% > loop through sub - param - contrast
% - - - - - - - - - - - - - - - - - - - - - - - - - -
eval(['cd ' EXPT.studydir])
%cd ..

for mysub = 1:length(EXPT.subjects)
    
    eval(['cd ' EXPT.subjects{mysub}])
    
    
    try
        
    for p = 1:length(paramlist);
           
        % - - - - - - - - - - - - - - - - - - - - - - -
        % > for each contrast
        % - - - - - - - - - - - - - - - - - - - - - - -

        for i = 1:size(EXPT.DX.contrasts)

	        if i < 10, myz = '000';, elseif i < 100, myz = '00';, elseif i < 1000, myz = '0';, else myz = [];, end
	
	        Q = ['nlcon_' myz num2str(i) '_' paramlist{p}];
            
	        V(i) = contrast_image(P{p},Q,EXPT.DX.contrasts(i,:));
    
            
            % save output file name in nlconP (input to rfx)
            % EXPT.DX.nlconP{myparam}{mycon};
            if mysub == 1
                which(Q)
                EXPT.DX.nlconP{p}{i} = which([Q '.img']);
            else
                EXPT.DX.nlconP{p}{i} = str2mat(EXPT.DX.nlconP{p}{i},which([Q '.img']));
            end
            
            % more output
            if mysub == 1 & p == 1
                fprintf(1,'contrast %s%3.0f\t%s\n',myz,i,EXPT.DX.connames{i})
            end
                
        end
        
    end
    
    fprintf(1,'%s ',EXPT.subjects{mysub})
    
    catch
        fprintf(1,'%s(problem) ',EXPT.subjects{mysub})
    end
    
    cd ..
end
    
fprintf(1,'\n')  

return
