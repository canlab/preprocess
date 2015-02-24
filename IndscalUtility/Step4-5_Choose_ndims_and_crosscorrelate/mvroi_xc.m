function [DATA,xc,xl,resid] = mvroi_xc(DATA,varargin)
% [DATA,xc,xl,residuals] = mvroi_xc(DATA,varargin)
%
% compute cross-correlations among a set of regions for each subject
% given c.data, a time x region x subject data matrix
%
% var args:
% 'doresid', ['yes'] or 'no'                 = calculate residuals from model
% 'doglobal, 'none' 'scale' or ['regression']
%           = way of removing global mean across timeseries
% residuals: a matrix of residuals following the removal of the FIR
% estimates.  same size as the data matrix.
%
% Outputs:
%   DATA.DATA.corr     cell vector of correlation matrices for each subject
%   DATA.DATA.lat      cell vector of latency matrices for each subject
%   DATA.DATA.xc       3-D matrix, reg x reg x state within subject
%                       suitable for input to decomposition 
    
numsub=length(DATA.DATA.resids);    %number of subjects


% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
shift_by = 4;
robustflag = 1;              % IRLS, robust regression
betaflag = 0;                % return beta weights rather than correlations

doresid = 'yes';
doglobal = 'regression';     % 'none' 'scale' or 'regression'


% ------------------------------------------------------------
% Check for required vars and process inputs
% ------------------------------------------------------------

if isfield(DATA.SPEC,'robustflag'), robustflag = DATA.SPEC.robustflag;,end               % robust IRLS regression
if isfield(DATA.SPEC,'betaflag'), betaflag = DATA.SPEC.betaflag;,end 

if isfield(DATA.SPEC,'shiftby'), shift_by = DATA.SPEC.shiftby;,end
DATA.SPEC.shiftby = shift_by;

if ~isfield(DATA.SPEC,'hanningwidth'), DATA.SPEC.hanningwidth = [];,end

if ~isfield(DATA.SPEC,'states'), 
    fprintf(1,'\n*** no states found; assuming single task condition.');,
    DATA.SPEC.states{1} = ones(size(DATA.DATA.dat{1},1),1);,
    for s = 1:numsub, DATA.SPEC.states{s} = DATA.SPEC.states{1};,end
end
    

% variable input arguments
% these should be pairs in the form 'field name', value
% goes through all var args and enters field names as vars in the workspace

% doresid and doglobal are inputs
for i = 1:2:length(varargin)
    
    N = varargin{i};    % get name of variable
    eval([N ' = varargin{i+1};']);
end

% save defaults; these are replaced by input arguments
% don't do it, because we have turned these off locally, but on globally,
% so they're done in other parts of the program, not in this fucntion.
%DATA.SPEC.doresid = doresid;
%DATA.SPEC.doglobal = doglobal;


% ------------------------------------------------------------
% Main function
% ------------------------------------------------------------


fprintf(1,'Cross correlation stage \n\t')

fprintf(1,'Number of task states: %3.0f\n',max(cat(1,DATA.SPEC.states{:})));

resp = {'Off' 'On'}; resp2 = {'Correlations' 'Regression slopes (betas)'};

fprintf(1,'Robust IRLS is: %s, Saving in xc and corr fields for stats: %s\n',resp{robustflag+1},resp2{betaflag+1});

if isempty(DATA.SPEC.hanningwidth)
    fprintf(1,'\tModel removal at this stage = %s\n\tGlobal removal at this stage = %s\n\tNo Hanning win\n\tCross-correlating (lag is %3.0f)\n\t Subject...',doresid,doglobal,shift_by);
else
    fprintf(1,'\tModel removal at this stage = %s\n\tGlobal removal at this stage = %s\n\tHanning (width = %3.0f)\n\tCross-correlating (lag is %3.0f)\n\t Subject...',doresid,doglobal,DATA.SPEC.hanningwidth,shift_by);
end

for i=1:numsub            %size(DATA.dat,3);
    
    %m= c.data(:,:,i);  
    m = DATA.DATA.resids{i};
    %m = DATA.dat{i};
    
    % ------------------------------------------------------------
    % Model fitting and residuals
    % ------------------------------------------------------------

    if strcmp(lower(doresid),'yes')
        b = pinv(DATA.DX{i}) * m;       % betas
        f = DATA.DX{i} * b;             % fit
        r = m - f;             % residuals
        
    elseif strcmp(lower(doresid),'no')
        r = m;
    else
        error('Choose yes or no for doresid in scan_xc')
    end
       
    resid{i} = r;
    
    % ------------------------------------------------------------
    % Global mean across regions -- scaling or regress out?
    % ------------------------------------------------------------
   
    if strcmp(lower(doglobal),'none')
        % do nothing
        
    elseif strcmp(lower(doglobal),'scale')
        % z-score columns and center rows (double-center)
        r=scale(r);r=scale(r')';
    
        
    elseif strcmp(lower(doglobal),'regression')
        % remove mean across regions by regression
        gm = mean(r,2); gm(:,end+1) = 1;
        
        gb = pinv(gm) * r;
        gf = gm * gb;
        r = r - gf;             % residuals for this subject
        
    else
        error('Choose yes or no for doresid in scan_xc')
    end

       
    fprintf(1,'%3.0f ',i);
    
    % DATA.SPEC.hanningwidth    contains width of hanning filter for
    % blocks, or empty for no Hanning window
    %
    % DATA.SPEC.states{s}       vector for each subject s of what the block
    % state is.  Correlations are computed within states.
    
    % loop through state types (coded by integers)
    for n=1:max(DATA.SPEC.states{i})    
        
        % residuals for this subj in this block state
        wh = find(DATA.SPEC.states{i}==n);
        wh(wh > length(r)) = [];    % in case state is too long for this subject
        cr = r(wh,:);  
        
            
        %%% make windowed residuals %%%%%%%%%%%
        % only if DATA.SPEC.hanningwidth is not empty (should contain width in
        % elements)
      
        if n == 1,
            [cr,h] = hanning_taper(DATA.SPEC.hanningwidth,cr);
        else
            [cr] = hanning_taper(h,cr);
        end
     
        %%%% cross-correlate %%%%%%%%%%%%%%
        warning off     % for robustfit iteration limit
            
        for j = 1:size(m,2)-1
            for k = j + 1 : size(m,2)
                % now can handle shift by 0 in shift_correl
                %if shift_by>0;
                
                    [sval,myc,nxl(j,k),nxc(j,k)] = shift_correl(cr(:,j),cr(:,k),shift_by,robustflag,betaflag);
                    
                %else
                %    rr=corrcoef(cr(:,j),cr(:,k));    %if there's no shift by
                %    nxc(j,k)=rr(1,2);
                %    nxl=[];
                %end
            end
        end   
        
        warning on
        
        nxc(end+1,:) = 0;   
        nxc=nxc+nxc'+eye(size(nxc,1));
		
        xc(n,:,:,i)=nxc;
        DATA.DATA.corr{i}(n,:,:) = nxc;
        clear nxc;
        
        nxl(end+1,:) = 0;   
        nxl=nxl+nxl';
		
        xl(n,:,:,i)=nxl;
        DATA.DATA.lat{i}(n,:,:) = nxl;
        clear nxl;
        
    end         % state loop
 

    
end             % subject loop


	%DATA.corf{s}=corf;
	%DATA.corb{s}=corb;
    fprintf(1,'\n')
    
% -----------------------------------------------------------------------
% get cross-correlations into the right format for indscal
% should be regions x regions x (states within subjects) 
% -----------------------------------------------------------------------
numstates = size(xc,1);
ind = 1; clear xc2

% xc2 has task states nested within subjects, along 3rd dim
for i = 1:length(DATA.DATA.corr), xc2(ind:ind+numstates-1,:,:) = DATA.DATA.corr{i};, ind = ind + numstates;, end
xc2 = permute(xc2,[2 3 1]);
xc = xc2;
DATA.DATA.xc = xc2;
DATA.DATA.xc_descrip = ['Cross-correls, lag ' num2str(DATA.SPEC.shiftby) ', format for input to decomposition.'];



return




