function Pw = ind_model_nlfit(P,hrf,DX,varargin)
% Pw = ind_model_nlfit(P,hrf,DX,[opt] vb,[opt] mask)
% single-subject shape-free beta extraction and nonlinear model fitting
% fits height and delay parameters with canonical hrf function, 
% and writes height and delay output imgs for the subject
%
% P     = image file names (str matrix) - should be dx_beta*imgs
% hrf   = impulse-response function to fit using lsqcurvefit - doesn't work right now.  standard hrf.
% DX    = structure of deconvolution params, as created in tor_build_deconv_design_ui
% vb    = [optional] verbose output level: 0 none, 1 some, 2 lots
%
% dims  = dimensions of image files in data
% cols  = no. of columns of DX, also no. of beta images to write
% betas = 4-D array of x,y,z,column
%
% Pw    = string matrix of output file names

diary on

vb = 2; if length(varargin) > 0, vb = varargin{1};, end
if vb > 0, t1 = clock;, fprintf(1,'\n\t\tSubject setup...'),end

V = spm_vol(P(1,:));

global dims
global cols

dims = V.dim(1:3);
cols = 3;   % how many parameters to fit for each condition
colnames = {'height' 'delay' 'intercept'};
for m = 1:length(DX.regsofinterest)
    betas{m} = NaN * zeros([dims cols]);
end

mask = ones(dims); 
if length(varargin) > 1, 
    mask = varargin{2};, 
    if isstr(mask), Vm = spm_vol(mask);, mask = spm_read_vols(Vm);,end
end

if vb > 0, fprintf(1,'\n\t\tFinished in %3.0f s',etime(clock,t1)),end

figure; 

% -------------------------------------------------------------------
% * for each slice...
% -------------------------------------------------------------------

for slicei = 1:dims(3)
    
    if vb == 2, t1 = clock;, fprintf(1,'\n\t\tSlice %3.0f ',slicei),end
    
    slicebetas =  process_slice(slicei,P,DX,mask(:,:,slicei));
    
    % store slice betas for each condition in 4-D matrices
    % x, y, z, beta = 4D.  Matrix indicates condition number
    for m = 1:length(DX.regsofinterest)
        betas{m}(:,:,slicei,:) = slicebetas{m};
    end
        
    if vb == 2, t1 = clock;, fprintf(1,'%6.0f s',etime(clock,t1)),end
end

    
% -------------------------------------------------------------------
% * write beta images
% -------------------------------------------------------------------    
 
if vb > 0, t1 = clock;, fprintf(1,'\n\t\tWriting %3.0f beta images.',cols * length(betas)),end

for bi = 1:length(betas)
    Pwb = write_beta(bi,V,betas{bi},colnames);
    if bi == 1, Pw = Pwb;, else, Pw = str2mat(Pw,Pwb);, end
end
    
if vb > 0, fprintf(1,'\t%3.0f s to dir: %s ',etime(clock,t1),fileparts(Pw(1,:))),end

return
    
    
    
    
    
    
    
    
    
function betas = process_slice (slicei,P,DX,varargin)
% ind_model_fit, P should be image names for raw images
% ind_model_nlfit, P should be dx_beta*imgs

mask = []; 
if length(varargin) > 0, 
    mask = varargin{1};, 
    subplot 221;
    imagesc(mask); axis image; axis off; colormap gray,drawnow
    title('mask')
end
 

% -------------------------------------------------------------------
% * nonlinear fitting options - hardcoded
% -------------------------------------------------------------------
myoptions = optimset('LevenbergMarquardt','on');
startestimates = [1 0 0];	% height, onset, intercept (baseline h)
xdata = zeros(DX.numframes,1); xdata(1) = 1;
%xdata(end+1) = DX.TR;

lb = [-Inf 3 -Inf];
ub = [Inf 8 Inf];

%lb = [-Inf 0 -Inf];
%ub = [Inf 3 Inf];


global dims
global cols
global TR
TR = DX.TR;

% -------------------------------------------------------------------
% * load the slice
% -------------------------------------------------------------------
O.z = slicei;
sl = timeseries2('slice',P,O);
if ~isempty(mask), sl(:,:,1) = sl(:,:,1) .* mask;, end

subplot 222;
imagesc(sl(:,:,1)); axis image; axis off; colormap gray;title('slice tp 1'),drawnow

% -------------------------------------------------------------------
% * fit the model to each nonzero voxel
% -------------------------------------------------------------------
for m = 1:length(DX.regsofinterest), betas{m} = NaN * zeros([dims(1:2) 3]);,end
wvox = find(sl(:,:,1) ~= 0 & ~isnan(sl(:,:,1)));
[i,j] = ind2sub(size(sl(:,:,1)),wvox);

fprintf(1,'\t %3.0f voxels ',length(i))
diary off

for k = 1:length(i)
    
    b = squeeze(sl(i(k),j(k),:));
    % scale betas to resting baseline and break into conditions
    % ---------------------------------------------------------------
    bc = beta2conditions(b,DX);
    
    if mod(k,5) == 0
        disp(['------------------------------------------------------'])
        disp(' ')
        mydir = pwd;
        disp(['...' mydir(end-10:end) '  Voxel ' num2str(k) ' of ' num2str(length(i))])
        disp(' ')
        disp(['------------------------------------------------------'])
    end
    
    if ~(any(isnan(cell2mat(bc(:)))))
    	for m = 1:length(bc)
            
            y = bc{m};
            
            if any(y ~= 0)
         
                % custom starting estimates

                t = abs(y - mean(y)); s1 = y(t == max(t)) - mean(y); s1 = s1(end);
                s2 = find(t == max(t)); s2 = s2(end);
                s3 = mean(y);
            
        	    [beta,resnorm,residual,exitflag,output,lambda,J]= ...
        	    lsqcurvefit('nlhrf3',[s1 s2 s3],xdata,y,lb,ub,myoptions);
    
    		    if exitflag < 0, beta = beta * Inf / Inf;, end
        	    betas{m}(i(k),j(k),:) = beta;
                
            end
    	end
    else
	for m = 1:length(bc)
		betas{m}(i(k),j(k),:) = NaN;
	end
    end

end

subplot 223;
imagesc(betas{1}(:,:,1)); axis image; axis off; colormap jet;title('cond 1 height'),drawnow

subplot 224;
imagesc(betas{2}(:,:,2)); axis image; axis off; colormap jet;title('cond 1 delay'),drawnow

diary on
return



function Pw = write_beta(bi,V,betas,colnames)

if bi < 10, myz = '000';, elseif bi < 100, myz = '00';, else myz = '000';,end

for nm = 1:length(colnames)
    V.fname = ['cond_' myz num2str(bi) '_' colnames{nm} '.img'];
    V.descrip = [colnames{nm} ' of ' num2str(length(colnames)) ' nonlinear fit parameter estimates for cond ' num2str(bi)];

    spm_write_vol(V,betas(:,:,:,nm));
    Pw = which(V.fname);
end

return
