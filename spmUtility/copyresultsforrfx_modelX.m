
oh forget it.
myImgDir = [myMainResDir

if ~(exist(
   eval(['!mkdir ' fmriDIR filesep 'RESULTS/model' num2str(model) '/imgs'])
catch
end

eval(['cd fmriDIR filesep 'RESULTS/model' num2str(model) '/imgs'])

disp(['Copying img files for all subjects '])
   
% copy files
% =====================================================================================
%fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'con*');
%fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'spmT*');
%fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'spmF*');

try
fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'ncon*');
%fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'nspmT*');
%fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'nspmF*');
catch
end

try
fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'sncon*');
%fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'snspmT*');
%fixedeffects('copy','/data',drive,['/intext2/RESULTS/model' model '/sub'],ss,[],'snspmF*');
catch 
end

%str = ['!cp ' fmriDIR filesep 'RESULTS' filesep 'model' num2str(model) filesep 'sub*/con* ' fmriDIR filesep 'RESULTS' filesep 'model' num2str(model) filesep 'imgs']
%eval(str)

%str = ['!cp ' fmriDIR filesep 'RESULTS' filesep 'model' num2str(model) filesep 'sub*/spmT* ' fmriDIR filesep 'RESULTS' filesep 'model' num2str(model) filesep 'imgs']
%eval(str)

%str = ['!cp ' fmriDIR filesep 'RESULTS' filesep 'model' num2str(model) filesep 'sub*/spmF* ' fmriDIR filesep 'RESULTS' filesep 'model' num2str(model) filesep 'imgs']
%eval(str)
