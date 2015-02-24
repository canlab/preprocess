function tor_spm2_normalize(subj,Template)
% function tor_spm2_normalize(subj,Template)
%
% subj(i).P = file to normalize
% subj(i).PP = images to apply to
%
% M = smooth_and_mask(P,0,-Inf);    %mask image
%subj(1).P = Q;
%tor_spm2_normalize(subj,M)

if isempty(Template)
    Template = which('scalped_single_subj_T1.img');
end

if ~isfield(subj,'PP')
    subj(1).PP = [];
end

donorm(subj,Template)

return




function donorm(subj,Template)

global defaults
defs = defaults.normalise;



% -----------------------------------------------------------------
% From SPM2
% -----------------------------------------------------------------

for i = 1:length(subj)
    
    subj(i).objmask = '';
    subj(i).matname = [spm_str_manip(subj(i).P,'sd') '_sn.mat'];

    spm_normalise(Template, subj(i).P, subj(i).matname,...
			defs.estimate.weight, subj(i).objmask,defs.estimate);
end
   
a1 = 3;
if a1 == 2 | a1 == 3,
	for i=1:length(subj),
            
        if isempty(subj(i).PP)
            subj(i).PP = subj(i).P;
        else
            subj(i).PP = str2mat(subj(i).P,subj(i).PP);
		    %spm('FigName',['Normalising (write) subj ' num2str(i)],...
		    %	Finter,CmdLine);
        end
		    
            spm_write_sn(subj(i).PP,subj(i).matname,defs.write);
	end;
end;

%fprintf('\n\n');
%spm('FigName','Normalise: done',Finter,CmdLine);
disp(['Normalize subject ' num2str(i) ': Done.'])

return

