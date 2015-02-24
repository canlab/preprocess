function batch_set_origins(dirs, varargin)
%
% Method to set all origins at once. Useful if scnlab_preproc_part2 needs
% to be run on multiple subjs. Afterwards, use
% scnlab_preproc_part2_noorigins instead.
%
% E.g.:
% dirs = filenames('rea*');
% obj = 'structural/T1.img';
% targ = 'structural/T1inplane.img';
% func = 'r*/ravol*img';
% batch_set_origins(dirs, 'func', func, 'obj', obj, 'targ', targ);


obj = 'structural/T1.img';
targ = 'structural/T1inplane.img';
func = 'r*/ravol*img';

%process var args
if length(varargin) > 0
    for i = 1:length(varargin)
        if strcmp(varargin{i},'targ'), targ = varargin{i+1};,end
        if strcmp(varargin{i},'obj'), obj = varargin{i+1};,end
        if strcmp(varargin{i},'func'), func = varargin{i+1};,end
    end
end

currentdir = pwd();
for i=1:length(dirs)
    cd(dirs{i});
    
    % set the origin of the hi-res SPGR
    set_hdr_current_coords(obj);

    % get list of functional files
    p = get_filename2(func);

    % set the origin of the in-plane T1 and adjust all functionals
    set_hdr_current_coords(targ,p);
    
    
    cd(currentdir);
end

return


