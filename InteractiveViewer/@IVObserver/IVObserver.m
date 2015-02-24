% IVObserver (InteractiveViewer observer) is a class that waits for the selected voxel in the orthviews to change
%   and then does something with the information, typically plotting something or printing out info. It implements the 
%   Observer pattern.
%
% Usage:
%   obs = IVobserver(fh); % where fh is a function handle for a function that takes in an InteractiveViewer, mm coordinates, and voxel coordinates
%
% E.g.
%   function fun(iv, mmpos, voxpos)
%       disp(mmpos);
%       disp(voxpos);
%   end
%
%   obs = IVObserver(@fun);
%   attach(iv, obs);
%   obs2 = IVObserver(@(iv, mmpos, voxpos) disp(mmpos));
%   attach(iv, obs2);
%
% For more on the Observer pattern, see:
%   Erich Gamma, Richard Helm, Ralph Johnson, and John Vlissides.
%   Design Patterns: Elements of Reusable Object-Oriented Software. Addison-Wesley, 1995.

function obs = IVObserver(fh)
    if(nargin == 0)
        error('First argument must be a function handle or []');
    elseif(isempty(fh))
        obs = class(struct('iv',{},'updateFun',{}), 'IVObserver');
    elseif(isa(fh, 'IVObserver'))
        obs = fh;
    else
        obs = initFields();
        obs = class(obs, 'IVObserver');

        if(~isa(fh, 'function_handle')), error('Argument is not a function handle in %s\n', mfilename()); end

        if(nargin(fh) ~= 3), error('Update function must take in three parameters, the InteractiveViewer, the mmpos and the voxpos, in %s\n', mfilename()); end

        obs.updateFun = fh;
    end
end

function obs = initFields()
    obs.iv = [];
    obs.updateFun = [];
end