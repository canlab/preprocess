% SET Set InteractiveViewer properties

function set(iv, varargin)
    if(mod(length(varargin),2) ~= 0)
        error('There appears to be an unmatched property/value pair in %s', mfilename);
    end
    iv.set(varargin{:});
end