% VOI = fakeVOI(ts, SPM, VOI_title, varargin)
function VOI = fakeVOI(ts, SPM, VOI_title, varargin)
    mm_location = [0 0 -50]';
    
    cl.timeseries = ts(:);
    cl.shorttitle = VOI_title;
    
    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'mm_location'
                    mm_location = varargin{i+1}(:);
            end
        end
    end
    
    cl.mm_center = mm_location;
    cl.XYZmm = mm_location;
    
    VOI = cltoVOI(cl, SPM, varargin{:});
end