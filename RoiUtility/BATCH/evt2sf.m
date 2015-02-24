function sf = evt2sf(evt,varargin)
% function sf = evt2sf(evt,[opt] sesslength in TRs,[opt] timebins)
%
% This function takes event times (cell array)
% Where event times are in TRs
% And the first scan is coded as 0
%
% And transforms them into stick (delta) functions
%
% If 1st optional argument is session length in TRs, sf vector
% will equal length sesslength * time bins
%
% Time resolution is TR / 16 ms.
% Unless you input time bins as 2nd argument
% 
% Tor Wager, 10/20/01

timebins = 16;
if nargin > 1, sesslen = varargin{1};, end
if nargin > 2, timebins = varargin{2};, end

if exist('sesslen') == 1, 
    sesslen = sesslen .* timebins;, 
    insl = zeros(sesslen,1);
end

for i = 1:length(evt)
    
    a = evt{i};
    
    % express in time bins
    a = a .* timebins;
    
    % round to nearest time bin
    a = round(a);
    
    % add 1 because time 0 is 1st time bin
    a = a + 1;
    
    % define length of sf output vector
    sl = zeros(max(a),1);
    if exist('sesslen') == 1,
        if length(sl) > length(insl),
	    disp(['Last TR with event is ' num2str(length(sl))])
             error(['Input vector ' num2str(i) ' contains times that happen after session is over!'])
        else
            sl = insl;
        end
    end
    
    sf{i} = sl;
    
    sf{i}(a) = 1;
    
end

return



