function [O] = tor_build_sf(e,varargin);
% function [O] = tor_build_sf(e,O [opt]);
%
% special function for building deconv matrix
% and getting event onset times
%
% e is cell array of all lists from excel
% each list has one column per condition
% each cell contains data from one session
%
% evts contains the events of interest only, to make deconvolution mtx
% allevts contains the same, but without modeling each session separately.
% DX is the deconv matrix for evts
% oDX is the deconv matrix for allevts
%
% Tor Wager 11/26/01

evts = [];
allevts = [];
sfres = 16;		% number of time bins per TR for stick function - default of 16

% ----------------------------------------------------------------------------------
% Process function inputs
% ----------------------------------------------------------------------------------

if nargin > 1, 
	O = varargin{1};
	N = fieldnames(O);
end

nsess = length(e);
disp(['Detected ' num2str(nsess) ' sessions in excel text input.']);

if ~any(strcmp(N,'nscans')),
	nscans = input(' Enter number of scans per session: ');
	O.nscans = nscans;
else nscans = O.nscans;
end

if ~any(strcmp(N,'whichCols')),
	whichCols = input(['excel text input contains ' num2str(size(e{1},2)) ' conditions. Enter nums of interest:']);
	O.whichCols = whichCols;
else	
	whichCols = O.whichCols;
end

if ~any(strcmp(N,'dxTRs')),
	dxTRs = input(' Enter number of TRs in deconvolution matrix: ');
	O.dxTRs = dxTRs;
else
	dxTRs = O.dxTRs;
end


% ----------------------------------------------------------------------------------
% Make event times
% ----------------------------------------------------------------------------------

for i = 1:length(e)
    
   % add session length (360 s = 360 TRs) to events lists
   e{i} = e{i} + (nscans * (i-1));
   
   myEvts = excel2evt(e{i});
   
   % evts contains the events of interest only, to make deconvolution mtx
   % in a scan-specific way
   evts = [evts myEvts(whichCols)];
   
   % allevts contains the same, but without modeling each session separately.
   allevts = [allevts myEvts(whichCols)];
   

   
end

% concatenate allevts into one big list including all sessions
ae{1} = cell2mat(allevts');
allevts = ae;

O.e = e;
O.evts = evts;
O.allevts = allevts;

% ----------------------------------------------------------------------------------
% make stick functions
% ----------------------------------------------------------------------------------

disp('Making stick functions')
try
    evts_sf = evt2sf(evts,nscans*nsess);
    evts_totalsf = evt2sf(allevts,nscans*nsess);
catch
    disp('Check evts  - longer than session?')
	keyboard
    return
end

O.evts_sf = evts_sf;

% ----------------------------------------------------------------------------------
% make deconv matrix
% ----------------------------------------------------------------------------------
% make deconv matrix in 1 TR resolution
disp('Making deconvolution matrix')
try
    [O.DX,shortSf] = tor_make_deconv_mtx2(evts_sf,dxTRs,sfres);
    % O.DX(:,end+1) = 1;
    % now done in function
    O.shortSf = shortSf;
catch
    disp('Problem making DX - check memory?')
    return
end


% make deconv matrix in 1 sec (1 TR) resolution
disp('Making overall deconv matrix')
try
    [O.oDX,O.evts_osf] = tor_make_deconv_mtx2(evts_totalsf,dxTRs,sfres);
    % O.oDX(:,end+1) = 1;
catch
    disp('Problem making oDX - check memory?')
    return
end


return



   