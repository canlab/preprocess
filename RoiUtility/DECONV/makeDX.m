function [DX,sf,hires_sf] = makeDX(Sess,tp,eres,varargin)
% function [DX,sf,hires_sf] = makeDX(Sess,tp,eres,varargin)
%
% varargin 1 = colsOfInterest
%
if nargin > 3
	colsOfInterest = varargin{1};
	if isempty(colsOfInterest), colsOfInterest = 1:size(Sess{1}.sf,2);, end
else
	colsOfInterest = 1:size(Sess{1}.sf,2);
end

% -------------------------------------------------------------------
% * concatenate session stick functions into overall one
% -------------------------------------------------------------------
sf = [];
for i = 1:length(Sess)
	a = cell2mat(Sess{1}.sf(:)');
	sf = [sf;a];
end


% -------------------------------------------------------------------
% * choose columns of interest and convert back to cell array
% -------------------------------------------------------------------

sf = sf(:,colsOfInterest);
try
	sf = mat2cell(sf,size(sf,1),ones(size(sf,2),1));
catch
	whos sf
	error('Can''t convert mat to cell array.')
end

% -------------------------------------------------------------------
% * make deconvolution matrix
% -------------------------------------------------------------------
hires_sf = sf;
[DX,sf] = tor_make_deconv_mtx2(sf,tp,eres);