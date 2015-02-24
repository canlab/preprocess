function [clusters,HP] = analyze_cluster_rois(clusters,xX,HP)
% function [clusters,HP] = analyze_cluster_rois(clusters,xX,HP filter cutoff [in s])
%
% inputs:
%
% clusters are as defined by Talairach Space Utility add-in
% vector of structures, each element is an spm cluster structure.
%
% xX is design matrix struct from SPM
% HP filter is high-pass filter cutoff in s
%
% outputs:
%	clusters(i).b = betas from design matrix.
%	clusters(i).adjustedy = filtered and adjusted data.
%	clusters(i).HPlength = HP filter length used. 
%
% Right now no low-pass filtering.  Use O.doLP = 1 to do this.
%
% by Tor Wager, 10/19/01

if nargin < 3
	HP = input('Enter high-pass filter cutoff in s, or return for none: ');
end

O.TR = xX.RT;
O.scanadjust = 1;
O.nscans = length(xX.K);
if ~isempty(HP),O.HP = HP;,O.doHP = 1;,end
O.percent = 1;
O.filtertype = 'spm';

fprintf(1,['Cluster '])

for i = 1:length(clusters)

	fprintf(1,[num2str(i) ' '])
	% disp('_____________________________________________________') 


	O.y = clusters(i).timeseries;
	y = filterAdjust(O);

	clusters(i).b = pinv(xX.X) * y;
	clusters(i).adjustedy = y;
	clusters(i).HPlength = O.HP;

end

return