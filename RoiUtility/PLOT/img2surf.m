function [cl] = img2surf(mmdist,mycolors,P,P2,P3,varargin)
% [varargout] = img2surf(mmdist,mycolors,P,P2,P3,varargin)
% by Tor Wager
%
%       mmdist = distance in mm for coloring surface vertices from nearest voxel
%
%       mycolors = list of colors in cell array {[0 0 1] [1 0 0] ...etc.}
%       if number of colors specified is greater than number of clusters 
%       structures entered, n+1 and n+2 colors are overlap of 2 and overlap
%       of all clusters, respectively.
%
%       P, P2, P3: surface mat files containing Faces and Vertices for the
%       brain rendering.  P2 = left medial surface, P3 = R med, P = lateral surf.
%
%       Output argument is cl{:} cell array of cluster structures.

str = [];

if length(varargin) == 0, disp('No masks entered as input.'),return, end
    
for i = 1:length(varargin)
    cl{i} = mask2clusters(varargin{i});
    str = [str ',cl{' num2str(i) '}'];
end

if isempty(P)
    P = 'C:\tor_scripts\3DheadUtility\canonical_brains\surf_scalped_single_subj_T1.mat';
    if ~(exist(P) == 2), P = spm_get(1,'*img','Choose lateral surface mat file');, end
end

if isempty(P2)
    P2 = 'C:\tor_scripts\3DheadUtility\canonical_brains\surf_single_subj_grayL.mat';
    if ~(exist(P) == 2), P = spm_get(1,'*img','Choose L medial surface mat file');, end
end

if isempty(P3)
    P3 = 'C:\tor_scripts\3DheadUtility\canonical_brains\surf_single_subj_grayR.mat';
    if ~(exist(P) == 2), P = spm_get(1,'*img','Choose R medial surface mat file');, end
end

  sstr = ['cluster_surf(P,mmdist,mycolors' str ');'];
  disp(sstr)
  eval(sstr)
  lightfollowview; set(gcf,'Color','w');
  
    sstr = ['cluster_surf(P2,mmdist,mycolors' str ');'];
  disp(sstr)
  eval(sstr)
  view(90,0); lightfollowview; set(gcf,'Color','w');
  
    sstr = ['cluster_surf(P3,mmdist,mycolors' str ');'];
  disp(sstr)
  eval(sstr)
  view(270,0); lightfollowview; set(gcf,'Color','w');
  
 return
  
  