function [bIsEmpty,clusters] = tor_talairach_clusters(clusters)
% [bIsEmpty,clusters] = tor_talairach_clusters(clusters)
%
% Modified by Tor Wager from the original:
%function bIsEmpty = ihb_GetClusterSet
%----------------------------------------------------------------------------------
% FORMAT bIsEmpty = ihb_GetClusterSet
%----------------------------------------------------------------------------------
% Function to select set of cluster (SPM.mat) and select one cluster
% from this set
% Return value == 1 if cluster set is nonempty and one cluster selected
%                 0 otherwise
%----------------------------------------------------------------------------------
%   10.04.01    Sergey Pakhomov
%==================================================================================
%----------------------------------------------------------------------------------
% Get clusters from SPM.mat
%----------------------------------------------------------------------------------
bIsEmpty = 1;
hFigMain = findobj('Type', 'figure', 'Tag', 'ihb_TalSpaceMain_fig');
ihb_HideShowFigureAll('off');

% Tor changed this: clusters are input now.
% clusters = ihb_GetClusters;

% but we still need to transform to talairach space, so:
%====================================================================
% Talariach volume related data
%====================================================================
tclusters = [];
for i = 1:length(clusters)
    cl = clusters(i);
    cl.isSpmCluster = 0;
    try,cl.hThreshold = cl.threshold;,catch, cl.hThreshold = NaN;, end
    cl.pVoxelLev = NaN;
    cl.pClustLev = NaN;
    clOut = ihb_UpdateClusterTalVoxSize(cl, 10, 10);
    try,nTotalProcessed = nTotalProcessed + clOut.numVox; ,catch,nTotalProcessed = NaN;,end 
    tclusters = [tclusters, clOut];
end
clusters = tclusters;



if isempty(clusters)
    bIsEmpty = 0;
    h = msgbox('No voxels above threshold', 'Warning!', 'warn', 'modal');
    waitfor(h);
    ihb_HideShowFigureAll('on');
    return;
end
setappdata(hFigMain, 'clusters', clusters);
%----------------------------------------------------------------------------------
% Select one of the clusters to view
%----------------------------------------------------------------------------------
[cl, indClusterToView] = ihb_SelectClusters(clusters, 'single');

if isempty(cl)
    bIsEmpty = 0;
    h = msgbox('No clusters selected','Warning','warn', 'modal');  
    waitfor(h);
    ihb_HideShowFigureAll('on'); 
    return; 
end
%----------------------------------------------------------------------------------
% Set clusters UI combo box
%----------------------------------------------------------------------------------
hPopUp = findobj('Type', 'uicontrol', 'Style', 'popupmenu', 'Tag', 'ihb_ClusterPopUp');
set(hPopUp, 'String', {clusters(:).name}, 'Value', indClusterToView);
%----------------------------------------------------------------------------------
% Load 
%----------------------------------------------------------------------------------
ihb_PrepareClusterToView;
%----------------------------------------------------------------------------------
% Make main figure visible
%----------------------------------------------------------------------------------
ihb_HideShowFigureAll('on');
