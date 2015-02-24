function [clout,cl] = cluster_graphic_select(cl,varargin)
% function cluster_graphic_select(cl,[clout])
% 
% Click on orthviews to select clusters of interest to save (as output)
% Tor Wager, jan 13 2006
%
% 
% *** You MUST run this with clout as the name of the output variable, and
% cl as the name of the input variable ***
% When you exit, cl will be replaced with clout (the clusters you have
% saved.)
%
% Usage:
% [clout,cl] = cluster_graphic_select(cl);
%
% Initialize orthviews:
% [clout,cl] = cluster_graphic_select;
%
% To save the currently-clicked cluster:
% [clout,cl] = cluster_graphic_select(cl,clout);


%global cl
%global clout

if length(varargin) > 0  
    % run in "pick" mode
    clout = varargin{1};

else 
    % initialize viewer
    clout = []; 
    cluster_orthviews(cl,'unique');
    set(gcf,'WindowButtonUpFcn','[clout,cl] = cluster_graphic_select(cl,clout);')
end




pos = spm_orthviews('Pos')';


% check to see if we're in a cluster
wh = 0; 
centers = cat(1,cl.mm_center);

% find closest cluster, based on center
d = distance_euclid(pos,centers); wh = find(d == min(d)); wh = wh(1);

% only accept if cursor is w/i 2 mm of a voxel in the cluster
d = distance_euclid(pos,cl(wh).XYZmm'); d = min(d); if d > 2, wh = 0; end

if wh
    %cluster_table(cl(wh));
    fprintf(1,'Cl. %3.0f, Voxels: %3.0f, Coords: %3.0f, %3.0f, %3.0f\n',wh,cl(wh).numVox,cl(wh).mm_center(1), cl(wh).mm_center(2),cl(wh).mm_center(3));
    
    saveit = input('Save this cluster? (y/n/x to exit) ','s');
    
    switch saveit,
        case 'y'
            clout = [clout cl(wh)]; 
            cluster_orthviews(cl(wh), {[1 0 0]}, 'add')
            
        case 'n'
        case 'x'
            set(gcf,'WindowButtonUpFcn','');
            cl = clout;
            
    end
end

return

