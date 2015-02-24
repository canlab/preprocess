function montage_clusters(ovl,clusters,varargin)
% montage_clusters(ovl,clusters,varargin)
% NOW PLOTS ALL CLUSTERS IN RED!!  2nd...nth clusters are plotted in respective colors, as are overlap
% using image patch method.
%
% by Tor Wager  last edit 12/11/02
%
% varargin (in any order) = 
%   a) additional clusters structures 
%   b) cell array of colors (text format), must be ROW vector {'r' 'g'} etc...
%       if length of color string is longer than number of clusters inputs,
%       additional colors are interpreted as '2-intersection' and 'all-intersection' 
%       colors, in that order.  This overrides single string argument (c, below) for
%       color input 
%   c) single string color argument for overlaps (intersections)
%       plots intersections of ANY TWO clusters right now.
%       also color for plotting points, if entered.
%   d) [n x 3] matrix of points to plot on map 
%   e) text labels for points, must be cell array in COLUMN vector
%   f) single number, 1/0 for whether to plot overlapping coordinates in overlap colors
%       default is 1.
%
% Intersections of 2 colors are magenta, and ALL colors are yellow
% unless otherwise specified

% ----------------------------------------
% * defaults
% ----------------------------------------

myc = {'b' 'r' 'g' 'c' 'm' 'w'};
customcolors = 0;
bcolor = 'm';
acolor = 'y';
XYZmm_both = [];
clindx = 1;
XYZpts = [];
myplab = [];
cl = [];
plotovl = 1;
doverb = 1;

XYZmm = cat(2,clusters.XYZmm);
XYZmm_all = XYZmm;

% ----------------------------------------
% * process input arguments
% ----------------------------------------

if length(varargin) > 0
    for i = 1:length(varargin)
        if isstruct(varargin{i})
            % struct inputs interpreted as clusters variables
            cl{clindx} = varargin{i};
            clXYZmm = cat(2,cl{clindx}.XYZmm);
            clindx = clindx + 1;
            a = XYZmm'; b = clXYZmm';
            XYZmm_both = [XYZmm_both intersect(a,b,'rows')'];
            if ~isempty(XYZmm_both), XYZmm_all = intersect(XYZmm_all',b,'rows')';, else, XYZmm_all = [];, end
            XYZmm = [XYZmm clXYZmm];
           
        elseif isstr(varargin{i})
            % string arguments interpreted as overlap colors
            bcolor = varargin{i};
        elseif iscell(varargin{i})
            % cell array vectors with one row interpreted as color inputs
            if size(varargin{i},1) == 1
                customcolors = 1;
                myc = varargin{i};
            elseif size(varargin{i},2) == 1
                % cell array vectors with one column interpreted as point coordinates to plot
                myplab = varargin{i};
            else
                error('cell array must be row (for colors) or column (for text labels) vector.')
            end
        elseif prod(size(varargin{i})) == 1
            % single integers interpreted as 'do overlap plot' flag, can be 1 or 0
            % default is 1
            plotovl = varargin{i};
        elseif any(size(varargin{i})) == 3
            % any 3-vector matrix interpreted as list of points to plot on figure
            XYZpts = varargin{i};
            XYZmm = [XYZmm XYZpts];
            % then it's a list of points to plot
            % not implemented yet.
        else
            error('Unknown input argument type.')
        end
    end
end

if length(cl) < 1, plotovl = 0;, end

if customcolors
    if length(myc) > length(cl)+2, acolor = myc{length(cl) + 3};, end
    if length(myc) > length(cl)+1, bcolor = myc{length(cl) + 2};, end
end
  
if doverb
    disp([num2str(length(cl) + 1) ' Clusters found.'])
    if plotovl
        disp(['Two cluster overlap: ' bcolor ', found ' num2str(length(XYZmm_both)) ' voxels.'])
        disp(['All cluster overlap: ' acolor ', found ' num2str(length(XYZmm_all)) ' voxels.'])
    else
        disp(['No overlap plotting.'])
    end
end

% ----------------------------------------
% * overlay image
% ----------------------------------------
V = spm_vol(ovl);
oimg = spm_read_vols(V);
V.M = V.mat;

textx = size(oimg,1) - 50;
texty = 6; %size(oimg,2) - 6;

%[array,hdr,h,whichslices,rows,cols,figh] = readim2(ovl,'p');

% how many slices, which ones

XYZ = mm2voxel(XYZmm,V,2)';
whsl = unique(XYZ(3,:));
nsl = length(whsl) + 1;
rc = ceil(sqrt(nsl));
h = [];

f1=figure; colormap gray; set(gcf,'Color','w')
index = 1;
for z = whsl
    
    h(index) = subplot(rc,rc,index); hold on;
    
    if ~isempty(oimg)
        imagesc(oimg(:,:,z)')
        set(gca,'YDir','normal');
        hold on; axis image; axis off
        zmm = voxel2mm([1 1 z]',V.mat);
        text(textx,texty,['z = ' num2str(zmm(3))],'Color','w')
    else
        hold on
    end
    
    index = index + 1;
end
colormap gray

% ----------------------------------------
% * plot first cluster structure (now plot ALL at once!)
% ----------------------------------------
mask = zeros(size(oimg));
mask(sub2ind(size(oimg),XYZ(1,:),XYZ(2,:),XYZ(3,:))) = 1;
mask(:,:,sum(sum(mask)) == 0) = [];
[h,ph,fh] = plotmask4(mask,1,h,f1);




    
% ----------------------------------------
% * plot additional cluster structures
% ----------------------------------------

if length(cl) > 0
    for i = 1:length(cl)
        index = plot_cluster(cl{i},[],rc,V,whsl,myc{i+1},textx,texty,i+1,size(oimg));
        for j = 1:length(cl{i}), 
            if cl{i}(j).numVox == 1, plot_points(cl{i}(j).XYZmm,rc,V,whsl,myc{i+1},myplab);, end
        end

    end
end


% ----------------------------------------
% * plot overlap areas
% ----------------------------------------

if plotovl
    
    if ~isempty(XYZmm_both)
        bcl.XYZmm = XYZmm_both;
        plot_cluster(bcl,[],rc,V,whsl,bcolor,textx,texty,i+1,size(oimg),1);
    end

    if ~isempty(XYZmm_all)
        bcl.XYZmm = XYZmm_all;
        plot_cluster(bcl,[],rc,V,whsl,acolor,textx,texty,i+1,size(oimg),1);
    end

end

% fix bug at end - replot XYZ slice text
%for z = whsl
%    subplot(rc,rc,index); 
%    zmm = voxel2mm([1 1 z]',V.mat);
%    text(textx,texty,['z = ' num2str(zmm(3))],'Color','w')
%    zi(z) = zmm(3);
%end
%zi(zi > 0)

try
    enlarge_axes(gcf)
catch
    disp('Error enlarging axes. Skipping.')
end

return




% ----------------------------------------
%
% * Sub-functions
%
% ----------------------------------------

function index = plot_cluster(clusters,oimg,rc,V,whsl,myc,textx,texty,clind,odims,varargin)
% varargin suppresses end text

XYZmm = cat(2,clusters.XYZmm);
XYZ = mm2voxel(XYZmm,V)';

% surface patch method
% ----------------------------------------------------------------------------------------
vol = voxel2mask(XYZ',odims);
vol = smooth3(vol);


index = 1;
for z = whsl
    
    subplot(rc,rc,index); 
    
    if ~isempty(oimg)
        imagesc(oimg(:,:,z)')
        set(gca,'YDir','normal');
        hold on; axis image; axis off
        zmm = voxel2mm([1 1 z]',V.mat);
        text(textx,texty,['z = ' num2str(zmm(3))],'Color','w')
    else
        hold on
    end
    

    if z>1,
        mvol = vol(:,:,z-1:z); for i = 1:size(mvol,3),myvol(:,:,i) = mvol(:,:,i)';,end
        %FVC = isocaps(vol(:,:,z-1:z)',0,'zmax');
    else 
        mvol = vol(:,:,z:z+1); for i = 1:size(mvol,3),myvol(:,:,i) = mvol(:,:,i)';,end
        %FVC = isocaps(vol(:,:,z:z+1)',0,'zmax');
    end
    FVC = isocaps(myvol,0,'zmax');
    
    try
	patch(FVC,'EdgeColor','none','FaceColor',myc,'FaceAlpha',1)
    catch
	patch(FVC,'EdgeColor','none','FaceColor',myc)
    end
    
    % plot method
    
    %myXYZ = XYZ(1:2,XYZ(3,:) == z);
    %plot(myXYZ(2,:),myXYZ(1,:),[myc 's'],'MarkerFaceColor',myc,'MarkerSize',3)
    
    index = index + 1;
    drawnow
    
end

if length(varargin) == 0
    subplot(rc,rc,index)
    a = pwd; a = a(end-6:end);
    b = num2str(clusters(1).threshold);

    c = num2str(length(clusters));
    text(0,clind-1,[myc ': ' clusters(1).title ' ' a ' u = ' b ', ' c ' clusters'])
    axis off
    axis([0 1 -1 clind])
end



return




function ph = plot_points(XYZmm,rc,V,whsl,myc,myplab)

XYZ = mm2voxel(XYZmm,V,1)';     % suppress unique voxel output
index = 1;
phind = 1;

for z = whsl
    
    subplot(rc,rc,index);
    hold on
    myXYZ = XYZ(:,XYZ(3,:) == z);
    if ~isempty(myplab), myplz = myplab(XYZ(3,:) == z);, end
    
    for i = 1:size(myXYZ,2)
        ph(phind) = plot3(myXYZ(1,i),myXYZ(2,i),100,[myc(1) '.'],'MarkerFaceColor',myc(1),'MarkerSize',8);
        
        if ~isempty(myplab)
            text(myXYZ(1,i),myXYZ(2,i),100,myplz{i},'Color',myc(1))
        end
        
        phind = phind + 1;    
    end
    
    index = index + 1;   
    %view(0,90)
    %plot(0,0,'kd')
end