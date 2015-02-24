function montage_clusters(ovl,clusters,varargin)
% montage_clusters(ovl,clusters,varargin)
% 
% by Tor Wager
% varargin = 
%   a) additional clusters structures {'r' 'g'} etc...
%   b) cell array of colors (text format), must be ROW vector
%   c) single string color argument for overlaps (intersections)
%       plots intersections of ANY TWO clusters right now.
%       also color for plotting points, if entered.
%   d) [n x 3] matrix of points to plot on map 
%   e) text labels for points, must be cell array in COLUMN vector
%

myc = {'b' 'r' 'y' 'g' 'm' 'c'};
bcolor = 'w';
acolor = 'w';
XYZmm_both = [];
XYZmm_all = [];
clindx = 1;
XYZpts = [];
myplab = [];
cl = [];

XYZmm = cat(2,clusters.XYZmm);

if length(varargin) > 0
    for i = 1:length(varargin)
        if isstruct(varargin{i})
            cl{clindx} = varargin{i};
            clXYZmm = cat(2,cl{clindx}.XYZmm);
            clindx = clindx + 1;
            a = XYZmm'; b = clXYZmm';
            XYZmm_both = [XYZmm_both intersect(a,b,'rows')'];
            if ~isempty(XYZmm_both), XYZmm_all = [XYZmm_all intersect(XYZmm_both',b,'rows')'];, end
            XYZmm = [XYZmm clXYZmm];
           
        elseif isstr(varargin{i})
            bcolor = varargin{i};
        elseif iscell(varargin{i})
            if size(varargin{i},1) == 1
                myc = varargin{i};
            elseif size(varargin{i},2) == 1
                myplab = varargin{i};
            else
                error('cell array must be row (for colors) or column (for text labels) vector.')
            end
        elseif length(size(varargin{i})) == 2
            XYZpts = varargin{i};
            XYZmm = [XYZmm XYZpts];
            % then it's a list of points to plot
            % not implemented yet.
        else
            error('Unknown input argument type.')
        end
    end
end

% ----------------------------------------
% * overlay image
% ----------------------------------------
V = spm_vol(ovl);
oimg = spm_read_vols(V);
V.M = V.mat;

textx = size(oimg,2) - 50;
texty = size(oimg,1) - 6;

%[array,hdr,h,whichslices,rows,cols,figh] = readim2(ovl,'p');

% how many slices, which ones

XYZ = mm2voxel(XYZmm,V)';
whsl = unique(XYZ(3,:));
nsl = length(whsl) + 1;
rc = ceil(sqrt(nsl));
h = [];


figure; colormap gray; set(gcf,'Color','w')


% plot first cluster structure

index = plot_cluster(clusters,oimg,rc,V,whsl,myc{1},textx,texty,1,size(oimg));
if ~isempty(XYZpts), ph = plot_points(XYZpts,rc,V,whsl,bcolor,myplab);, end

% plot additional cluster structures

if length(cl) > 0
    for i = 1:length(cl)
        index = plot_cluster(cl{i},[],rc,V,whsl,myc{i+1},textx,texty,i+1,size(oimg));
    end
end

if ~isempty(XYZmm_both)
    bcl.XYZmm = XYZmm_both;
    plot_cluster(bcl,[],rc,V,whsl,bcolor,textx,texty,i+1,size(oimg),1);
end

if ~isempty(XYZmm_all)
    bcl.XYZmm = XYZmm_all;
    plot_cluster(bcl,[],rc,V,whsl,acolor,textx,texty,i+1,size(oimg),1);
end

% fix bug at end - replot XYZ slice text
%for z = whsl
%    subplot(rc,rc,index); 
%    zmm = voxel2mm([1 1 z]',V.mat);
%    text(textx,texty,['z = ' num2str(zmm(3))],'Color','w')
%    zi(z) = zmm(3);
%end
%zi(zi > 0)

return



% sub-functions
%

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
        set(gca,'YDir','reverse');
        imagesc(oimg(:,:,z))
        hold on; axis image; axis off
        zmm = voxel2mm([1 1 z]',V.mat);
        text(textx,texty,['z = ' num2str(zmm(3))],'Color','w')
    else
        hold on
    end
    

    if z>1,FVC = isocaps(vol(:,:,z-1:z),0,'zmax');
    else FVC = isocaps(vol(:,:,z:z+1),0,'zmax');
    end
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
    myplz = myplab(XYZ(3,:) == z);
    
    for i = 1:size(myXYZ,2)
        ph(phind) = plot3(myXYZ(2,i),myXYZ(1,i),100,[myc(1) '.'],'MarkerFaceColor',myc(1),'MarkerSize',8);
        
        if ~isempty(myplab)
            text(myXYZ(2,i),myXYZ(1,i),100,myplz{i},'Color',myc(1))
        end
        
        phind = phind + 1;    
    end
    
    index = index + 1;   
    view(0,90)
    plot(0,0,'kd')
end