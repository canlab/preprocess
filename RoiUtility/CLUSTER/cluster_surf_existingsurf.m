function cluster_surf(varargin)
% cluster_surf(varargin)
% Tor Wager
% surface plot of clusters on a standard brain
%
% inputs, in any order:
%   clusters structures, as created in tor_extract_rois.m
%   cell array of colors for each cluster: {[1 0 0] [0 1 0] [0 0 1]}
%       if number of colors specified is greater than number of clusters 
%       structures entered, n+1 and n+2 colors are overlap of 2 and overlap
%       of all clusters, respectively.
%   file name of mat file containing brain surface vertices and faces
%       as created with isosurface.
%       OR special string: 'bg' 'hipp' (hcmp,thal,amy)
%   number of mm to plot from surface (mmdeep)
%   optional string: 'colorscale'.  This scales colors by Z-scores of voxels
%       if used, Z scores should be in ROW vector
%   optional string: 'heatmap'.  This is used WITH or instead of 'colorscale', 
%   and specifies
%   that surface colors should be heatmapped to Z-scores, rather than a single
%   color with varying hues.
%
% color [0 1 1] (cyan) is reserved for the overlap color btwn cluster sets.
%
% see also img2surf.m
%
% example:
% P = 'C:\tor_scripts\3DheadUtility\canonical_brains\surf_single_subj_T1_gray.mat';
% cluster_surf(tcl,acl,P,10,{[0 1 0] [1 0 0]},'colorscale','heatmap')

% -------------------------------------------------------------------------
% * set up input arguments
% -------------------------------------------------------------------------
mmdeep = 10;
cscale = 0;
heatm = 0;

clind = 1;
for i = 1:length(varargin)
    
    if isstruct(varargin{i}), cl{clind} = varargin{i};, clind = clind+1;, 
    
    elseif iscell(varargin{i}), mycolors = varargin{i};, 
    
    elseif isstr(varargin{i}), 
        if strcmp(varargin{i},'colorscale'), cscale = 1;
        elseif strcmp(varargin{i},'heatmap'), heatm = 1;
        else, P = varargin{i};, 
        end
    
    elseif ishandle(varargin{i}), p = varargin{i};, % handle for existing surface
        
    else, mmdeep = varargin{i};
        
    end
    
end

if ~exist('mycolors') == 1
    for i = 1:length(cl)
        mycolors{i} = rand(1,3);
    end
end

if ~exist('P') == 1
    %P = spm_get(1,'*mat','Choose brain surface file');
    P = which('surf_single_subj_T1_gray.mat');
end
  
disp('cluster_surf')
disp('___________________________________________')
fprintf(1,'\t%3.0f cluster structures entered\n',length(cl))
disp('  Colors are:')
for i = 1:length(mycolors)
    disp([' ' num2str(mycolors{i})])
end
if length(mycolors) > length(cl)
    disp([' overlap color is ' num2str(mycolors{length(cl)+1})])
    ovlc = ['[' num2str(mycolors{length(cl)+1}) ']'];
else
    ovlc = '[0 1 1]';
end

if length(mycolors) > length(cl)+1 & length(cl) > 2
    disp([' all overlap color is ' num2str(mycolors{length(cl)+2})])
    aovlc = ['[' num2str(mycolors{length(cl)+2}) ']'];
else
    aovlc = '[1 1 1]';
end

disp([' Surface stored in: ' P])
    
% -------------------------------------------------------------------------
% * build xyz list
% -------------------------------------------------------------------------
for i = 1:length(cl)
    xyz{i} = cat(2,cl{i}.XYZmm)';
    
    if cscale | heatm, 
        for j = 1:length(cl{i}), if size(cl{i}(j).Z,1) > size(cl{i}(j).Z,2),cl{i}(j).Z = cl{i}(j).Z'; ,end, end
        Z{i} = cat(2,cl{i}.Z)';,
        
        % order voxels from lowest to highest, so that peak colors appear
        tmp = [xyz{i} Z{i}]; 
        tmp = sortrows(tmp,4);
        xyz{i} = tmp(:,1:3); Z{i} = tmp(:,4);
        
        if ~cscale % if heat map only, set mycolor{1} = [1 1 1]
            mycolors{1} = [1 1 1];
        end
    
    end
    
end

% -------------------------------------------------------------------------
% * build function call
% -------------------------------------------------------------------------

if length(cl) > 2 & exist('aovlc') == 1
    str = ['[c,alld] = getVertexColors(xyz{1},p,mycolors{1},[.5 .5 .5],' num2str(mmdeep) ',''ovlcolor'',' ovlc ',''allcolor'',' aovlc];
else
    str = ['[c,alld] = getVertexColors(xyz{1},p,mycolors{1},[.5 .5 .5],' num2str(mmdeep) ',''ovlcolor'',' ovlc];
end

if heatm, 
    str = [str ',''colorscale'',actcolors{1}'];, 
elseif cscale
    str = [str ',''colorscale'',Z{1}'];, 
end


for i = 2:length(cl)
    str = [str ',''vert'',xyz{' num2str(i) '},mycolors{' num2str(i) '}'];
    if heatm, 
        str = [str ',''colorscale'',actcolors{' num2str(i) '}'];, 
    elseif cscale, 
        str = [str ',''colorscale'',Z{' num2str(i) '}'];, 
    end
end

str = [str ');'];


% ------------------------------------------------------------
% for heatmap option: get actcolors
% -------------------------------------------------------------
if heatm
    actcolors = get_actcolors(Z);    
end

% -------------------------------------------------------------------------
% * run brain surface
% -------------------------------------------------------------------------
%[dtmp,ftmp,etmp]=fileparts(P);

%if strcmp(etmp,'.mat')
    
%eval(['load ' P])
%figure
% p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
%  'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
% lighting gouraud;camlight right
% axis image; myLight = camlight(0,0);set(myLight,'Tag','myLight');
% set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightfollowview
% drawnow


% -------------------------------------------------------------------------
% * run color change
% -------------------------------------------------------------------------
disp([' eval: ' str])
eval(str)


% this for subcortex stuff
elseif strcmp(P,'bg')
    P = which('Tal_Cau.img');
    figure('Color','w');[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
    if findstr(P,'Tal_Cau.img'), str(49:58) = '[.5 .6 .6]'; delete(p(1)); p = p(2);,eval(str),end
    
    P = which('Tal_Glo.img');
    [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
    if findstr(P,'Tal_Glo.img'), str(49:58) = '[.5 .6 .5]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
    
    P = which('Tal_Put.img');
    [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
    if findstr(P,'Tal_Put.img'), str(49:58) = '[.5 .5 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
    
    [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -15],'filename','brain_render_T1');
    set(p(1),'FaceColor',[.6 .4 .3]); colormap copper;material dull;axis off
    h = findobj('Type','Light'); delete(h); [az,el]=view;lightangle(az,el); lightangle(az-180,el-60);
    
elseif strcmp(P,'hipp')   
    P = which('Tal_Hip.img');
    figure('Color','w');[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
    if findstr(P,'Tal_Hip.img'), str(49:58) = '[.5 .6 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
    
    P = which('Tal_Amy.img');
    [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
    if findstr(P,'Tal_Amy.img'), str(49:58) = '[.4 .4 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
    
    P = which('Tal_Tha.img');
    [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
    if findstr(P,'Tal_Tha.img'), str(49:58) = '[.5 .6 .5]'; eval(str); end
    
    [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -20],'filename','brain_render_T1');
    set(p(1),'FaceColor',[.6 .4 .3]); colormap copper;material dull;axis off
    h = findobj('Type','Light'); delete(h); [az,el]=view;lightangle(az,el); lightangle(az-180,el-60);
    
elseif strcmp(P,'amy')
    P = which('rICBM_amygdala.img');
    [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
    if findstr(P,'rICBM_amygdala.img'), str(49:58) = '[.4 .4 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
     
    [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -20],'filename','brain_render_T1');
    set(p(1),'FaceColor',[.6 .4 .3]); colormap copper;material dull;axis off
    h = findobj('Type','Light'); delete(h); [az,el]=view;lightangle(az,el); lightangle(az-180,el-60);
    
else
    error('Must input mat surf file or img file to convert to surf')
end



disp('Finished!')
disp('___________________________________________')

return






function actcolor = get_actcolors(sz)

% ------------------------------------------------------------
% for heatmap option: define color maps - biscale hot/cool
% -------------------------------------------------------------

	% color map - hot
	% --------------------------------------------
	h1 = (0:1/99:1)';
	h2 = ones(size(h1)); 
	h3 = zeros(size(h1));
	h = [h1 h3 h3; h2 h1 h3; h2 h2 h1];
	h(1:75,:) = []; % take only red values to start
	% in new matlab: h = colormap(hot(300));

	% color map - winter
	% --------------------------------------------
	h1 = (0:1/249:1)';
	h2 = (1:-1/(249*2):.5)';
	h3 = zeros(size(h1));
	hc = [h3 h1 h2];
    
    % -------------------------------------------------------------
	% determine overall z-score range
    % -------------------------------------------------------------
    
	zrange = cat(2,sz{:}); 
	tmp = zrange(zrange > 0);
	tmpc = zrange(zrange < 0);

	if ~isempty(tmp)
		zrange = [min(tmp) max(tmp)];
		zh = zrange(1):(zrange(2)-zrange(1))./224:zrange(2);
		zh = round(zh*100);
        
        if isempty(zh), zh = [1 1 0], end   % only one element?
	end

	if ~isempty(tmpc)
		zrangec = [min(tmpc) max(tmpc)];
		zhc = zrangec(1):(zrangec(2)-zrangec(1))./249:zrangec(2);
		zhc = round(zhc*100);
        if isempty(zhc), zhc = [0 0 1], end   % only one element?
	end
    
    % -------------------------------------------------------------
	% loop through sets of input coordinates
    % -------------------------------------------------------------
    
	for i = 1:length(sz)
        
        % -------------------------------------------------------------
	    % find color for each xyz 
        % -------------------------------------------------------------        
		clear h2,clear wh
        myz = sz{i};
        
		for j = 1:length(myz)
			if myz(j) > 0, docool = 0; else, docool = 1;, end

			if docool,
				tmp = find(round(myz(j)*100) == zhc);
				if isempty(tmp), 
					tmp = find((zhc-round(myz(j)*100)).^2 == min((zhc-round(myz(j)*100)).^2));
				end
			else
				tmp = find(round(myz(j)*100) == zh);
				if isempty(tmp), 
					tmp = find((zh-round(myz(j)*100)).^2 == min((zh-round(myz(j)*100)).^2));
				end
			end

			wh(j) = tmp(1);

			if docool
				actcolor{i}(j,:) = hc(wh(j),:);
			else
				actcolor{i}(j,:) = h(wh(j),:);
			end

		end

    end
    
    
    
    
    % -------------------------------------------------------------
	% color scale bar - we must create by hand
    % -------------------------------------------------------------
    
    try
    
	zrange = cat(2,sz{:}); 
	tmp = zrange(zrange > 0);
	tmpc = zrange(zrange < 0);

    	if ~isempty(tmp)
    		figure('Color','w'); hold on;
		zh2 = zh./100;
            
            for i = 2:size(h,1), fill([zh2(i-1) zh2(i-1) zh2(i) zh2(i)],[0 1 1 0],h(i,:),'EdgeColor','none');, end
		    set(gca,'YTickLabel',''); 
            xlabel('Z-score','FontSize',14)
            
    		docolbar = 0;
    	end

    	if ~isempty(tmpc)
    		figure('Color','w'); hold on;
		zh2 = zhc./100;
    		axis([0 .3 zh2(1) zh2(end)]),hold on
		for i = 1:size(h,1), plot([0 1],[zh2(i) zh2(i)],'Color',hc(i,:));, end
		set(gca,'XTickLabel',''); % ylabel('Z-score')
		h3 = get(gcf,'Position');
    		set(gcf,'Position',[h3(1:2) h3(3)*.3 h3(4)*.5])
    		docolbar = 0;
    	end

    catch 
        figure; disp('Cannot make colorbar.  Only one voxel?')
    end
   return
        