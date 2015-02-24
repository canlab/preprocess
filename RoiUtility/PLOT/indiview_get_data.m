function [coord,con,mm] = indiview_get_data;


global f
global f2
global VOL
global con

% for HTW images
global NLCON
if ~isstruct(NLCON), NLCON = [];, end

% get coordinate
coord = spm_orthviews('Pos');
coord = mm2voxel(coord',VOL);
mm = round(spm_orthviews('Pos')');



    % -------------------------------------------------------------------
    % * load images, if necessary, and get data
    % -------------------------------------------------------------------

    if isstr(con{1}), disp('Loading images (this may take a few minutes).');,end
   
    % if we already have volumes, P is image data; otherwise, names, and
    % vols are loaded
   
    for i = 1:length(con)
        [ts(i),con{i}] = timeseries3(coord,con{i});
        
        dat(:,i) = ts(i).indiv;
    end
 
    titlestr = (['mm: ' num2str(mm) ' vox: ' num2str(coord)]);

    figure(f); cla; barplot_columns(dat,titlestr,[],'nofig');
    
    
    
    
    % -------------------------------------------------------------------
    % * H, T, W view
    % ------------------------------------------------------------------- 
    if ~isempty(NLCON)
        if isstr(NLCON.height{1}), disp('Loading H, T, W images from NLCON (this may take a few minutes).');,end
        disp('H, T, W image view'); disp('-------------------------');
        
        for i = 1:length(NLCON.height)
            [tsh(i),NLCON.height{i}] = timeseries3(coord,NLCON.height{i});
            [tsd(i),NLCON.delay{i}] = timeseries3(coord,NLCON.delay{i});
            [tsw(i),NLCON.width{i}] = timeseries3(coord,NLCON.width{i});
            
            dath(:,i) = tsh(i).indiv;
            datd(:,i) = tsd(i).indiv;
            datw(:,i) = tsw(i).indiv;
        end

        figure(f2); 
        subplot(1,3,1); cla; barplot_columns(dath,'Height',[],'nofig','noind');
        if isfield(NLCON,'names'), set(gca,'XTickLabel',NLCON.names),end
        
        subplot(1,3,2); cla; barplot_columns(datd,'Time to Peak',[],'nofig','noind');
        if isfield(NLCON,'names'), set(gca,'XTickLabel',NLCON.names),end
        
        subplot(1,3,3); cla; barplot_columns(datw,'Width',[],'nofig','noind');
        if isfield(NLCON,'names'), set(gca,'XTickLabel',NLCON.names),end
    end
    
return