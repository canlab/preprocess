function [clusters,cmatx] = simple_correl_clusters(clusters,interf,varargin)
% [clusters,cmatx] = simple_correl_clusters(clusters,interf,[opt] doplot or point labels)
% Tor Wager

doplot = 0;
if nargin > 2, 
	doplot = varargin{1};,
	if iscell(doplot), mylabels = doplot;, doplot = 1;,end
end

    
if nargin > 3, M = varargin{2};, 
elseif isfield(clusters,'M'),M = clusters(1).M;
end

for i = 1:length(clusters)
    
  try
   
    if size(interf,2) > size(interf,1), interf = interf';,end
    
    a = corrcoef(interf(~isnan(interf)),clusters(i).timeseries(~isnan(interf)));
    
    if doplot,
        
        X = [interf ones(length(interf),1)];
	X(isnan(interf),:) = [];

        % transform to z-scores for plotting comparability
        %y = (clusters(i).timeseries - mean(clusters(i).timeseries)) ./ std(clusters(i).timeseries);
 	
	% here, don't transform       
	y = (clusters(i).timeseries);
	y(isnan(interf)) = [];
	b = X \ y;
        
        figure; hold on; grid on;
        plot(interf,clusters(i).timeseries,'ko','MarkerSize',6,'MarkerFaceColor','k')
        x = 0:max(interf);
        plot(x,b(1) * x + b(2),'k-','LineWidth',2)
	if exist('mylabels') == 1
		for j = 1:length(interf)
			text(interf(j),clusters(i).timeseries(j),mylabels{j})
		end
	end
        ylabel('Contrast beta value','FontSize',14)
        xlabel('Behavioral score','FontSize',14)
        set(gcf,'Color','w')
        dr = pwd;
        title([dr(end-7:end) ' cl = ' num2str(i) ' r = ' num2str(a(1,2))])
    end
    
    cmatx(i,1) = a(1,2);
    cmatx(i,2) = clusters(i).numVox;
    cmatx(i,3) = max(clusters(i).Z);
    clusters(i).correl = a(1,2);
    
    % get range of individual voxel correlations
    cv = [];
    for j = 1:size(clusters(i).all_data,2)
        cm = corrcoef(interf,clusters(i).all_data(:,j));
        cv(j) = cm(1,2);
    end
    clusters(i).corr_range = [min(cv) max(cv)];
    
    % not working yet
    %clusters(i).cor_stat = tor_r2z(a(1,2),length(clusters(i).timeseries-1);
    
    % for use with tor_get_spheres2.m, breaking up cluster into spheres
    if isfield(clusters,'center') & exist('M') == 1 & isfield(clusters,'from_cluster')
        cmatx(i,4) = clusters(i).from_cluster;
        centers{i} = (round(voxel2mm(clusters(i).XYZ(:,clusters(i).Z == max(clusters(i).Z)),M)'));
        cmatx(i,5) = centers{i}(1);
        cmatx(i,6) = centers{i}(2);
        cmatx(i,7) = centers{i}(3);
    end
    
  catch
    a = corrcoef(interf,clusters(i).timeseries);
    cmatx(i,1) = a(1,2);
    disp(['Problem with cluster ' num2str(i)])
    cmatx(i,2) = clusters(i).numVox;
     cmatx(i,3) = max(clusters(i).Z);
     
    % for use with tor_get_spheres2.m, breaking up cluster into spheres
    if isfield(clusters,'center') & exist('M') == 1 & isfield(clusters,'from_cluster')
        cmatx(i,4) = clusters(i).from_cluster;
        centers{i} = (round(voxel2mm(clusters(i).XYZ(:,clusters(i).Z == max(clusters(i).Z)),M)'));
        cmatx(i,5) = centers{i}(1);
        cmatx(i,6) = centers{i}(2);
        cmatx(i,7) = centers{i}(3);
    end
    
 end
    
end

%if isfield(clusters,'center') & exist('M') == 1 & isfield(clusters,'from_cluster')
    % sort by which cluster its from
%    cmatx = sortrows(cmatx,4);
%    fprintf(1,'corr\tvoxels\tmaxZ\tfrom_clust\tmax_coords\n')
%    for i = 1:size(cmatx,1)
%        fprintf(1,'%3.2f\t%3.0f\t%3.2f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t\n',cmatx(i,1),cmatx(i,2),cmatx(i,3),cmatx(i,4),cmatx(i,5),cmatx(i,6),cmatx(i,7))
%    end
%else
%    disp(' ')
%    fprintf(1,'corr\tvoxels\tmaxZ\n')
%    for i = 1:size(cmatx,1)
%        fprintf(1,'%3.2f\t%3.0f\t%3.2f\n',cmatx(i,1),cmatx(i,2),cmatx(i,3))
%    end
%end

return