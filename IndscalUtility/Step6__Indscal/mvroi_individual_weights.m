function indX = mvroi_individual_weights(DATA)
%indX = mvroi_individual_weights(DATA)
% Requires INDSCAL field with the following fields:
%W = DATA.INDSCAL.W;
%Gs = DATA.INDSCAL.Gs;
%tv = DATA.INDSCAL.taskvector;
%nms = DATA.SPEC.tasknames;

W = DATA.INDSCAL.W;
Gs = DATA.INDSCAL.Gs;
tv = DATA.INDSCAL.taskvector;
nms = DATA.SPEC.tasknames;

colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
while size(Gs,1) > length(colors), colors = [colors colors];,end


% stimulus coordinates (reconstructed) for each subject (3rd dim)
for i = 1:size(W,1), indX(:,:,i) = Gs * diag(W(i,:));, end

figure; hold on;
for i = 1:size(W,1), 
    for j = 1:size(indX,1), 
       % plot(indX(j,1,i),indX(j,2,i),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:));, 
    end 
end

pts = permute(indX,[3 2 1]);    % put subjects x task on the 1st dim, regions on the 3rd       
    
% plot group points
%for j = 1:size(Gs,1), hold on; plot(Gs(j,1),Gs(j,2),'ks','MarkerFaceColor',colors{j}(1),'MarkerSize',8);,end

% plot confidence regions

%for i = 1:size(pts,3)
%    [ax,ci,cen,ub,lb,S,e,lam] = conf_region(pts(:,1:2,i),colors{i});
%end


if ~isempty(tv)
    
    figure('Color','w');
    % separate by task states
    for i = 1:max(tv)
        
        ptstate{i} = pts(find(tv==i),:,:);
        
        h(1) = subplot(1,max(tv)+1,1);hold on;set(gca,'FontSize',14);
        for j = 1:size(pts,3)
            hh(j) = plot(mean(ptstate{i}(:,1,j)),mean(ptstate{i}(:,2,j)),colors{i},'MarkerFaceColor',colors{i}(1));
        end

        if i == max(tv), legend(hh(1:size(pts,3):end),nms);,end
        
        
        h(i+1) = subplot(1,max(tv)+1,i+1); hold on; set(gca,'FontSize',14);
        
        for j = 1:size(pts,3)
            [ax,ci,cen,ub,lb,S,e,lam] = conf_region(ptstate{i}(:,1:2,j),colors{i});
            plot(mean(ptstate{i}(:,1,j)),mean(ptstate{i}(:,2,j)),colors{i},'MarkerFaceColor',colors{i}(1));
        end
        
        title([nms{i}],'FontSize',18);
    end
    
    equalize_axes(h)
    
end


    try, saveas(gcf,'step6_indiv_weights','fig');, 
        saveas(gcf,'step6_indiv_weights','tif');,close,
    catch, disp('cannot save figure');,
    end
    
return