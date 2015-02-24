clear; close all
load roc

figure('color','white');
set(gca,'box','on')
hold on;
for s=1:size(roc,1);
        data1=[0 squeeze(roc(s,1,:))' 1];
        data2=[0 squeeze(roc(s,2,:))' 1];

    for p=1:4
        point=squeeze(roc(s,:,p));
        plot(point(1),point(2),'ok','MarkerSize',5,'MarkerFaceColor','k');
        %plot(0,0,'ok','MarkerSize',5,'MarkerFaceColor','k');
        %plot(1,1,'ok','MarkerSize',5,'MarkerFaceColor','k');
    end
    line(data1,data2,'Linewidth',2,'color','k');
    Axis([0 1 0 1]);
end