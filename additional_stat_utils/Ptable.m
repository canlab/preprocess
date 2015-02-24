function [P]=Ptable(D,W)


if ~iscell(D)
    for i=1:size(D,1)
        for j=1:size(D,2)
            P(i,j)=sum((D(i,j,:)./sqrt((sum(D(i,:,:),2)))).*(W./(sum(W.*sqrt(sum(D(i,:,:),2))))));
        end
    end
else

end
            