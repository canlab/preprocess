function out=testsigW(W,cmat);

ncond=max(cmat);
nsub=size(W,1)/ncond;
r=size(W,2);

%if r~=2;
%    error('not yet working for more than 2d')
%end

% reshape into condition space (conditions on 3rd dim)
for n=1:max(cmat);
    tW(:,n,:)=W(find(cmat==n),:);  
end

%%%%%%NOW DEFUNCT: ANOVA  needs to be repeated measures
%y=reshape(tW,1,nsub*r*ncond)';
%cond=cat(1,repmat(cmat(1:nsub),ncond,1),repmat(cmat(nsub+1:nsub*2),ncond,1));       %conditions
%dim=repmat(cmat,r,1);      %dimensions  
%[out.P,out.T,out.STATS]=anovan(y,{cond dim},'full',[],['cond'; 'dim ']);

%we need to do MANOVA on difference scores (only works for ncond=2)
%y=tW(:,:,1)-tW(:,:,2);
y=reshape(tW,nsub*ncond,r);
[D,out.P,out.STATS] = manova1(y,cmat);
for n=1:size(out.P)
disp(['p-value for dim',num2str(n),'...',num2str(out.P(n))]);
end

mtW=squeeze(mean(tW));
figure;bar(mtW')
xlabel('dimension');


