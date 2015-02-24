function ccOUT=testsigL(ldat,oXc,cmat,X,names,compmat,numperm);

% test significance BETWEEN clusters
counter=0;
for c=1:max(oXc);
    for c1=c+1:max(oXc);
    which=find(c==oXc);
    which1=find(c1==oXc);    
    diff=squeeze(sum(sum(ldat(which,which1,:))));
    [h t]=ttest(diff);
    p=1-tcdf(t,length(cmat)-1);    
    disp(['between ',num2str(c),' and ',num2str(c1),' t=',num2str(t),' p<',num2str(p)]);
    end
end
    

return













ncond=max(cmat);
nsub=size(ldat,3)/ncond;
nchan=size(ldat,1);

ldat=reshape(ldat,nchan,nchan,nsub,ncond);
cd=shiftdim(ldat,2);

%reshape cd to put it in same order as oXc
totwhich=0;
totwhich1=1;
for c=1:max(oXc)
    which=find(c==oXc);
    totwhich1=0;
    for c1=1:max(oXc)
        which1=find(c1==oXc);
        cd1(:,:,totwhich+1:totwhich+length(which),totwhich1+1:totwhich1+length(which1))=cd(:,:,which,which1);
        totwhich1=totwhich1+length(which1);
    end
    totwhich=totwhich+length(which);
end

% do random permutation testing on all channels
[m p]=randpermeeg_nocor(cd1,compmat,numperm,0.05);
pm=((1000*0.05)-p)/1000;
viz(m);colormap jet;
viz(pm);colormap(flipud(hot));
sc=numperm/1000;
figure;
for n=1:size(compmat,1);
condnames{n}=num2str(compmat(n,:));
end
for n=1:size(p,1);
    sigmat=squeeze(pm(n,:,:)<0.05);
    subplot(1,size(p,1),n);
    nmdsfig1(X,oXc,names,sigmat,(squeeze(p(n,:,:)))/sc);
    title([char(condnames(n))]);
end



% test significance within each cluster
for c=1:max(oXc);
    which=find(c==oXc);
    if length(which)>1;
    clust_data=squeeze(cdat(which,which,:,:));  
    for s=1:size(clust_data,3);
        for cond=1:size(clust_data,4);
            cd=squeeze(clust_data(:,:,s,cond));
            acd(s,cond)=squeeze(mean(cd(cd~=1)));
                   
        end
    end
        Fw=RManova(acd,compmat);
        pFw=1-fcdf(Fw,size(acd,2)-1,size(acd,1));
        for n=1:size(Fw);
        disp(['comp ',num2str(n),' within ',num2str(c),' F ',num2str(squeeze(Fw(n))),' p<',num2str(squeeze(pFw(n)))]);
    end
    end
 end
 
 
% test significance BETWEEN clusters
counter=0;
for c=1:max(oXc);
    for c1=c+1:max(oXc);
    which=find(c==oXc);
    which1=find(c1==oXc);    
    clust_data=mean(mean(cdat(which,which1,:,:),1),2);  
    bcd=squeeze(clust_data);
   
    Fb=RManova(bcd,compmat);
    pFb=1-fcdf(Fb,size(acd,2)-1,size(acd,1)); 

        for n=1:size(Fb);
        disp(['comp ',num2str(n),' between ',num2str(c),' and ',num2str(c1),' F ',num2str(Fb(n)),' p<',num2str(pFb(n))]);
        end
    end
end
    

ccOUT.Fb=Fb;
ccOUT.Fw=Fw;