function z6=getindscal;

close all
fid=fopen('c:\matlab6p5\work\indscal\indscalx\justdata.txt','r');
z3=zeros(224,16);
counter=0;
for n=1:210;
clear z1;clear z2;
z1=fgetl(fid);
z2=z1(18:2:length(z1));
if (n-1)/15==fix((n-1)/15);
    counter=counter+1;
end
for z=1:length(z2);    
    z3(n+counter,z)=str2num(z2(z));
end

end

fclose(fid);

z4=shiftdim(reshape(z3,16,14,16),1);
for n=1:16;
    z5(:,n,:)=squeeze(z4(:,:,n));
end
z5=shiftdim(z5,1);

%make symmetric
for n=1:14;
        z5(:,:,n)=z5(:,:,n)+rot90(fliplr(tril(squeeze(z5(:,:,n)))),1);
end

%similarity to dissimilarity
z5=10-z5;

for n=1:14;    
    a=z5(:,:,n);
    ssa=sum(a(:).^2);
    nssa=sqrt((a.^2)./ssa);
    z6(:,:,n)=nssa;
end





