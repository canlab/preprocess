function dat1=loadjakobo;

load('jakobo');
dat=reshape(dat,15,16,15);
dat=shiftdim(dat,2);

for n=1:16;
    d=dat(:,:,n);
    d1=d*0;
    d1=d1+tril(d);
    d1=d1+rot90(fliplr(tril(d)),1);
    dat1(:,:,n)=d1;
end


%this line writes the data in case you need to
%fid=fopen('jakoff.txt','w');for s=1:16;for n=1:15;fprintf(fid,' %3.0f
%%3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f\n',squeeze(dat(n,:,s)));end;end;fclose(fid);
%fid=fopen('jakoff1.txt','w');for s=1:16;for n=1:15;for n1=1:n-1;outy=squeeze(dat(n,:,s));fprintf(fid,' %3.0f',outy(n1));endfprintf(fid,'\n');end;end;fclose(fid);fclose(fid);