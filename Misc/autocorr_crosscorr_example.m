for i= 1:1000,
n = noisevector(100,[1 .5 .3 .2 .1 0],1);
n2 = noisevector(100,[1 0],1);
c(i,:) = xcorr(n,n2);
end
figure;imagesc(c);
colorbar
for i= 1:1000,
n = noisevector(100,[1 0],1);
n2 = noisevector(100,[1 0],1);
c2(i,:) = xcorr(n,n2);
end
figure;imagesc(c);
figure;imagesc(c2);colorbar