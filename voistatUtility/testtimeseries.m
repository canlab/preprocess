a = zeros(8,6,4);
for i = 1:8
for j = 1:6
for k = 1:4
a(i,j,k) = m; m = m+1;
end
end
end
hdr.xdim = 8;hdr.ydim=6;hdr.zdim=4;
hdr.datatype = 4;
make_img('test',a,hdr)
[b,hdr] = readim2('test');

figure;imagesc(a(:,:,2))

voxel = [6 7 4]
ts = timeseries('voxel','test',1,[voxel(1) voxel(2) voxel(3)]),a(voxel(2),voxel(1),voxel(3)),b(voxel(2),voxel(1),voxel(3))
voxel = [6 2 4]
ts = timeseries('voxel','test',1,[voxel(1) voxel(2) voxel(3)]),a(voxel(2),voxel(1),voxel(3)),b(voxel(2),voxel(1),voxel(3))


