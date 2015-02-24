%this script takes clusters as specified by the struct.mat file
%then uses timeseries information from the contrast images of each subject
%-->gets timeseries for specific clusters

p1 = spm_get(1,'*struct.mat','Select struct file');
[d,f,e] = fileparts(p1);
eval(['cd(''' d ''')'])
eval(['load ' p1])

P = spm_get(Inf,'*.img','Select contrast images to extract');

clusters = tor_extract_rois(P,V,V);

eval(['save ' [f e] ' V clusters'])