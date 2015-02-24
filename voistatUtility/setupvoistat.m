%eval(['addpath /data1/intext/sub' num2str(snum) '/task'])
eval(['cd /data1/intext/sub' num2str(snum) '/task'])

if snum == 0
addpath /data1/intext/sub1/task
addpath /data1/intext/sub2/task
addpath /data1/intext/sub4/task
addpath /data1/intext/sub5/task
addpath /data1/intext/sub6/task
addpath /data2/intext/sub8/task
addpath /data2/intext/sub9/task
addpath /data4/intext/sub10/task
addpath /data4/intext/sub11/task
addpath /data4/intext/sub12/task
end

disp('Your file is...');drawnow
!ls sravol*.img | grep 1000
