
n   = 39;
Lab = {'Ext:Obj','Ext:Att','Int:Obj','Int:Att'};
save WK Lab -append

ViewCov('init','Covar.img',Lab,n)

load WK
ViewCovY('init',Y,'Covar.img',Lab)
