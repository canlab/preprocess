%
%  Note, for this to work 'Multivolume Analyze' must be turned on
%
% $Id: WK_CovEst.m,v 1.1 2005/12/13 20:35:55 nichols Exp $

spm_defaults
global defaults
defaults.analyze.multivol = 1;

cd /group/nichols/Tom/MvarMFX
cd d:\nichols\TENlab\Tom\MvarMFX

n = 39;
p = 4;
Y = spm_get(n*p,'IMAGE',sprintf('Select %d images',n*p));

CovEst(Y,n);
Cov = spm_get(p*(p+1)/2,'Cov*IMAGE',sprintf('Select %d images',p*(p+1)/2));

save WK Y n p Cov

load WK
CovLRT(Y,n,Cov)


/group/nichols/data/WagerSwitch/020711ab/con_0008.img
/group/nichols/data/WagerSwitch/020711ab/con_0009.img
/group/nichols/data/WagerSwitch/020711ab/con_0011.img
/group/nichols/data/WagerSwitch/020711ab/con_0012.img
/group/nichols/data/WagerSwitch/020711dm/con_0008.img
/group/nichols/data/WagerSwitch/020711dm/con_0009.img
/group/nichols/data/WagerSwitch/020711dm/con_0011.img
/group/nichols/data/WagerSwitch/020711dm/con_0012.img
/group/nichols/data/WagerSwitch/020711mt/con_0008.img
/group/nichols/data/WagerSwitch/020711mt/con_0009.img
/group/nichols/data/WagerSwitch/020711mt/con_0011.img
/group/nichols/data/WagerSwitch/020711mt/con_0012.img
/group/nichols/data/WagerSwitch/020725js/con_0008.img
/group/nichols/data/WagerSwitch/020725js/con_0009.img
/group/nichols/data/WagerSwitch/020725js/con_0011.img
/group/nichols/data/WagerSwitch/020725js/con_0012.img
/group/nichols/data/WagerSwitch/020726ag/con_0008.img
/group/nichols/data/WagerSwitch/020726ag/con_0009.img
/group/nichols/data/WagerSwitch/020726ag/con_0011.img
/group/nichols/data/WagerSwitch/020726ag/con_0012.img
/group/nichols/data/WagerSwitch/020726yc/con_0008.img
/group/nichols/data/WagerSwitch/020726yc/con_0009.img
/group/nichols/data/WagerSwitch/020726yc/con_0011.img
/group/nichols/data/WagerSwitch/020726yc/con_0012.img
/group/nichols/data/WagerSwitch/020808aw/con_0008.img
/group/nichols/data/WagerSwitch/020808aw/con_0009.img
/group/nichols/data/WagerSwitch/020808aw/con_0011.img
/group/nichols/data/WagerSwitch/020808aw/con_0012.img
/group/nichols/data/WagerSwitch/020827mk/con_0008.img
/group/nichols/data/WagerSwitch/020827mk/con_0009.img
/group/nichols/data/WagerSwitch/020827mk/con_0011.img
/group/nichols/data/WagerSwitch/020827mk/con_0012.img
/group/nichols/data/WagerSwitch/020829jh/con_0008.img
/group/nichols/data/WagerSwitch/020829jh/con_0009.img
/group/nichols/data/WagerSwitch/020829jh/con_0011.img
/group/nichols/data/WagerSwitch/020829jh/con_0012.img
/group/nichols/data/WagerSwitch/020903lb/con_0008.img
/group/nichols/data/WagerSwitch/020903lb/con_0009.img
/group/nichols/data/WagerSwitch/020903lb/con_0011.img
/group/nichols/data/WagerSwitch/020903lb/con_0012.img
/group/nichols/data/WagerSwitch/020910rb/con_0008.img
/group/nichols/data/WagerSwitch/020910rb/con_0009.img
/group/nichols/data/WagerSwitch/020910rb/con_0011.img
/group/nichols/data/WagerSwitch/020910rb/con_0012.img
/group/nichols/data/WagerSwitch/020912am/con_0008.img
/group/nichols/data/WagerSwitch/020912am/con_0009.img
/group/nichols/data/WagerSwitch/020912am/con_0011.img
/group/nichols/data/WagerSwitch/020912am/con_0012.img
/group/nichols/data/WagerSwitch/020923ap/con_0008.img
/group/nichols/data/WagerSwitch/020923ap/con_0009.img
/group/nichols/data/WagerSwitch/020923ap/con_0011.img
/group/nichols/data/WagerSwitch/020923ap/con_0012.img
/group/nichols/data/WagerSwitch/020924mr/con_0008.img
/group/nichols/data/WagerSwitch/020924mr/con_0009.img
/group/nichols/data/WagerSwitch/020924mr/con_0011.img
/group/nichols/data/WagerSwitch/020924mr/con_0012.img
/group/nichols/data/WagerSwitch/021017sb/con_0008.img
/group/nichols/data/WagerSwitch/021017sb/con_0009.img
/group/nichols/data/WagerSwitch/021017sb/con_0011.img
/group/nichols/data/WagerSwitch/021017sb/con_0012.img
/group/nichols/data/WagerSwitch/021029bn/con_0008.img
/group/nichols/data/WagerSwitch/021029bn/con_0009.img
/group/nichols/data/WagerSwitch/021029bn/con_0011.img
/group/nichols/data/WagerSwitch/021029bn/con_0012.img
/group/nichols/data/WagerSwitch/021029ju/con_0008.img
/group/nichols/data/WagerSwitch/021029ju/con_0009.img
/group/nichols/data/WagerSwitch/021029ju/con_0011.img
/group/nichols/data/WagerSwitch/021029ju/con_0012.img
/group/nichols/data/WagerSwitch/021104se/con_0008.img
/group/nichols/data/WagerSwitch/021104se/con_0009.img
/group/nichols/data/WagerSwitch/021104se/con_0011.img
/group/nichols/data/WagerSwitch/021104se/con_0012.img
/group/nichols/data/WagerSwitch/021105ny/con_0008.img
/group/nichols/data/WagerSwitch/021105ny/con_0009.img
/group/nichols/data/WagerSwitch/021105ny/con_0011.img
/group/nichols/data/WagerSwitch/021105ny/con_0012.img
/group/nichols/data/WagerSwitch/021105rh/con_0008.img
/group/nichols/data/WagerSwitch/021105rh/con_0009.img
/group/nichols/data/WagerSwitch/021105rh/con_0011.img
/group/nichols/data/WagerSwitch/021105rh/con_0012.img
/group/nichols/data/WagerSwitch/021111ks/con_0008.img
/group/nichols/data/WagerSwitch/021111ks/con_0009.img
/group/nichols/data/WagerSwitch/021111ks/con_0011.img
/group/nichols/data/WagerSwitch/021111ks/con_0012.img
/group/nichols/data/WagerSwitch/021211js/con_0008.img
/group/nichols/data/WagerSwitch/021211js/con_0009.img
/group/nichols/data/WagerSwitch/021211js/con_0011.img
/group/nichols/data/WagerSwitch/021211js/con_0012.img
/group/nichols/data/WagerSwitch/021211km/con_0008.img
/group/nichols/data/WagerSwitch/021211km/con_0009.img
/group/nichols/data/WagerSwitch/021211km/con_0011.img
/group/nichols/data/WagerSwitch/021211km/con_0012.img
/group/nichols/data/WagerSwitch/021212cd/con_0008.img
/group/nichols/data/WagerSwitch/021212cd/con_0009.img
/group/nichols/data/WagerSwitch/021212cd/con_0011.img
/group/nichols/data/WagerSwitch/021212cd/con_0012.img
/group/nichols/data/WagerSwitch/021220kw/con_0008.img
/group/nichols/data/WagerSwitch/021220kw/con_0009.img
/group/nichols/data/WagerSwitch/021220kw/con_0011.img
/group/nichols/data/WagerSwitch/021220kw/con_0012.img
/group/nichols/data/WagerSwitch/030206cp/con_0008.img
/group/nichols/data/WagerSwitch/030206cp/con_0009.img
/group/nichols/data/WagerSwitch/030206cp/con_0011.img
/group/nichols/data/WagerSwitch/030206cp/con_0012.img
/group/nichols/data/WagerSwitch/030206es/con_0008.img
/group/nichols/data/WagerSwitch/030206es/con_0009.img
/group/nichols/data/WagerSwitch/030206es/con_0011.img
/group/nichols/data/WagerSwitch/030206es/con_0012.img
/group/nichols/data/WagerSwitch/030212bs/con_0008.img
/group/nichols/data/WagerSwitch/030212bs/con_0009.img
/group/nichols/data/WagerSwitch/030212bs/con_0011.img
/group/nichols/data/WagerSwitch/030212bs/con_0012.img
/group/nichols/data/WagerSwitch/030212vb/con_0008.img
/group/nichols/data/WagerSwitch/030212vb/con_0009.img
/group/nichols/data/WagerSwitch/030212vb/con_0011.img
/group/nichols/data/WagerSwitch/030212vb/con_0012.img
/group/nichols/data/WagerSwitch/030313kh/con_0008.img
/group/nichols/data/WagerSwitch/030313kh/con_0009.img
/group/nichols/data/WagerSwitch/030313kh/con_0011.img
/group/nichols/data/WagerSwitch/030313kh/con_0012.img
/group/nichols/data/WagerSwitch/030320cg/con_0008.img
/group/nichols/data/WagerSwitch/030320cg/con_0009.img
/group/nichols/data/WagerSwitch/030320cg/con_0011.img
/group/nichols/data/WagerSwitch/030320cg/con_0012.img
/group/nichols/data/WagerSwitch/030325ts/con_0008.img
/group/nichols/data/WagerSwitch/030325ts/con_0009.img
/group/nichols/data/WagerSwitch/030325ts/con_0011.img
/group/nichols/data/WagerSwitch/030325ts/con_0012.img
/group/nichols/data/WagerSwitch/030326mp/con_0008.img
/group/nichols/data/WagerSwitch/030326mp/con_0009.img
/group/nichols/data/WagerSwitch/030326mp/con_0011.img
/group/nichols/data/WagerSwitch/030326mp/con_0012.img
/group/nichols/data/WagerSwitch/030326my/con_0008.img
/group/nichols/data/WagerSwitch/030326my/con_0009.img
/group/nichols/data/WagerSwitch/030326my/con_0011.img
/group/nichols/data/WagerSwitch/030326my/con_0012.img
/group/nichols/data/WagerSwitch/030327jk/con_0008.img
/group/nichols/data/WagerSwitch/030327jk/con_0009.img
/group/nichols/data/WagerSwitch/030327jk/con_0011.img
/group/nichols/data/WagerSwitch/030327jk/con_0012.img
/group/nichols/data/WagerSwitch/030402af/con_0008.img
/group/nichols/data/WagerSwitch/030402af/con_0009.img
/group/nichols/data/WagerSwitch/030402af/con_0011.img
/group/nichols/data/WagerSwitch/030402af/con_0012.img
/group/nichols/data/WagerSwitch/030402st/con_0008.img
/group/nichols/data/WagerSwitch/030402st/con_0009.img
/group/nichols/data/WagerSwitch/030402st/con_0011.img
/group/nichols/data/WagerSwitch/030402st/con_0012.img
/group/nichols/data/WagerSwitch/030408ma/con_0008.img
/group/nichols/data/WagerSwitch/030408ma/con_0009.img
/group/nichols/data/WagerSwitch/030408ma/con_0011.img
/group/nichols/data/WagerSwitch/030408ma/con_0012.img
/group/nichols/data/WagerSwitch/030410ca/con_0008.img
/group/nichols/data/WagerSwitch/030410ca/con_0009.img
/group/nichols/data/WagerSwitch/030410ca/con_0011.img
/group/nichols/data/WagerSwitch/030410ca/con_0012.img



d:\nichols\data2\tmp\Wager\Inhib\020711ab\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020711ab\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020711ab\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020711ab\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020711dm\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020711dm\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020711dm\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020711dm\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020711mt\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020711mt\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020711mt\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020711mt\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020725js\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020725js\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020725js\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020725js\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020726ag\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020726ag\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020726ag\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020726ag\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020726yc\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020726yc\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020726yc\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020726yc\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020808aw\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020808aw\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020808aw\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020808aw\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020827mk\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020827mk\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020827mk\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020827mk\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020829jh\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020829jh\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020829jh\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020829jh\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020903lb\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020903lb\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020903lb\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020903lb\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020910rb\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020910rb\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020910rb\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020910rb\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020912am\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020912am\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020912am\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020912am\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020923ap\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020923ap\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020923ap\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020923ap\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\020924mr\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\020924mr\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\020924mr\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\020924mr\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021017sb\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021017sb\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021017sb\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021017sb\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021029bn\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021029bn\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021029bn\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021029bn\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021029ju\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021029ju\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021029ju\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021029ju\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021104se\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021104se\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021104se\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021104se\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021105ny\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021105ny\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021105ny\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021105ny\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021105rh\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021105rh\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021105rh\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021105rh\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021111ks\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021111ks\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021111ks\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021111ks\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021211js\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021211js\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021211js\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021211js\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021211km\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021211km\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021211km\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021211km\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021212cd\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021212cd\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021212cd\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021212cd\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\021220kw\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\021220kw\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\021220kw\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\021220kw\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030206cp\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030206cp\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030206cp\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030206cp\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030206es\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030206es\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030206es\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030206es\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030212bs\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030212bs\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030212bs\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030212bs\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030212vb\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030212vb\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030212vb\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030212vb\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030313kh\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030313kh\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030313kh\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030313kh\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030320cg\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030320cg\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030320cg\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030320cg\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030325ts\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030325ts\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030325ts\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030325ts\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030326mp\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030326mp\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030326mp\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030326mp\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030326my\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030326my\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030326my\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030326my\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030327jk\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030327jk\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030327jk\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030327jk\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030402af\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030402af\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030402af\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030402af\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030402st\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030402st\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030402st\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030402st\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030408ma\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030408ma\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030408ma\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030408ma\con_0012.img
d:\nichols\data2\tmp\Wager\Inhib\030410ca\con_0008.img
d:\nichols\data2\tmp\Wager\Inhib\030410ca\con_0009.img
d:\nichols\data2\tmp\Wager\Inhib\030410ca\con_0011.img
d:\nichols\data2\tmp\Wager\Inhib\030410ca\con_0012.img


Find Bcov

cd d:\nichols\data2\tmp\Wager\Inhib\020711ab

clear
load SPM xX
load xCon
C = [xCon([8 9 11 12]).c];
Bcov = C'*xX.Bcov*C;

Bcor = diag(sqrt(1./diag(Bcov)))*Bcov*diag(sqrt(1./diag(Bcov)));

% Bcor =
% 
%     1.0000   -0.1086   -0.0000    0.0000
%    -0.1086    1.0000    0.0000   -0.0000
%    -0.0000    0.0000    1.0000   -0.1203
%     0.0000   -0.0000   -0.1203    1.0000


Good!  Nearly identity
