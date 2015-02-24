
load DATA
DATA.SPEC.ndims = 'choose'; DATA.SPEC.clustsolutions = 'choose'; 

%DATA.SPEC.robustflag = 1; DATA.SPEC.betaflag = 0;


try,DATA.DATA = rmfield(DATA.DATA,'b');,catch,end

try,DATA.DATA = rmfield(DATA.DATA,'filtered_dat');,catch,end
try,DATA.DATA = rmfield(DATA.DATA,'resids');,catch,end
try,DATA.DATA = rmfield(DATA.DATA,'fits');,catch,end
try,DATA.DATA = rmfield(DATA.DATA,'fir');,catch,end
try,DATA.DATA = rmfield(DATA.DATA,'corr');,catch,end
try,DATA.DATA = rmfield(DATA.DATA,'xc');,catch,end

try,DATA = rmfield(DATA,'NDIMS');,catch,end
try,DATA = rmfield(DATA,'INDSCAL');,catch,end
try,DATA = rmfield(DATA,'CLUSTER');,catch,end
try,DATA = rmfield(DATA,'CORRELS');,catch,end
try,DATA = rmfield(DATA,'APPLY_CLUSTER');,catch,end

DATA = mvroi(DATA,DATA.SPEC);