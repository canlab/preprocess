%	to open physiological data file and do some plotting of
%	of the ecg and respiratory data.  
%	
%	rev 0 	11/3/99		original from motplot
%	rev 1 	12/7/99		does the resp phase from peaks
%	rev 2 	12/14/99	calculates running HR, resp rate, and
%				breath volume
%				makes an output for each acquisition
%	rev 3 	1/25/00		makes an output every desired TR
%	rev 11 	2/21/04		to read file from excite physio
%	rev 12 	2/25/04		file has trig time pts and then resp waveform
%				trig table is terminated by -9999 entry
%		2/14/05		changed the resp threshold and made the
%				output in resp rate rather then interval
%		7/26/05		fix odd factor of two in the 'ignore time'

clear; 
kern2 = [19:-1:-19];

FS = 40.;		% sampling rate of data logger, Hz
dt = 1/FS;		% in seconds

%  open the input file

fname = input('gimme file = ', 's');
fid = fopen(fname, 'r');
[dat nf] = fscanf(fid, '%g %g\n');
fclose(fid);

% parse the triggers and resp waveform

ntrig = find(dat == -9999)-1;
etrig = dat(1:ntrig)*dt;	% trigger times
respr = dat(ntrig+2:nf);	% resp waveform
n = length(respr);

tig = input('gimme number of seconds to ignore [0] = ');
if(isempty(tig))
  tig = 0; 
end
n1 = 1 + fix(tig/dt);	% start point
respr = respr(n1:n);
n = length(respr);
etrig = etrig - tig;
etrig = nonzeros(etrig.*(etrig>0));
ne = length(etrig);
aveecg = (etrig(ne) - etrig(1))/ne;
avehr = 60/aveecg;		% BPM
fprintf('average HR = %d BPM\n', fix(avehr));

ts = n*dt;		% sampled time
t = (1:n)*dt;		%

Tseq = input('output sample interval, s [1s]) = ');
if(isempty(Tseq))
  Tseq = 1;
end
nout = fix(ts/Tseq);
time = (0:nout-1)*Tseq;
fprintf('num output samples = %d\n', nout);

% find ecg intvl for each cardiac cycle and thence hrate

ecgtime = diff(etrig);
necgt = ne - 1;
tecgtime = (etrig(1:necgt) + etrig(2:ne))*.5;
tecgintvl = spline(tecgtime, ecgtime, time);
hrate = 60./tecgintvl;
subplot(2,2,1); 
plot(time, hrate);grid
ylabel('hrate, BPM');
xlabel('time, s');
title(fname);

%  now do the resp.  find the peaks

respx = max(respr);
respn = min(respr);
resp = 100*(respr - respn)/(respx - respn);
subplot(2,2,2); 
plot(t, resp);grid
ylabel('respiration'); 
xlabel('time, s');

drdt = conv(resp, kern2);
drdt = drdt(19:n+18);
d2rdt2 = conv(drdt, kern2);
d2rdt2 = d2rdt2(19:n+18);
rpeak = (d2rdt2 > 0.5e-6);  	% nice threshold

% find the resp trigs

nr = 0; 
for (j=2:n)
  if (rpeak(j)==1 & rpeak(j-1)==0)  % first only
    nr = nr + 1;
    rtrig(nr) = j*dt;
  end
end
averesp = (rtrig(nr) - rtrig(1))/nr;
averrate = 60/averesp;		% breaths/min
fprintf('average resp = %d breaths/min\n', fix(averrate));

% find resp intvl for each breath

resptime = diff(rtrig);
nrespt = nr - 1;
tresptime = (rtrig(1:nrespt) + rtrig(2:nr))*.5;
respintvl = spline(tresptime, resptime, time);
rrate = 60./respintvl;
subplot(2,2,3); 
plot(time, rrate);grid
ylabel('resp rate, BrPM');
xlabel('time, s');

% save the output

fnout = input('output file [cr for default, s for none]= ', 's');
if(isempty(fnout))
  fnout = sprintf('%s.out',fname);
  elseif(strcmp(fnout, 's'))
  return;
end
fout = fopen(fnout, 'w');
fprintf(fout, '%d %f\n', [nout Tseq]);
fprintf(fout, '%f %f\n', [hrate', rrate']');
fclose(fout); 
fprintf('wrote file  %s\n', fnout);

