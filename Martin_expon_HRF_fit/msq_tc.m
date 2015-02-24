function [m] = msq_tc(V,t,tc,delta)

dim = size(delta);
len = dim(1);
cond = dim(2);

timecourse = zeros(len,1);

for i=1:cond,
    
    s = (i-1)*8;
    HR = il_hdmf_tw(t,V((s+1):(s+8)));
    timecourseA = conv(delta(:,i), HR);
    timecourseA = timecourseA(1:len);
    
    timecourse = timecourse + timecourseA;                  % The combined time course for stimuli A and B

end;

m=sum((tc-timecourse).^2);
