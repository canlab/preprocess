function m=fitfun_logit(V,delta, t, tc)


dim = size(delta);
len = dim(1);
cond = dim(2);

timecourse = zeros(len,1);

for i=1:cond,
    
    s = (i-1)*7;
    HR = il_hdmf_tw2(t,V((s+1):(s+7)));
    timecourseA = conv(delta(:,i), HR);
    timecourseA = timecourseA(1:len);
    
    timecourse = timecourse + timecourseA;                  % The combined time course for stimuli A and B

end;

m=(tc-timecourse);