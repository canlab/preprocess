function F=SSE(p,func,t,tp)

%from Trisha vand Zandt's chapter

%SSE objective function for minimizing SSE
%Call as,e.g., X = fminsearch('SSE',p,[],func,t,tp)

F=norm(feval(func,t,p)-tp)^2;