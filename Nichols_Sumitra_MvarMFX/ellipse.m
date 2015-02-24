function hh=ellipse(varargin)
%ELLIPSE   Draws an ellipse.
%   ELLIPSE(F1,F2,R,[,LINESPEC])
%   ELLIPSE(O,A,B[,LINESPEC])
%
%   H = ELLIPSE(...)
%
%   See also CIRLCE, FCIRCLE.

% Copyright (c) 2001-09-15, B. Rasmus Anthin.

ni=nargin;
opt='';
if nargin & ischar(varargin{end})
   opt=varargin{end};
   ni=ni-1;
end
error(nargchk(3,3,ni))
if length(varargin{2})>1
   f1=varargin{1};
   f2=varargin{2};
   r=varargin{3}*2;
   if abs(f1(1)-f2(1))<abs(f1(2)-f2(2))
      m=mean([f1(1) f2(1)]);
      f1(1)=m;f2(1)=m;
   else
      m=mean([f1(2) f2(2)]);
      f1(2)=m;f2(2)=m;
   end
   o=[mean([f1(1) f2(1)]) mean([f1(2) f2(2)])];
   c=abs(f2-o);
   if ~c(1)
      b=max([r c(2)]);
      a=sqrt(b^2-c(2)^2);
   elseif ~c(2)
      a=max([r c(1)]);
      b=sqrt(a^2-c(1)^2);
   end
   %hold on,plot(f1(1),f1(2),'x'),plot(f2(1),f2(2),'x')
else
   o=varargin{1};
   a=varargin{2};
   b=varargin{3};
end
phi=linspace(0,2*pi);
h=plot(o(1)+a*cos(phi),o(2)+b*sin(phi),opt);
set(h,'user',{'ellipse',abs(a),abs(b),o});
if nargout,hh=h;end
