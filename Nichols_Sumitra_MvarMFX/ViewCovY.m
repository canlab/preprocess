function varargout = ViewCovY(action,varargin)
%
% ViewCovY('init','Dat.img','Corr.img',{labels})
%    Initialize
%
% ViewCovY('setloc',mm[,fig])   ViewCovY('setmm',mm[,fig])
%    Set location in mm
%
% ViewCovY('update'[,f])
%    Update; only useful if VC structure has been modified manually.
%
% ViewCovY('SetCoords',mm,h,hReg)
%    See spm_XYZreg; used for interacting with spm_orthviews
%
%
% $Id: ViewCovY.m,v 1.8 2005/12/15 18:36:13 nichols Exp $

if (nargin<1), action = 'init'; end

switch lower(action)

 case 'init'
 %------------------------------------------------------------------------
  if (nargin >= 2)
    Dat = varargin{1};
  else
    Dat = spm_get(Inf,'IMAGE','Select all n*p images');
  end
  if (nargin >= 3), 
    Cor = varargin{2};
  else
    Cor = spm_get(1,'IMAGE','Select Correlation Image'); 
  end
  if (nargin >= 4)
    Lab = varargin{3};
  else
    Lab = GetLabels(p_from_k(Nvol(Cor)));
  end

  f = CreateVC(Dat,Cor,Lab);
  Update(f);

 case 'setcoords'
 %------------------------------------------------------------------------
  mm = varargin{1};
  f = varargin{2};
  hReg = varargin{3};
  % f = GoodFig([],0);

  SetLoc(f,'mm',mm);
  Update(f);

 case 'setmm', 'setloc'
 %------------------------------------------------------------------------
  mm = varargin{1};
  if nargin>=3
    f = varargin{2};
  else
    f = GoodFig([],1); 
  end

  SetLoc(f,'mm',mm);
  Update(f)
  
 case 'setvox'
 %------------------------------------------------------------------------
  vox = varargin{1};
  if nargin>=3
    f = varargin{2};
  else
    f = GoodFig([],1); 
  end

  SetLoc(f,'vox',vox);
  Update(f)

 case 'orthview_update'
 %------------------------------------------------------------------------
  vox = varargin{1};
  if nargin>=3
    f = varargin{2};
  else
    f = GoodFig([],1); 
  end

  mm = spm_orthviews('pos')';
  SetLoc(f,'mm',mm);
  Update(f)

 case 'update'
 %------------------------------------------------------------------------
  if nargin>=3
    f = varargin{2};
  else
    f = GoodFig([],1); 
  end

  Update(f);

 case 'togglereg'
 %------------------------------------------------------------------------
  f = varargin{1};

  ToggleReg(f);

 case 'register'
 %------------------------------------------------------------------------
  Register(varargin{1});

 otherwise
 %------------------------------------------------------------------------
  error(['Unknown action: ''' action ''''])
end

return



function f = CreateVC(Dat,Cor,Lab)
% Cor - Correlation image filename
% Lab - Length p cell array of condition labels
% n   - Number of subjects

% Preliminaries
%   Find p, map correlation elements
%-------------------------------------------------------------------------

np      = size(Dat,1);
k       = Nvol(Cor);
[n,p]   = np_from_jk(np,k);
Corr    = eye(p);
CorrP   = repmat(NaN,p,p);
Data    = repmat(NaN,p,n);
Mean    = repmat(NaN,p,1);
Stdev   = ones(1,p);
Q       = BuildQ(p);

if ~isempty(findstr(Cor,','))
  Cor(min(find(Cor==',')):end) = [];
end

for i=1:np
  Vdat(i) = spm_vol(Dat(i,:));
end

for i=1:k
  Vcor(i) = spm_vol(sprintf('%s, %d',Cor,i));
end


% Create default/empty Figure
%-------------------------------------------------------------------------

f = figure('position',[100 100 700 600],...
	   'MenuBar','none',...
	   'Tag','ViewCovY');
colormap(winter);

% Set button
hui = uicontrol(f,'String','',...
		'Style','PushButton',...
		'Position',[10/2 7/2 100 20 ],...
		'CallBack','ViewCovY(''ToggleReg'',gcbf)',...
		'Visible','on');
	  

% Image
hi = imagesc(zeros(p,p),[-1 1]);
ha = gca;
axis image
set(hi,'ButtonDownFcn',['VC=get(gcbf,''UserData'');'...
	    'Stdev=VC.Stdev,Cor=VC.Corr,CorP=VC.CorrP']);

set(ha,'Position',[.15 .21 0.65 0.65]/2,...
       'Xtick',1:p,'Ytick',1:p,...
       'XtickLabel',Lab,'YtickLabel',Lab);
hc = colorbar;
set(hc,'Position',[.85 .21 0.06 0.65]/2);

% Value labels
hvt = zeros(p,p);
for x=1:p
  for y=1:p
    hvt(x,y) = text(x,y,'',...
		    'HorizontalAlign','Center',...
		    'VerticalAlign','middle',...
		    'FontWeight','bold');
    if (x>y)
      set(hvt(x,y),'Color','white');
    elseif (x==y)
      set(hvt(x,y),'Color','Black');
    else
      set(hvt(x,y),'Color','yellow');
    end
  end
end

% Status Text: Location
ha_text = axes('Position',[0 0 1 1],'Visible','off');
ht_loclab = text(.4/2,.01/2,...
	      {'mm:','vox:'},...
	      'VerticalAlign','Bottom');
ht_loc = text(.5/2,.01/2,...
	      {sprintf('%g %g %g',[NaN NaN NaN]),...
	       sprintf('%d %d %d',[NaN NaN NaN])},...
	      'VerticalAlign','Bottom');

% Status Text: LR value
X2 = NaN;
X2p = NaN;
ht_x2lab  = text(.03/2,.96/2,'LRT \chi^2 =');
ht_x2     = text(.15/2,.96/2,sprintf('%g',X2));
ht_x2plab  = text(.03/2,.90/2,'LRT P  =');
ht_x2p     = text(.15/2,.90/2,sprintf('%g',X2p));

% Legend 
ThBonf = 0.05/(p*(p-1)/2);
LegAtt = {'fontweight','bold','Color'};
hLeg(1) = text(.30/2,.97/2,'Correlation',LegAtt{:},'White');
hLeg(2) = text(.30/2,.93/2,'Standard Deviation',LegAtt{:},'Black');
hLeg(3) = text(.30/2,.89/2,'P-value for Corr',LegAtt{:},'Yellow');
hLeg(4) = text(.6/2,.89/2,sprintf('Bonf Thresh = %f',ThBonf));



% Make scatterplot matrix
%-------------------------------------------------------------------------
hSM = MyScatPltMtx('init',Data);

h  = struct('i',hi,...            % Image
	    'a',ha,...            % Image axis
	    'c',hc,...            % Colorbar
	    'ui',hui,...          % Pushbutton
	    'tv',hvt,...          % Image value text
	    'ta',ha_text,...      % Whole-figure axis, for text
	    'tll',ht_loclab,...   % mm/vox label text
	    'tl',ht_loc,...       % mm/vox value text
	    'tx2l',ht_x2lab,...   % LR label text
	    'tx2',ht_x2,...       % LR stat text
	    'tx2pl',ht_x2plab,... % LR P label text
	    'tx2p',ht_x2p,...     % LR P stat text
	    'leg',hLeg,...        % Misc legend text
	    'SM',hSM...           % Scatterplot handles
	    );

VC = struct('p',p,...
	    'n',n,...
	    'mat',Vcor(1).mat,...
	    'Vdat',Vdat,...
	    'Vcor',Vcor,...
	    'Data',Data,...
	    'Mean',Mean,...
	    'Corr',Corr,...
	    'CorrSig',CorrP,...
	    'CorrSigTh',ThBonf,...
	    'Stdev',Stdev,...
	    'X2',X2,...
	    'X2p',X2p,...
	    'mm',[NaN NaN NaN],...
	    'vox',[NaN NaN NaN],...
	    'h',h,...
	    'Q',{Q},...
	    'Reg',[]);

SetVC(f,VC);

SetLoc(f,'vox',round(Vcor(1).dim(1:3)/2));

VC = GetVC(f);

Register(f);

return

function hReg = Register(f)
%
%  Register to link-up with spm_orthviews (possibly creating a registry
%  in the process).
%
VC = GetVC(f);

global st

mm = [];
if isempty(st)
  % No orthviews running
  VC.Reg = [];
else
  mm = spm_orthviews('pos')';
  % Create a registry if needed
  if isempty(VC.Reg)
    hReg = [];
    if isfield(st,'registry') 
      % Orthviews already knows of a Registry, use it
      hReg = st.registry.hReg; 
      try 
	tmp=get(hReg);
      catch 
	% Orthviews had a dud-registry... kill it
	st = rmfield(st,'registry');
	hReg = [];
      end
    end
    if isempty(hReg)
      % Orthviews not registered (or had dud registry); create one and
      % let orthviews know 
      hReg = VC.h.ta; % Use our text axis as registry
      [hReg,xyz] = spm_XYZreg('InitReg',hReg,...
			      VC.Vcor(1).mat,VC.Vcor(1).dim(1:3)',mm');
      spm_orthviews('register',hReg);
    end
  else
    hReg = VC.Reg.hReg;
  end
  % Tell the registry we're here
  hMe = f;
  spm_XYZreg('Add2Reg',hReg,hMe,'ViewCovY');
  VC.Reg = struct('hReg',hReg,'hMe',hMe);
end
SetVC(f,VC);
if ~isempty(mm)
  SetLoc(f,'mm',mm);  % Done this way to make sure vox is set correctly
end
return


function UnRegister(f)
%
%  Un-register ourselves.  Doesn't kill registry
%
%
VC = GetVC(f);

if ~isempty(VC.Reg)
  spm_XYZreg('Del2Reg',VC.Reg.hReg,VC.Reg.hMe);
  VC.Reg.hMe = [];
end

SetVC(f,VC);

return


function SetLoc(f,type,loc)
VC = GetVC(f);

M = VC.mat;

switch(type)
 case 'mm'
  mm = loc(:)';
  vox = M\[mm';1];
  vox = vox(1:3)';
 case 'vox'
  vox = loc(:)';
  mm = M*[vox';1];
  mm = mm(1:3)';
end

VC.mm = mm;
VC.vox = vox;

SetVC(f,VC);

return


function Update(f)

VC    = GetVC(f);
Q     = VC.Q;
p     = VC.p;
n     = VC.n;
vox   = VC.vox;
Vcor  = VC.Vcor;
Vdat  = VC.Vdat;
hvt   = VC.h.tv;
PTh   = VC.CorrSigTh;

%
% Read stdev, corr
%
if any(isnan(vox))
  Corr = eye(p);
  Stdev = ones(1,p);
  Data = zeros(p,n);
else
  Stdev   = spm_get_data(Vcor(1:p),vox')';
  CorDat  = spm_get_data(Vcor(p+1:end),vox');
  Data    = reshape(spm_get_data(Vdat,vox'),p,n);
  Mean    = mean(Data,2);

  Corr = eye(p);
  CorrP = eye(p); 
  for i=p+1:length(Q)
    Corr = Corr + Q{i}*CorDat(i-p);
    % corr->z Ref: Morrison, 1990, pg 119
    z = sqrt(n-3)*atanh(CorDat(i-p));
    CorrP = CorrP + Q{i}*(1-spm_Ncdf(abs(z)))*2;
  end
  for i=1:p
    CorrP(i,i) = NaN;
  end
end

SetRegBut(f);

%
% Update image
%
set(VC.h.i,'CData',CombCorStd(Corr,Stdev));

% Display image values
tmp = tril(Corr)+diag(Stdev)+triu(CorrP,+1);
for x=1:p
  for y=1:p
    set(hvt(x,y),'String',sprintf('%.3f',tmp(x,y)));
    if (x~=y)
      if (CorrP(x,y)<=PTh)
	set([hvt(x,y),hvt(y,x)],'FontWeight','Bold')
      else
	set([hvt(x,y),hvt(y,x)],'FontWeight','Normal')
      end
    end
  end
end

% Update scatter plot matrix
MyScatPltMtx('update',VC.h.SM,Data,Mean,Stdev,Corr);

%
% Compute LR
%

VC.X2 = -(n-1-(2*p+5)/6)*log(det(Corr)); % Bartlett's correction
VC.X2p = 1-spm_Xcdf(VC.X2,p*(p-1)/2);


% Display LR
set(VC.h.tx2,'String',sprintf('%g',VC.X2));
set(VC.h.tx2p,'String',sprintf('%.5f',VC.X2p));

% Display location
set(VC.h.tl,'String',...
	    {sprintf('%.2f %.2f %.2f',VC.mm),...
	     sprintf('%.2f %.2f %.2f',round(VC.vox))});


VC.Corr = Corr;
VC.CorrP = CorrP;
VC.Stdev = Stdev;
VC.Data = Data;
VC.Mean = Mean;

SetVC(f,VC)

return



function VC = GetVC(f)
% Basic info stored in figure 
VC = get(f,'UserData');
return

function SetVC(f,VC)
% Basic info stored in figure 
set(f,'UserData',VC);
return


function n = Nvol(V)

if ~isstruct(V), 
  V = spm_vol(V); 
  spm_close_vol(V);
end

if length(V)>1
  n = length(V);
else
  fp   = fopen(V.fname);
  fseek(fp,0,'eof');
  Len  = ftell(fp);
  fclose(fp);
  n    = Len/(prod(V.dim(1:3))*spm_type(V.dim(4),'bits')/8);
end

return



function p = p_from_k(k)
% Computes p, dimension of covariance matrix, from k, number of unique
% elements. 

p = (-1+sqrt(1+8*k))/2;

return

function [n,p] = np_from_jk(j,k)
% Computes n, number of subjects, and p, dimension of covariance matrix,
% from j, total number of data points (n*p) and k, number of unique 
% covariance elements. 

p = (-1+sqrt(1+8*k))/2;
n = j/p;

return



function Lab = GetLabels(p)

Lab = cell(1,p);
for i = 1:p
  str = sprintf('Enter label %d (of %d)',i,p); 
  def = sprintf('Cond %d',i); 
  Lab{i} = spm_input(str,'+1','s',def);
end

return



function g = GoodFig(f,crash)
%
% Check that f is a ViewCovY figure, and optionally crash out
%
% If f is empty, try to return a valid figure handle
%
if nargin<2, crash = 0; end

fs = findobj('Tag','ViewCovY');


if isempty(fs)
  g = 0;
else
  if isempty(f)
    g = fs(1);
  else
    g = any(f==fs);
  end
end

if crash & g==0 
  error('Selected figure is not a ViewCovY figure!')
end

return


function Q = BuildQ(p)
%
%  Need "n=2" trick spm_non_sphericity into creating Vi.
%
%  Maybe should read this in from a Matfile
%
n = 2;

I = repmat(1:n,p,1);
I = [I(:) repmat(1:p,1,n)' zeros(n*p,2)];
var = [0 1 0 0];
dep = [0 1 0 0];
xVi = struct('I',I,'var',var,'dep',dep);
xVi = spm_non_sphericity(xVi);

nQ = length(xVi.Vi);
Q = cell(1,nQ);
for i = 1:nQ
  tmp = xVi.Vi{i};
  Q{i} = tmp(1:p,1:p);
end

return


function CS = CombCorStd(Corr,Stdev,Crng)
p = size(Corr,1);
if nargin<3
  Crng = [min(Corr-eye(p)) max(Corr-eye(p))];
  if diff(Crng)==0, 
    Crng = [-1 1];
  end
end
Srng = [min(Stdev) max(Stdev)];
if diff(Srng)==0,
  if Srng(1)==0
    Srng = [0 1];
  else
    Srng = max(Srng+[-0.5 0.5],[0 0]);
  end
end

CS = Corr - eye(p); % + (diag(Stdev)-Srng(1))/(Srng(2)-Srng(1))*(Crng(2)-Crng(1));

return




function SetRegBut(f)
VC = GetVC(f);

if ~isempty(VC.Reg) & ~isempty(VC.Reg.hMe)
  set(VC.h.ui,'String','Un-link');
else
  set(VC.h.ui,'String','Link');
end

return

function ToggleReg(f)
VC = GetVC(f);

if ~isempty(VC.Reg) & ~isempty(VC.Reg.hMe)
  UnRegister(f)
  SetRegBut(f);
else
  Register(f);
  Update(f);  % Need to both change button *and* update data
end


return


function [hSM] = MyScatPltMtx(action,varargin)
% hSM - Struct of handles: h,ha,hsb,hz
%  .h matrix of plotting point handles
%  .hs matrix of axes handles
%  .hz matrix of the zero lines
%  .hsb the 'big' axis handle
%
% MyScatPltMtx('init',Data)
% MyScatPltMtx('update',hSM,Data,Mean,Std,Corr)  Struct of handles

switch lower(action)

 case 'init'

  Dat = varargin{1};
  p = size(Dat,1);

  hold on
  [h,ha,hsb] = gplotmatrix(Dat');
  f = get(hsb,'parent');
  axis(ha(:),'auto');
  hold off
  for j = 1:p
    for i = (j):p
      delete([h(i,j) ha(i,j)]);
      h(i,j) = 0;
      ha(i,j) = 0;
    end
  end

  set(ha(1,2),'YtickLabelMode','auto')
  set(ha(p-1,p),'XtickLabelMode','auto')
  hz = repmat(NaN,[p,p,2]);  % zero lines
  he = repmat(NaN,[p,p]);  % ellipses
  for i = 1:p
    for j = (i+1):p
      PlotHere(ha(i,j));
      he(i,j) = myellipse('init',ha(i,j));
      ZerAtt = {'LineStyle',':','Color',[.6 .6 .6]};
      hz(i,j,:) = [ line(NaN*[1 1],[0 0],ZerAtt{:})
		    line([0 0],NaN*[1 1],ZerAtt{:}) ];
      hh = get(ha(i,j),'children');
      set(ha(i,j),'Children',hh([4 1 2 3]));
    end
  end

  hSM = struct('h',h,'ha',ha,'hz',hz,'he',he,'hsb',hsb);

 case 'update'
  hSM = varargin{1};
  h = hSM.h;
  ha = hSM.ha;
  hz = hSM.hz;
  he = hSM.he;
  Dat  = varargin{2};
  Mean = varargin{3};
  Std  = varargin{4};
  Corr = varargin{5};
  p = size(Dat,1);
  f = get(hSM.hsb,'parent');

  Rng = [min(Dat(:)) max(Dat(:))];
  Rng = Rng+0.05*diff(Rng)*[-1 1];
  if any(isnan(Rng)) | all(Rng==0)
    Rng = [-.5 .5]; 
  end

  for i = 1:p
    for j = (i+1):p
      set(hz(i,j,1),'Xdata',Rng);
      set(hz(i,j,2),'Ydata',Rng);
      myellipse('update',he(i,j),Mean([j i]),Std([j i]),Corr(j,i));
      set(ha(i,j),'Xlim',Rng,'Ylim',Rng);
      set(h(i,j),'Xdata',Dat(j,:),'Ydata',Dat(i,:));
    end
  end
  if any(isnan([Rng])), Rng = [0 1]; end
%%% For histograms
%  for i=1:p
%    set(ha(i,i),'Xlim',Rng,'Ylim',Rng);
%  end

 otherwise 
  error(['Unknown MyScatPltMtx action: ''' action ''''])
end

return




function varargout = myellipse(action,varargin)
%
%  Have X ~ N(mu, Sigma)
% 
%  Want A s.t. P(X \in A) = 1-alpha
%
%  Note that (X-mu)'*inv(Sigma)*(X-mu) ~ \chi^2_p
%  So A = {x : (x-mu)'*inv(Sigma)*(x-mu) < F^{-1}_{\chi^2_p}(1-alpha)}
%  is what we want.  Let r2 = F^{-1}_{\chi^2_p}(1-alpha)
%
%  Also, make life easier by standardizing.  Note that
%
%       (X-mu)'*inv(Sigma)*(X-mu) =  (X-mu)'*inv(S)*inv(C)*inv(S)*(X-mu) 
%
%  where S = sqrt(diag(Sigma)) and C = correlation matrix (i.e. S*C*S =
%  Sigma).  So let Z = inv(S)*(X-mu), and then the definition of the
%  elipse that defines the (1-alpha) confidence region is
%
%       Z'*inv(C)*Z = r2
%
%  For p = 2, C=[1 rho;rho 1] and it's inverse is simple.  So with the
%  quadratic formula, if you fix z2, you can solve for z1; I get
%
%            z1 = \rho z2 +/- sqrt(1-rho^2)*sqrt(r2-z2^2)
%
%  Notably, we only will have real solutions when |z|<=sqrt(r2)
%

hnPoint = 50;  % number of points in one half of the ellipse
alpha = 0.05;  % 1-alpha confidence

switch(lower(action))

 case 'init'

  ha  = varargin{1};

  x1 = repmat(NaN,1,hnPoint*2+1);
  x2 = repmat(NaN,1,hnPoint*2+1);

  hl  = line(x1,x2,'color','cyan');

  varargout{1} = hl;

 case 'update'

  hl  = varargin{1};
  mu  = varargin{2};
  sig = varargin{3};
  rho = varargin{4};
  
  r2 = spm_invXcdf(1-alpha,2);
  
  z2 =linspace(-sqrt(r2),sqrt(r2),hnPoint);
  
  z1a = rho * z2 + sqrt(1-rho^2)*sqrt(r2-z2.^2);
  z1b = rho * z2 - sqrt(1-rho^2)*sqrt(r2-z2.^2);
  
  % Now de-normalize (i.e. shift and scale)
  x1a = mu(1) + z1a*sig(1);
  x1b = mu(1) + z1b*sig(1);
  x2  = mu(2) + z2 *sig(2);
  
  set(hl,'Xdata',[x1a,NaN,x1b],'Ydata',[x2,NaN,x2])
  
end

return


function PlotHere(ha)

f = get(ha,'Parent');
set(0,'CurrentFigure',f);
set(f,'CurrentAx',ha);

return
