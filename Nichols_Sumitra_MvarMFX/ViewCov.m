function varargout = ViewCov(action,varargin)
%
% ViewCov('init','Corr.img',n,{labels})
%    Initialize
%
% ViewCov('setloc',mm[,fig])   ViewCov('setmm',mm[,fig])
%    Set location in mm
%
% ViewCov('update'[,f])
%    Update; only useful if VC structure has been modified manually.
%
% ViewCov('SetCoords',mm,h,hReg)
%    See spm_XYZreg; used for interacting with spm_orthviews
%
%
% $Id: ViewCov.m,v 1.1 2005/12/13 20:34:33 nichols Exp $

if (nargin<1), action = 'init'; end

switch lower(action)

 case 'init'
 %------------------------------------------------------------------------
  if (nargin >= 2), 
    Cor = varargin{1};
  else
    Cor = spm_get(1,'IMAGE','Select Cov Image'); 
  end
  if (nargin >= 3)
    Lab = varargin{2};
  else
    Lab = GetLabels(p_from_k(Nvol(Cor)));
  end
  if (nargin >= 4)
    n = varargin{3};
  else
    n = spm_input('Number of subj?','+1','n');
  end

  f = CreateVC(Cor,Lab,n);
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



function f = CreateVC(Cor,Lab,n)
% Cor - Correlation image filename
% Lab - Length p cell array of condition labels
% n   - Number of subjects

% Preliminaries
%   Find p, map correlation elements
%-------------------------------------------------------------------------

k       = Nvol(Cor);
p       = p_from_k(k);
Corr    = eye(p);
CorrP   = repmat(NaN,p,p);
Stdev   = ones(1,p);
Q       = BuildQ(p);

if ~isempty(findstr(Cor,','))
  Cor(min(find(Cor==',')):end) = [];
end

for i=1:k
  Vcor(i) = spm_vol(sprintf('%s, %d',Cor,i));
end

% Create default/empty Figure
%-------------------------------------------------------------------------

f = figure('position',[100 100 350 300],...
	   'MenuBar','none',...
	   'Tag','ViewCov');
colormap(winter);

% Set button
hui = uicontrol(f,'String','',...
		'Style','PushButton',...
		'Position',[10 7 100 20 ],...
		'CallBack','ViewCov(''ToggleReg'',gcbf)',...
		'Visible','on');
	  

% Image
hi = imagesc(zeros(p,p),[-1 1]);
ha = gca;
axis image
set(hi,'ButtonDownFcn',['VC=get(gcbf,''UserData'');'...
	    'Stdev=VC.Stdev,Cor=VC.Corr,CorP=VC.CorrP']);

set(ha,'Position',[.15 .21 0.65 0.65],...
       'Xtick',1:p,'Ytick',1:p,...
       'XtickLabel',Lab,'YtickLabel',Lab);
hc = colorbar;
set(hc,'Position',[.85 .21 0.06 0.65]);

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
ht_loclab = text(.4,.01,...
	      {'mm:','vox:'},...
	      'VerticalAlign','Bottom');
ht_loc = text(.5,.01,...
	      {sprintf('%g %g %g',[NaN NaN NaN]),...
	       sprintf('%d %d %d',[NaN NaN NaN])},...
	      'VerticalAlign','Bottom');

% Status Text: LR value
X2 = NaN;
X2p = NaN;
ht_x2lab  = text(.03,.96,'LRT \chi^2 =');
ht_x2     = text(.15,.96,sprintf('%g',X2));
ht_x2plab  = text(.03,.90,'LRT P  =');
ht_x2p     = text(.15,.90,sprintf('%g',X2p));

% Legend 
LegAtt = {'fontweight','bold','Color'};
Leg(1) = text(.30,.97,'Correlation',LegAtt{:},'White');
Leg(2) = text(.30,.93,'Standard Deviation',LegAtt{:},'Black');
Leg(3) = text(.30,.89,'P-value for Corr',LegAtt{:},'Yellow');
Leg(4) = text(.6,.89,sprintf('Bonf Thresh = %f',0.05/(p*(p-1)/2)));

h  = struct('i',hi,...            % Image
	    'a',ha,...            % Image axis
	    'c',hc,...            % Colorbar
	    'ui',hui,...            % Colorbar
	    'tv',hvt,...          % Image value text
	    'ta',ha_text,...      % Whole-figure axis, for text
	    'tll',ht_loclab,...   % mm/vox label text
	    'tl',ht_loc,...       % mm/vox value text
	    'tx2l',ht_x2lab,...   % LR label text
	    'tx2',ht_x2,...        % LR stat text
	    'tx2pl',ht_x2plab,... % LR P label text
	    'tx2p',ht_x2p,...      % LR P stat text
	    'leg',Leg...          % Misc legend text
	    );

VC = struct('p',p,...
	    'n',n,...
	    'mat',Vcor(1).mat,...
	    'Vcor',Vcor,...
	    'Stdev',Stdev,...
	    'Corr',Corr,...
	    'CorrSig',CorrP,...
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
	get(hReg);
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
  spm_XYZreg('Add2Reg',hReg,hMe,'ViewCov');
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
hvt   = VC.h.tv;

%
% Read stdev, corr
%
if any(isnan(vox))
  Corr = eye(p);
  Stdev = ones(1,p);
else
  Stdev = spm_get_data(Vcor(1:p),vox')';
  CorDat  = spm_get_data(Vcor(p+1:end),vox');

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

% tor changed this so diagonals will be color-mapped same as others
% doesn't seem to change LRT values
set(VC.h.i,'CData',  Corr); %%CombCorStd(Corr,Stdev));

% Display image values
tmp = tril(Corr)+diag(Stdev)+triu(CorrP,+1);
% tor changed line above because we were getting Cor values of 2 on
% diagonal
tmp = tril(Corr,-1)+diag(Stdev)+triu(CorrP,+1);

for x=1:p
  for y=1:p
    set(hvt(x,y),'String',sprintf('%.3f',tmp(x,y)));
  end
end

%
% Compute LR
%

VC.X2 = -(n-1-(2*p+5)/6)*log(det(Corr));
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

fp   = fopen(V.fname);
fseek(fp,0,'eof');
Len  = ftell(fp);
fclose(fp);
n    = Len/(prod(V.dim(1:3))*spm_type(V.dim(4),'bits')/8);

return



function p = p_from_k(k)
% Computes p, dimension of covariance matrix, from k, number of unique
% elements. 

p = (-1+sqrt(1+8*k))/2;

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
% Check that f is a ViewCov figure, and optionally crash out
%
% If f is empty, try to return a valid figure handle
%
if nargin<2, crash = 0; end

fs = findobj('Tag','ViewCov');


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
  error('Selected figure is not a ViewCov figure!')
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
