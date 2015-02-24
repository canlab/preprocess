function CovLRT(VY,n,Vcov)
% CovLRT(VY,n)
%
% VY    - a n*p vector of mapped image files
% n     - Number of subjects
% Vcov  - Mapped covariance matrix
%
% Output
%   LRT.img    LRT of Ho vs Ha
%   LRT_P.img  -log10 P-value assuming asymtotic \chi^2 dist of LRT
%   F.img      F test for any effect
%   F_P.img    -log10 P-value for F test
%   CovLRT.mat  
%
%
% Uses Bartlett's test for correlation
%   Bartlett, M.S.
%   A Note on the Multiplying Factors for Various \ch^2 Approximations.
%   Journal of the Royal Statistical Society. Series B (Methodological),
%   Vol. 16, No. 2, 1954, pp. 296-298.
%   
%
% $Id: CovLRT.m,v 1.3 2005/12/13 20:34:44 nichols Exp $

if ~isstruct(VY), VY = spm_vol(VY); end
if ~isstruct(Vcov), Vcov = spm_vol(Vcov); end
p   = length(VY)/n;
nQ  = p*(p+1)/2;
Dim = VY(1).dim(1:3);
if nQ~=length(Vcov), 
  error(sprintf('Wrong number of covariance elements nQ=%d Vcov:%d',...
		nQ,length(Vcov)))
end

%
% Prep work
%
Vlrt   = MkOutImg(VY(1),'LRT.img','LRT of Corr');
Vlrtp  = MkOutImg(VY(1),'LRT_nlP.img','-log10 P LRT of Corr');
VF     = MkOutImg(VY(1),'F.img','LRT of Corr');
VFp    = MkOutImg(VY(1),'F_nlP.img','-log10 P LRT of Corr');

X     = spm_DesMtx(repmat(1:p,1,n)');
Xpi   = pinv(X);
Ctr   = eye(n*p)-X*Xpi;

[Qm,Q1m] = MakeQmat(n,p);

for z = 1:Dim(3)
  fprintf('Slice %2d: reading',z);
  [Y,Cov,Msk,nvox] = ReadSlice(VY,Vcov,Dim,z,n*p,nQ);
  
  fprintf('. Computing %d vox',nvox);

  LRT  = zeros(1,nvox);  LRTp = zeros(1,nvox);
  F    = zeros(1,nvox);  Fp   = zeros(1,nvox);
  for i=1:nvox

    Covi = reshape(Q1m*Cov(:,i),p,p);
    Vari = diag(Covi); 
    Cori = diag(sqrt(1./Vari))*Covi*diag(sqrt(1./Vari));
    
    % Bartlett's test
    LRT(i) = -(n-1 -(2*p+5)/6)*log(det(Cori));
    LRTp(i) = -log10(max(eps,1-spm_Xcdf(LRT(i),p*(p-1)/2)));

    % Do F-test for any effect at all
    betah = Xpi*Y(:,i);
    F(i) = n*betah'*inv(Covi)*betah/p;
    Fp(i) = -log10(max(eps,1-spm_Fcdf(F(i),p,n*p-p)));

  end
  
  fprintf(' (Mx=%g Mn=%g).  Writing',max(LRT),min(LRT));

  % Write LRT
  WriteS(LRT, Msk, Vlrt, z);  WriteS(LRTp, Msk, Vlrtp, z);
  WriteS(F,   Msk, VF,   z);  WriteS(Fp,   Msk, VFp,   z);

  fprintf('.\n');
  
end

% Close lambda
spm_close_vol(Vlrt);  spm_close_vol(Vlrtp);
spm_close_vol(VF);    spm_close_vol(VFp);

clear SIG iSIG SIG0 iSIG0 LRT LRTp F Fp
save CovLRT.mat 
return



%%%
%%%  Support functions
%%%
%%%

function [Q,Qp,B,R] = BuildQ(n,p);

I = repmat(1:n,p,1);
I = [I(:) repmat(1:p,1,n)' zeros(n*p,2)];
var = [0 1 0 0];
dep = [0 1 0 0];
xVi = struct('I',I,'var',var,'dep',dep);
xVi = spm_non_sphericity(xVi);

Q    = xVi.Vi;
nQ   = length(Q);
Qp   = cell(nQ,2);
Blk  = spm_DesMtx(I(:,1));
Rep  = spm_DesMtx(I(:,2));
for i = 1:nQ
  Qp{i,1} = triu(Q{i})'*Blk;
  Qp{i,2} = triu(Q{i})*Blk;
end
B = Blk;
R = Rep;
return

function V = MkOutImg(Vtemplate, fname, descrip,n)
if nargin<4, n=1; end
Type = 'double';

V         = Vtemplate;
V.dim(4)  = spm_type(Type);
V.fname   = fname;
V.descrip = descrip;
V.n       = n;
spm_create_vol(V);
V.pinfo(3)=prod(V.dim(1:3))*spm_type(Type,'bits')/8*(n-1);

return

function V = WriteS(Data, Msk, V, z)
Sl = repmat(NaN, V.dim(1:2));
Sl(Msk) = Data;
spm_write_plane(V,Sl,z);

return

function [Qm,Q1m] = MakeQmat(n,p)
% Matricies to help 'unwrap' the covariance elements into covariance
% matrices 

Q = BuildQ(n,p);
nQ = length(Q);
Qm = zeros((n*p)^2,nQ);
Q1m = zeros(p^2,nQ);
for i=1:nQ
  tmp     = Q{i};
  Qm(:,i) = tmp(:);
  tmp     = tmp(1:p,1:p);
  Q1m(:,i)= tmp(:);
  if i==p
    Qom = Qm(:,1:p);
    Qo1m = Q1m(:,1:p);
  end
end

%%Note used
%% % These matricies 'estimate' the lambda values for a particular
%% % covariance matrix
%%Q1e  = inv(Q1m'*Q1m)*Q1m';
%%Qo1e = inv(Qo1m'*Qo1m)*Qo1m';

return

function [Y,Cov,Msk,nvox] = ReadSlice(VY,Vcov,Dim,z,nY,nC)

Y = zeros(nY,prod(Dim(1:2)));
Cov = zeros(nC,prod(Dim(1:2)));
for i=1:nY
  tmp      = spm_slice_vol(VY(i),spm_matrix([0 0 z]),Dim(1:2),1);
  Y(i,:)   = tmp(:)';
end
for i=1:nC
  tmp      = spm_slice_vol(Vcov(i),spm_matrix([0 0 z]),Dim(1:2),1);
  Cov(i,:) = tmp(:)';
end
Msk = all(~isnan(Y));
Y    = Y(:,Msk);
Cov  = Cov(:,Msk);
nvox = size(Y,2);

return
  