function CovEst(VY,n)
% CovEst(VY,n)
%
% VY    - a n*p vector of mapped image files
% n     - Number of subjects
%
% Output
%   Lamda.img  Contains p*(p+1)/2 volumes, the coeffieicients of
%              the covariance structure
%   Corr.img   Contains p*(p+1)/2 volumes; first the p standard
%              deviations and then the p*(p-1)/2 correlation
%              coeffieicients.
%   CovEst.mat  Useful variables, include Q
%
% $Id: CovEst.m,v 1.1 2005/12/13 18:56:23 nichols Exp $

if ~isstruct(VY), VY = spm_vol(VY); end
p    = length(VY)/n;
Dim = VY(1).dim(1:3);


%
% Prep work
%
[Q,Qp,B,R] = BuildQ(n,p);
nQ     = length(Q);
Vlam   = MakeLam(VY(1),nQ);
Vcor   = MakeCor(VY(1),nQ);

Ctr    = eye(n*p)-R*pinv(R);

for z = 1:Dim(3)
  fprintf('Slice %2d: reading',z);
  Y = zeros(n*p,prod(Dim(1:2)));
  for i=1:n*p
    tmp      = spm_slice_vol(VY(i),spm_matrix([0 0 z]),Dim(1:2),1);
    Y(i,:)   = tmp(:)';
  end
  Msk = all(~isnan(Y));
  Y = Y(:,Msk);
  nvox = size(Y,2);

  fprintf('. Computing %d vox',nvox);

  % Center Y
  Yc = Ctr*Y;

  % Compute Lambda coefficients (covariance elements)
  Lam = zeros(nQ,nvox);
  for i=1:nQ
    Lam(i,:) = 1/(n-1) * sum((Qp{i,1}'*Yc).*(Qp{i,2}'*Yc));
  end
  % Compute correlation coefficients
  Cor = zeros(nQ-p,nvox);
  for i=1:nQ
    if i<=p
      Cor(i,:) = sqrt(Lam(i,:));
    else
      A  = full(Q{i}(1:p,1:p));
      ii = find(max(A,[],1));
      Cor(i,:) = Lam(i,:)./(sqrt(Lam(ii(1),:).*Lam(ii(2),:)));
    end
  end


  fprintf('.  Writing');
  % Write lambda
  for i=1:nQ
    WriteS(Lam(i,:), Msk, Vlam(i), z);
    WriteS(Cor(i,:), Msk, Vcor(i), z);
  end
  fprintf('.\n');
  
end

% Close lambda
spm_close_vol(Vlam);
spm_close_vol(Vcor);

clear Y Yc Msk nvox Lam
save CovEst.mat 
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

function Vlam = MakeLam(V,nQ)
V.dim(4) = spm_type('double');
for i = 1:nQ
  Vlam(i) = MkOutImg(V,'Covar.img','Unique var-covar elements',i);
end
return

function Vcor = MakeCor(V,p)
V.dim(4) = spm_type('double');
for i = 1:p
  Vcor(i) = MkOutImg(V,'Corr.img','Correlation coeffs',i);
end
return

function V = MkOutImg(Vtemplate, fname, descrip,n)
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
