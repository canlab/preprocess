 function [MAOV1] = MAOV1(X,alpha)
%Single-Factor Multivariate Analysis of Variances test.
%Computes a Multivariate Analysis of Variance for equal or unequal sample sizes.
%[Testing of the mean differences in several variables among several samples (groups).
%It consider the matrices of sum-of-squares between (H) and within (E) to compute the
%Wilks'L (lambda) and according to the number of data it chooses to approximate to the
%Chi-square distribution or through several multiple F-approximations such as Rao,
%Pillai, Lawley-Hotelling and Roy tests.]
%
%   Syntax: function [MAOV1] = MAOV1(X,alpha) 
%      
%     Inputs:
%          X - data matrix (Size of matrix must be n-by-(1+p); sample=column 1, variables=column 2:p). 
%      alpha - significance (default = 0.05).
%     Outputs:
%          - Complete multivariate analysis of variance table.
%
%    Example: From the example 6.8 of Johnson and Wichern (1992, pp. 249-251), to test the 
%             multivariate analysis of variance of data with a significance of 0.05.
%
%                               Data
%                      ----------------------
%                       Sample    X1      X2
%                      ----------------------
%                          1       9       3
%                          1       6       2
%                          1       9       7
%                          2       0       4
%                          2       2       0
%                          3       3       8
%                          3       1       9
%                          3       2       7
%                      ----------------------
%
%     Data matrix must be:
%             X=[1 9 3;1 6 2;1 9 7;2 0 4;2 2 0;3 3 8;3 1 9;3 2 7];
%
%     Calling on Matlab the function: 
%             maov1(X)
%
%     Answer is:
%
%   It is considering as a small sample problem (n < 25).
%
%   Multivariate Analysis of Variance Table.
%   --------------------------------------------
%   No. data    Samples     Variables       L
%   --------------------------------------------
%       8          3            2        0.0385
%   --------------------------------------------
% 
%   ------------------------------------------------------------------------------
%   Test                 Statistic     df1     df2         F       P    Conclusion
%   ------------------------------------------------------------------------------
%   Rao                    0.038         4       8       8.20   0.0062       S
%   Pillai                 1.541         4      10       8.39   0.0031       S
%   Lawley-Hotelling       9.941         4       6       7.46   0.0164       S
%   Roy                    8.076       2.0     5.0      20.19   0.0040       S
%   ------------------------------------------------------------------------------
%   With a given significance of: 0.05
%   According to the P-value, the sample mean vectors could be significant (S) or
%   not significant (NS).
%

%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
%  August 10, 2003.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2003). MAOV1: Single-Factor Multivariate Anaysis of 
%    Variance test. A MATLAB file. [WWW document]. URL http://www.mathworks.com/matlabcentral/fileexchange/
%    loadFile.do?objectId=3863&objectType=FILE
%
%  References:
% 
%  Boik, J. R. (2002). Lecture Notes: Statistics 537. Classical Multivariate Analysis.
%              Spring 2002. pp. 38-41. [WWW document]. URL http://math.montana.edu/~rjboik/
%              classes/537/notes_537.pdf
%  Johnson, D. E. (1998), Applied Multivariate Methods for Data Analysis.
%              New-York:Brooks Cole Publishing Co., an ITP Company. Chapter 11.
%  Johnson, R. A. and Wichern, D. W. (1992), Applied Multivariate Statistical Analysis.
%              3rd. ed. New-Jersey:Prentice Hall. pp. 246-251.
%  SAS/Stat User's Guide. Chap. 3. Introduction to Regression Procedures: Statistical
%              Background: Multivariate Tests. Software-Kategorie Statistik und Mathematik:
%              Datenanalyse mit SAS mit SAS Online Doc. Universität Zürich. [WWW document].
%              URL http://www.id.unizh.ch/software/unix/statmath/sas/sasdoc/stat/chap3/sect10.htm
%

P = [];PP = [];PLH = [];PR = [];
if nargin < 2,
   alpha = 0.05;
end; 

g = max(X(:,1));

N = [];
indice = X(:,1);
for i = 1:g
   Xe = find(indice==i);
   eval(['X' num2str(i) '= X(Xe,2);']);
   eval(['n' num2str(i) '= length(X' num2str(i) ') ;'])
   eval(['xn= n' num2str(i) ';'])
   N = [N,xn];
end;

[f,c] = size(X);
X = X(:,2:c);

[n,p] = size(X);
r = 1;
r2 = N(1);

for k = 1:g
   eval(['M' num2str(k) '= mean(X(r:r2,:));']);
   if k < g
      r = r+N(k);
      r2 = r2+N(k+1);
   end;
end;
M = mean(X);

dT = [];
for k = 1:p
   eval(['dT  = [dT,(X(:,k) - mean(X(:,k)))];']);
end;
T = dT'*dT;  %total sum of squares

r = 1;
r2 = N(1);
g = length(N);
for k = 1:g
   Md(k,:) = (mean(X(r:r2,:)) - mean(X));
   if k < g
      r = r+N(k);
      r2 = r2+N(k+1);
   end;
end;

H = [];  %between samples sum of squares
for k = 1:g
   h = N(k)*(Md(k,:)'*Md(k,:));
   if k == 1
      H = h;
   else
      H = H+h;
   end;
end;

E = T-H;  %within samples sum of squares
% Wilks'test (Wilks'U)
LW = det(E)/det(T);  %Wilks'lambda

if f >= 25;  %approximation to chi-square statistic
   disp('It is considering as a large sample problem (n >= 25).')
   v = p*(g-1);
   X2 = (-1)*(((f-1)-(.5*(p+g)))*log(LW));
   P = 1-chi2cdf(X2,v);
   disp(' ')
   disp('Multivariate Analysis of Variance Table.')
   fprintf('-----------------------------------------------------------------------------------\n');
   disp('No. data    Samples     Variables       L          Chi-sqr.         df          P')
   fprintf('-----------------------------------------------------------------------------------\n');
   fprintf('%5.i%11.i%13.i%14.4f%15.4f%12.i%13.4f\n',f,g,p,LW,X2,v,P);
   fprintf('-----------------------------------------------------------------------------------\n');
   fprintf('With a given significance of: %.2f\n', alpha);
   if P >= alpha;
      disp('Sample mean vectors results not significant.');
   else
      disp('Sample mean vectors results significant.');
   end;
else  %approximation to Rao's F-statistic
   disp('It is considering as a small sample problem (n < 25).')
   if p == 2 | (g-1) == 2
      s = 2;
   else
      s = sqrt((p^2*(g-1)^2-4)/(p^2+((g-1)^2)-5));
   end;
   v1 = p*(g-1);
   v2 = (s*((f-1)-((p+g)/2)))-(((p*(g-1))-2)/2);
   v = LW^(1/s);
   F = ((1-v)/v)*(v2/v1);
   P = 1-fcdf(F,v1,v2);
   
   %Pillai's test (Trace test) 
   V = trace(inv(T)*H);
   q = g-1;
   s = min(p,q);
   v = f-g;
   m = (abs(p-q)-1)/2;
   n = (v-p-1)/2;
   FP = (((2*n)+s+1)/((2*m)+s+1))*(V/(s-V));
   v1P = s*((2*m)+s+1);
   v2P = s*((2*n)+s+1);
   PP = 1-fcdf(FP,v1P,v2P);
   
   %Lawley-Hotelling's trace test
   U = trace(inv(E)*H);
   FLH = (2*((s*n)+1)*U)/(s^2*((2*m)+s+1));
   v1LH = s*((2*m)+s+1);
   v2LH = 2*((s*n)+1);
   PLH = 1-fcdf(FLH,v1LH,v2LH);
   
   %Roy's Union_intersection test (Largest root test) 
   R = max(eig(inv(E)*H));
   r = max(p,q);
   FR = R*((v-r+q)/r);
   v1R = r;
   v2R = v-r+q;
   % Because the numerator and denominator degrees of freedom could results a fraction,
   % the probability function associated to the F statistic is resolved by the Simpson's
   % 1/3 numerical integration method.
   x = linspace(.000001,FR,100001);
   DF = x(2)-x(1);
   y = ((v1R/v2R)^(.5*v1R)/(beta((.5*v1R),(.5*v2R))));
   y = y*(x.^((.5*v1R)-1)).*(((x.*(v1R/v2R))+1).^(-.5*(v1R+v2R)));
   N2 = length(x);
   PR = 1-(DF.*(y(1)+y(N2) + 4*sum(y(2:2:N2-1))+2*sum(y(3:2:N2-2)))/3.0);
   
   if P >= alpha;
      ds1 ='NS';
   else
      ds1 =' S';
   end;
   if  PP >= alpha;
      ds2 ='NS';
   else
      ds2 =' S';
   end;
   if  PLH >= alpha;
      ds3 ='NS';
   else
      ds3 =' S';
   end;
   if PR >= alpha;
      ds4 ='NS';
   else
      ds4 =' S';
   end;
   ;
   disp(' ')
   disp('Multivariate Analysis of Variance Table.')
   fprintf('--------------------------------------------\n');
   disp('No. data    Samples     Variables       L')
   fprintf('--------------------------------------------\n');
   fprintf('%5.i%11.i%13.i%14.4f\n',f,g,p,LW);
   fprintf('--------------------------------------------\n');
   disp(' ')
   fprintf('------------------------------------------------------------------------------\n');
   disp('Test                 Statistic     df1     df2         F       P    Conclusion')
   fprintf('------------------------------------------------------------------------------\n');
   fprintf('Rao            %13.3f%10i%8i%11.2f%9.4f      %s\n',LW,v1,v2,F,P,ds1);
   fprintf('Pillai         %13.3f%10i%8i%11.2f%9.4f      %s\n',V,v1P,v2P,FP,PP,ds2);
   fprintf('Lawley-Hotelling      %6.3f%10i%8i%11.2f%9.4f      %s\n',U,v1LH,v2LH,FLH,PLH,ds3);
   fprintf('Roy            %13.3f%10.1f%8.1f%11.2f%9.4f      %s\n',R,v1R,v2R,FR,PR,ds4);
   fprintf('------------------------------------------------------------------------------\n');
   fprintf('With a given significance of:% 3.2f\n', alpha);
   disp('According to the P-value, the sample mean vectors could be significant (S) or'); 
   disp('not significant (NS).');
end;

