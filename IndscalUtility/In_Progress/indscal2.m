function [X,W,sp,Brecon] = indscal(varargin)
% indscal(dat,r,names)
%
% Carroll & Chang's INDSCAL algorithm for 3-way Metric Weighted MDS (see pdf with same name)
% -----------------------------------------------
% By Howard Wang and Tor Wager
%
% Function call with no arguments gets sample data from U.S. Cities
%

% Programmers' notes
% 1 - Loop way to get scalar product matrices Bk
% for k=1:m
% 	Di(:,:,k) = (1 - tmp(:,:,k)) .^ .5;
% 	Ok = Di(:,:,k).^2;   % the subject k's matrix of squared dissimilarity data (n x n)
% 
% 	Bk = Ok-Ok*U*inv(U'*U)-inv(U'*U)*U'*Ok+inv(U'*U)*U'*Ok*U*inv(U'*U);    %this is the double centered scalar-product matrix for subejce k (n x n)
% 
% 	for i = 0:n-1
% 		B1(k,i*n:(i+1)*n) = Bk(n,:);    %  re-arrange matrices Bk into an (m x n^2) matrix B1 whose k'th row is formed from Bk by horizontally adjoining rows of Bk. 
% 	end
% 	B2((k-1)*n:k*n,:)=Bk(:,:);     % re-arrange matrices Bk into an (m*n x n) matrix B2 by vertically adjoining the m matrices Bk. 
% end
% 

% -----------------------------------------------
% Input arguments
% -----------------------------------------------
%names = []; 
%r = 2;
mnm=[];
if length(varargin) == 0    
    error('sample data no longer works')
    % load sample data
    %[D,names,dat,weightsactual,sp,r] = sample_data;,
    % D should be dissimilarity
else 
    dat = varargin{1};
    d=dat(:,:,1);
    del = 10*eps;  
    m=size(dat,3);
    % this is where we convert from sim to diss, if necessary
    if all(abs(diag(d) - 1) < del) & all(all(d < 1+del))        
        warning('converting similarity to dissimilarity')        
        for k=1:m
        d=dat(:,:,k);                    
        d = sqrt(1 - d);
        dat(:,:,k)=d;
        clear d;
        end    
    end    
end

if dat(1,1,1)~=0
    error('not dissimilarity matrix');
end

if length(varargin) > 1, r = varargin{2}; end
if length(varargin) > 2, names = varargin{3}; end
if length(varargin) > 3, cmat = varargin{4}; end
if length(varargin) > 4, mnm = varargin{5}; end

% -----------------------------------------------
% Input to test the program
% -----------------------------------------------
if isempty(r), disp('defaulting to 2d'),r = 2;, end
if isempty(mnm), disp('defaulting to metric'),mnm={'m'};, end
    
if isempty(dat),  end
  
% -----------------------------------------------
% Initial variables
% -----------------------------------------------

% tmp = ...  this is the correlation between subjects, this is the input from Matrix_Multivariate (n x n x m)
% r = 2; %  the number of relavant dimensions in the stimulus space, input from Matrix_Multivariate
% D = ...   distance matrix between the stimuli. (n x n)

%iteration = 100; % the number of iterations



% -----------------------------------------------
% Step 1 - Setting up scalar-product matrices (double centering)
% -----------------------------------------------

    % GET DISTANCE MATRIX
    % IF dat = distance, DD = dat
    % We'll replace the next 2 lines with code like in cmdscale that checks
    % if input is sim or dist
    %Di = tmp(:,:,k) .^ .5;     % convert correlation to distance
    
  % convert to scalar product
  [sp,DD,B1,B2] = dist2scalprod(dat); % sp is scal prod, DD is 3-D distances


% -----------------------------------------------
% Step 2, 3 - find X, the stimulus matrix (n x r)
% -----------------------------------------------

% D = mean(Di,3);
% X=cmdscale(D); X=X(:,1:r)    	%    X is the first 2 columns of cmdscale of D
%X=mean(sp,3);

% Get group stimulus coordinates
Dgrp = mean(DD,3);   % group average distance matrix
if char(mnm)=='m';
    X=cmdscale(Dgrp);
    X = X(:,1:r);     % Save only specified number of dimensions
else
    X = cnmds(Dgrp,r);  % Group average stimulus coordinates  
end

% -----------------------------------------------
% Step 4 - iterations to minimize Bk- XL*Wk*XR
% -----------------------------------------------


%if ndims(X)~=1;
%fn=nmdsfig(X,ones(size(X,1),1),names); % Make a figure to check
%title('group space');
%end

[fit,X,W,XL,XR] = indscalf(X,sp,B1,B2);

%W = W .^ 2;     % interpret squared weights ????, but do not square them now!

colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko'};
% plot final weights
%fwei = tor_fig; 

nmdsfig(W,cmat,cell(1,length(cmat)));
title('subject weights');
for c=1:max(cmat);
leg{c}=['cond',num2str(c)];
end
legend(leg);

%axes([0 1],[0 1]);
%axis([0 max(W(:,1)) 0 max(W(:,2))])

%%NORMALIZATION OF X AND W
%nX=[zscore(X(:,1)) zscore(X(:,2))];     %make variances equal to 1 for each dimension
%cX=nX(1,:)./X(1,:);   %constant scaling factor 
%cW=1./(cX.^2);
%nW=[W(:,1).*cW(1) W(:,2).*cW(2)];
%nmdsfig(nX,ones(size(nX,1),1),names); % Make a figure to check
%title('normalized group space');

%figure;
%for n=1:max(cmat);
%hold on; plot(W(find(cmat==n),1)./std(W(:,1)), W(find(cmat==n),2)./std(W(:,2)),colors{n})
%hold on; plot(nW(find(cmat==n),1), nW(find(cmat==n),2),colors{n});
%    for s=1:size(nW,1)
%    text(nW(s,1)+0.00002,nW(s,2),num2str(s));
%    end    
%end
%title('normalized subject weights');

for ii = 1:size(W,1), Brecon(:,:,ii) = XL * diag(W(ii,:)) * XR';, end
fit = (sp - Brecon).^2; fit = sqrt(sum(fit(:))./prod(size(Brecon)));
fprintf(1,'Residual error (rms) is %3.4f, fitness logit score (0 - 1) is %3.3f\n',fit,1./(1+fit))
return

% -----------------------------------------------
% Get Scalar Products
% -----------------------------------------------

function [sp,DD,B1,B2] = dist2scalprod(tmp)   
%input dissimilarity, output scalar products


n = size(tmp,1); % this the number of stimuli
m = size(tmp,3); % this the number of subjects
U = repmat(1,n,1);  % a column of n ones

P1 = eye(n) - repmat(1/n,n,n);  % dbl-centered identity matrix
B2 = [];
for k = 1:m     % loop through m subjects
    
    Di = tmp(:,:,k);            % already a dist matrix - n x n
    
    DD(:,:,k) = Di;             % save 3-D Distance Matrix
    Bk = P1 * (-.5 * Di .* Di) * P1;    % P1 mult. reverses the dot mult. of Di .* Di
    % SAME AS THIS - double-centers subject distances and squares each one.
    % Bkk = scale(Di,1); Bkk = scale(Bkk',1)'; Bkk = Bkk .* Bkk;
    
    Bk = (Bk+Bk')./2;       % minimize roundoff errors

    % This is where we normalize Bk to make its SS == 1, if desired
    % convert to ss=1
    ssBk=sum(Bk(:).^2);
    Bk=sqrt((Bk.^2)./ssBk);
           
    B1(k,:) = Bk(:)';
    %for i = 0:n-1
	%	B1(k,i*n+1:(i+1)*n) = Bk(n,:);    %  re-arrange matrices Bk into an (m x n^2) matrix B1 whose k'th row is formed from Bk by horizontally adjoining rows of Bk. 
    %end
	%B2((k-1)*n+1:k*n,:)=Bk; 
    B2 = [B2; Bk];
    sp(:,:,k)=Bk;
end

return





% -----------------------------------------------
% Alternating least squares - core of indscal
% -----------------------------------------------

function [fit,X,W,XL,XR] = indscalf(X,sp,B1,B2)
% [fit,X,W] = indscalf(X,dp,B1,B2)
%
% 

m=size(sp,3);

XR=X; XL=X;    			%    set XL and XR to be X
iteration=100;
r=size(X,2)

for i = 1:iteration
	
	% part 1 - use estimates of XL to estimate W
    % tmp = XL(:,1) * XL(:,1)';
    
    % B = WP is the (transposed) regression equation
    % each data vector is a row of B
    % B = P*w, for each subject, where w is a vector of weights (col. of W)
    % so P is like the model space, design space.
    % XL is latest estimate of stimulus coordinates in r dimensions
    % each row of P is crossproduct of all stim coordinates XL in dimension r
    % which is estimated 'scalar product'.  
    for j = 1:r
        tmp = XL(:,j) * XR(:,j)';   % outer product of dimension j
        P(j,:) = tmp(:)';           % strung out!  reproduction of data (B) ??
    end
	% _______________________________________
	%for j = 1:r
	%	for k = 0:n-1
	%		for o = 1:n
	%			P(j,k*n+o) = XL(k+1,j)*XR(o,j); 
                %		end
                %	end
                %end

    % Fitting model space P against data B1 to get parameter estimates W     
    % B1 = PW   solve for W
    % Or, sp = XL * W * XR, solve for W
    % W2 = inv(XL'*XL)*XL'*sp(:,:,1)*XR*inv(XR'*XR);    diagonals of this
    % are weights.
    % for j = 1:size(sp,3),tmp=inv(XL'*XL)*XL'*sp(:,:,j)*XR*inv(XR'*XR);,W2(j,:)=diag(tmp)';,end
    
	W=B1*P'*inv(P*P');		% m x r matrix, rows are subjects, columns are dimension weights
    rW(i,:)=ranger(W);
    % add noise
    %W = W + randn(size(W)) * .10*mean(mean(W)); W(W < 0) = 0;
    
	for k=1:m			% set of r x r diag matrices with weights on diagonal
		Wk(:,:,k) = diag(W(k,:));
	end
		

	
	% part 2 - use estimates of XL and W to estimate XR
	% _______________________________________
		
	% L is matrix of weighted XL, concatenated across subjects (rows)	
	L=[];
    for k=1:m
		L=[L; XL*Wk(:,:,k)];     % re-arrange matrices XL into an (m*n x n) matrix L by vertically adjoining the m matrices XL*Wk.
	end
	%XR=(pinv(L)*B2)'; % finding XR from XL
    XR = (inv(L' * L) * L' * B2)';
    
 	% part 3 - use estimates of XR and W to estimate XL
	% _______________________________________

	R=[];
    for k=1:m
		R=[R; XR*Wk(:,:,k)];     % re-arrange matrices XR into an (m*n x n) matrix R by vertically adjoining the m matrices XR*Wk.
    end
	%XL=B2'*pinv(R'); % finding XL from XR
    XL = B2' * R * inv(R' * R);

    % OR...solve XL*W*XR = B2' for XL, comes out the same
    % this seems more intuitive.  
    %R=[];
    %for k=1:m
    %    R=[R Wk(:,:,k) * XR']; 
    %end
    %tmp = B2' * R' * inv(R * R');
end

for ii = 1:size(W,1), Brecon(:,:,ii) = XL * diag(W(ii,:)) * XR';, end
fit = (sp - Brecon).^2; fit = sqrt(sum(fit(:))./prod(size(Brecon)));
fit1 = fit;
fit = 1./(1+fit);   % btwn 1 and 0 for positive error values, 1 = perfect fit (max fitness)
%fprintf(1,'Residual error (rms) is %3.3f, fitness logit score (0 - 1) is %3.3f\n',fit1,fit)

return

%WW = W'; WW = WW(:); WW = diag(WW);     % make super-subject weight matrix
%XL = repmat(XL,size(W,1));            % repeat all XL
%XR = repmat(XR,size(W,1));            % repeat all XL
%B = XL * WW * XR';
% doesn't work b.c 0 weights no invertible







% -----------------------------------------------
% Sample data - U.S. Cities
% -----------------------------------------------

function [D,names,dat,weightsactual,sp,r] = sample_data

  r = 2;
  
names = {'Birmingham, Ala.'
'Boston, Mass.'
'Buffalo, N. Y.'
'Chicago, Ill.'
'Cleveland, Ohio'
'Dallas, Tex.'
'Denver, Colo.'
'Detroit, Mich.'
'El Paso, Tex.'
'Houston, Tex.'
'Indianapolis, Ind.'
'Kansas City, Mo.'
'Los Angeles, Calif.'
'Louisville, Ky.'
'Memphis, Tenn.'
'Miami, Fla.'
'Minneapolis, Minn.'
'New Orleans, La.'
'New York, N. Y.'
'Omaha, Neb.'
'Philadelphia, Pa.'
'Phoenix, Ariz.'
'Pittsburgh, Pa.'
'St. Louis, Mo.'
'Salt Lake City, Utah'
'San Francisco, Calif.'
'Seattle, Wash.'
'Washington, D.C.'};

names = names';

D = [0	1052	776	578	618	581	1095	641	1152	567	433	579	1802	331	217	665	862	312	864	732	783	1456	608	400	1466	2013	2082	661
1052	0	400	851	551	1551	1769	613	2072	1605	807	1251	2596	826	1137	1255	1123	1359	188	1282	271	2300	483	1038	2099	2699	2493	393
776	400	0	454	173	1198	1370	216	1692	1286	435	861	2198	483	803	1181	731	1086	292	883	279	1906	178	662	1699	2300	2117	292
578	851	454	0	308	803	920	238	1252	940	165	414	1745	269	482	1188	355	833	713	432	666	1453	410	262	1260	1858	1737	597
618	551	173	308	0	1025	1227	90	1525	1114	263	700	2049	311	630	1087	630	924	405	739	360	1749	115	492	1568	2166	2026	306
581	1551	1198	803	1025	0	663	999	572	225	763	451	1240	726	420	1111	862	443	1374	586	1299	887	1070	547	999	1483	1681	1185
1095	1769	1370	920	1227	663	0	1156	557	879	1000	558	831	1038	879	1726	700	1082	1631	488	1579	586	1320	796	371	949	1021	1494
641	613	216	238	90	999	1156	0	1479	1105	240	645	1983	316	623	1152	543	939	482	669	443	1690	205	455	1492	2091	1938	396
1152	2072	1692	1252	1525	572	557	1479	0	676	1264	839	701	1254	976	1643	1157	983	1905	878	1836	346	1590	1034	689	995	1376	1728
567	1605	1286	940	1114	225	879	1105	676	0	865	644	1374	803	484	968	1056	318	1420	794	1341	1017	1137	679	1200	1645	1891	1220
433	807	435	165	263	763	1000	240	1264	865	0	453	1809	107	384	1024	511	712	646	525	585	1499	330	231	1356	1949	1872	494
579	1251	861	414	700	451	558	645	839	644	453	0	1356	480	369	1241	413	680	1097	166	1038	1049	781	238	925	1506	1506	945
1802	2596	2198	1745	2049	1240	831	1983	701	1374	1809	1356	0	1829	1603	2339	1524	1673	2451	1315	2394	357	2136	1589	579	347	959	2300
331	826	483	269	311	726	1038	316	1254	803	107	480	1829	0	320	919	605	623	652	580	582	1508	344	242	1402	1986	1943	476
217	1137	803	482	630	420	879	623	976	484	384	369	1603	320	0	872	699	358	957	529	881	1263	660	240	1250	1802	1867	765
665	1255	1181	1188	1087	1111	1726	1152	1643	968	1024	1241	2339	919	872	0	1511	669	1092	1397	1019	1982	1010	1061	2089	2594	2734	923
862	1123	731	355	630	862	700	543	1157	1056	511	413	1524	605	699	1511	0	1051	1018	290	985	1280	743	466	987	1584	1395	934
312	1359	1086	833	924	443	1082	939	983	318	712	680	1673	623	358	669	1051	0	1171	847	1089	1316	919	598	1434	1926	2101	966
864	188	292	713	405	1374	1631	482	1905	1420	646	1097	2451	652	957	1092	1018	1171	0	1144	83	2145	317	875	1972	2571	2408	205
732	1282	883	432	739	586	488	669	878	794	525	166	1315	580	529	1397	290	847	1144	0	1094	1036	836	354	833	1429	1369	1014
783	271	279	666	360	1299	1579	443	1836	1341	585	1038	2394	582	881	1019	985	1089	83	1094	0	2083	259	811	1925	2523	2380	123
1456	2300	1906	1453	1749	887	586	1690	346	1017	1499	1049	357	1508	1263	1982	1280	1316	2145	1036	2083	0	1828	1272	504	653	1114	1983
608	483	178	410	115	1070	1320	205	1590	1137	330	781	2136	344	660	1010	743	919	317	836	259	1828	0	559	1668	2264	2138	192
400	1038	662	262	492	547	796	455	1034	679	231	238	1589	242	240	1061	466	598	875	354	811	1272	559	0	1162	1744	1724	712
1466	2099	1699	1260	1568	999	371	1492	689	1200	1356	925	579	1402	1250	2089	987	1434	1972	833	1925	504	1668	1162	0	600	701	1848
2013	2699	2300	1858	2166	1483	949	2091	995	1645	1949	1506	347	1986	1802	2594	1584	1926	2571	1429	2523	653	2264	1744	600	0	678	2442
2082	2493	2117	1737	2026	1681	1021	1938	1376	1891	1872	1506	959	1943	1867	2734	1395	2101	2408	1369	2380	1114	2138	1724	701	678	0	2329
661	393	292	597	306	1185	1494	396	1728	1220	494	945	2300	476	765	923	934	966	205	1014	123	1983	192	712	1848	2442	2329	0];

% Create "individual subjects"
D = D ./ 1000;
%X = rand(28,2); D = squareform(pdist(X));

X = cmdscale(D);
X = X(:,1:r);

weightsactual=[1:1:5 0 1:1:5]';
weightsactual(:,2)=[1:5 0 5:-1:1]';

weightsactual = rand(11,2);

for k=1:size(weightsactual,1)			% set of r x r diag matrices with weights on diagonal
		weightsactualk(:,:,k) = diag(weightsactual(k,:));
        
        %dat(:,:,k)=squareform(pdist(X*weightsactualk(:,:,k)));       %
        %this is wrong b.c doesn't weight distance, but rather stim loc.
        %(which is not a ratio scale)
        
        dat(:,:,k) = weighted_pdist(X,weightsactual(k,:));          % w applied to squared dist along each dim 
        
        sp(:,:,k)=X * weightsactualk(:,:,k) * X';  % this is a scalar product matrx, must conv to orig
        
       
end

return


% INPUT SIMULATED SP MATRIX TO TEST
for k=1:size(weightsactual,1)			% set of r x r diag matrices with weights on diagonal
        sp(:,:,k)=X * weightsactualk(:,:,k) * X';  % this is a scalar product matrx, must conv to orig
       %DD = -2 * ( inv(P'*P)*P' * B * P' * inv(P * P') );  % this doesn't work.
end

% GET SIM X FROM SP

[V E] = eig(mean(sp,3)); % guard against spurious complex e-vals from roundoff
[e i] = sort(diag(E)); e = flipud(e); i = flipud(i); % sort descending
keep = find(e > eps^(3/4)*max(abs(e))); % keep only positive e-vals (beyond roundoff)
if isempty(keep)
    X = zeros(n,1);
else
    X = V(:,i(keep)) * diag(sqrt(e(keep)));
end

B2 = []; m = size(weightsactual,1);
for k = 1:m
    Bk = sp(:,:,k);
    B1(k,:) = Bk(:)';
    B2 = [B2; Bk];
end

% INDSCAL ON X (compare fit using sp)

[fit,X2,W,XL,XR] = indscalf(X,sp,B1,B2);
fn = nmdsfig(X,ones(size(X,1),1),names); % Make a figure to check
hold on; plot(X2(:,1), X2(:,2),'bo','MarkerFaceColor','b')




