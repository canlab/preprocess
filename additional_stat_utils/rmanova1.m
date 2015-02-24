function stats = rmanova1(data,alpha,bxplot,ttst);

%Repeated-measures one-way ANOVA
%
%USAGE: stats = rmanova1(data,[alpha],[bxplot],[ttst]);
%  inputs: 
%    data = samples x treatments matrix
%    alpha (optional): p-value threshold (default: 0.05)
%    bxplot (optional): if 1, will produce a boxplot (default: 0)
%    ttst (optional): if 1, will perform pairwise t-tests between all (default: 0)

if nargin < 2; alpha = 0.05; bxplot = 0; ttst = 0; end
if nargin < 3; bxplot = 0;  ttst = 0; end
if nargin < 4; ttst = 0; end

% remove nans
[wasnan, data] = nanremove(data);

k = size(data,2);           %number of treatment levels
n = size(data,1);           %sample size

%compute stats for each group
for g = 1:k;
    avg(g) = mean(data(:,g));
    SS(g) = sum((data(:,g)-avg(g)).^2);
    T(g) = sum(data(:,g));
end

AVGT = sum(sum(data.^2));
G = sum(sum(data));
NT = size(data,1) * size(data,2);
P = sum(data,2);

%calculate the F ratio
SST = AVGT - ((G^2)/NT);
SSBT = sum(((T.^2)./n)) -((G^2)/NT);
SSWT = sum(SS);
SSBS = sum(((P.^2)./k)) -((G^2)/NT);
SSE = SSWT - SSBS;
dft = NT - 1;
dfbt = k - 1;
dfwt = NT - k;
dfbs = n - 1;
DFE = dfwt - dfbs;
MSBT = SSBT / dfbt;
MSE = SSE / DFE;
F = MSBT / MSE;


stats.SS.between_treatments = SSBT;
stats.SS.within_treatments = SSWT;
stats.SS.between_subjects = SSBS;
stats.SS.error = SSE;
stats.SS.total = SST;

stats.df.between_treatments = dfbt;
stats.df.within_treatments = dfwt;
stats.df.between_subjects = dfbs;
stats.df.error = DFE;
stats.df.total = dft;

stats.MS.between_treatments = MSBT;
stats.MS.error = MSE;

stats.F = F;
stats.P = 1- fcdf(F,dfbt,DFE);
stats.alpha = alpha;

%decision rule
if stats.P < alpha;
    stats.significant = 'Yes';
else
    stats.significant = 'No';
end

% print output line

fprintf('Omnibus test: F(%3.0f, %3.0f) = %3.2f, MSE = %3.2f, p = %3.6f\n', stats.df.between_treatments, stats.df.error, stats.F, stats.MS.error, stats.P)

%make boxplot
if bxplot;
    boxplot(data);
end

%pairwise t-tests
if ttst;
    %disp('Performing post-hoc t-tests...');
    tcounter = 0;
    for t = 1:k;
        for t2 = 1:k;
            if t2 > t;
                tcounter = tcounter + 1;
                [h,p,ci,stat] = ttest(data(:,t),data(:,t2));
                stats.ttests(tcounter).P = p;
                stats.ttests(tcounter).comparison = [num2str(t),' > ',num2str(t2)];
                stats.ttests(tcounter).means = [mean(data(:,t)) mean(data(:,t2))];
                stats.ttests(tcounter).ci = ci';
            end
        end
    end
end