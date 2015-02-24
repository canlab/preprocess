function [xout,IDX] = mixture_model(tmp,varargin)
% [xout,IDX] = mixture_model(tmp,[verbose])
%
% tor wager
% this may be a silly way, not a real mixture model, but this function does
% a test for normality on data (tmp), and if it fails, uses k-means to
% cluster values into K classes.  if data within each class is normally
% distributed (via significance test at fixed p), then it stops, otherwise it proceeds
% up to max K = 5.  The function then subtracts class means.  This is
% designed for timeseries data in which a true signal is masked by overall
% 'cycling' of mean fMRI signal among two or more baseline values.  
% Needs the stats toolbox.
%
% normality test with arbitrary mean and variance is done using jbtest.m
                        % jbtest is very sensitive to outliers, so we
                        % windsorize first - but if the number trimmed is greater
                        % than 10% of the observations, we move
                        % automatically to clustering.
%                        
% Convergence criteria (any of these signals convergence - "OR"):
%   - all classes pass normal distribution test
%   - less than 10% of points in a given class  (if so, returns to previous
%       K value that meets this requirement)
%   - K = 5
%
% returns class ID, silhouette index for clusters (0 to 1, clustering
% goodness), and adjusted data.

verbose = 1; if length(varargin) > 0, verbose = varargin{1};,end

minp = length(tmp) .* .1;   % at least this many observations in each class
pthr = .01;                  % p-value for significance

[tmp2,ntrimmed] = trimts(tmp,3,[]);
[h,p] = jbtest(tmp2); % normality test with arbitrary mean and variance
                        % jbtest is very sensitive to outliers, so we trim
                        % first - but if the number trimmed is greater
                        % than 10% of the observations, we move
                        % automatically to clustering.
tmp = tmp2;           % use this, so kmeans is not influenced by outliers as well.
h = p < pthr;

if verbose, fprintf(1,'Initial test of normality: H = %3.0f, p = %3.4f, windsor = %3.0f\n',h,p,ntrimmed);, end

% gamma function test
%[tmp2,x] = hist(tmp,min(30,length(tmp))); tmp2 = tmp2 ./ max(tmp2);
%[phat] = gamfit(tmp2)
%gcdf = gamcdf(tmp2,phat(1),phat(2))
%h = kstest(tmp2',[tmp2' gcdf'])

           if h | ntrimmed > minp   % not normal - how many modes?
               
               stop = 0; K = 2; avgsil = [];
               while ~stop
                   try
                       [IDX, C,SUMD,D] = kmeans(tmp, K);
                   catch
                       try
                          [IDX, C,SUMD,D] = kmeans(tmp, K);
                      catch
                          IDX = ones(size(tmp));
                          stop = 1;
                      end
                  end
                  
                   clear h p n
                   for c = 1:max(IDX)
                       [tmp2,ntrimmed] = trimts(tmp(IDX==c),3,[]);
                       [h(c),p(c)] = jbtest(tmp2);        % normality test
                       n(c) = sum(IDX==c);                % number per class
                   end
                   
                   % re-threshold
                   h = p < pthr;
                   
                   if verbose
                       fprintf(1,'K = %3.0f, Hypothesis test for normality:\n',K)
                       disp(['H = ' num2str(h)]),disp(['p = ' num2str(p)])
                       disp(['Class n = ' num2str(n)]);
                   end
                   
                   kd = D(IDX==K,:);  % distances to all centers for class K
                  tmp2 = kd; tmp2(:,K)=[]; wh = find(sum(tmp2) == min(sum(tmp2))); % next closest cluster
                  dif = (tmp2(:,wh) - kd(:,K)) ./ max([tmp2(:,wh) kd(:,K)]')';   % silhouette indices
                  avgsil(K) = mean(dif);
                       
                   if ~any(h) | K == 5 | length(unique(IDX)) < K | any(n < minp), 
                       stop = 1;, 

                       if (length(unique(IDX)) < K | any(n < minp)) & K > 2
                           % We've gone too far - empty class or class <
                           % minp - go back 1
                           K = K - 1;
                           [IDX, C,SUMD,D] = kmeans(tmp, K);
                   
                           clear h p n
                            for c = 1:max(IDX)
                                [tmp2,ntrimmed] = trimts(tmp(IDX==c),3,[]);
                                [h(c),p(c)] = jbtest(tmp2);        % normality test
                                n(c) = sum(IDX==c);                % number per class
                            end
                   
                           % re-threshold
                            h = p < pthr;
                   
                            if verbose
                                fprintf(1,'Too few members in class!  Going back.\n')
                                fprintf(1,'K = %3.0f, Hypothesis test for normality:\n',K)
                                disp(['H = ' num2str(h)]),disp(['p = ' num2str(p)])
                                disp(['Class n = ' num2str(n)]);
                            end
                   
                            kd = D(IDX==K,:);  % distances to all centers for class K
                            tmp2 = kd; tmp2(:,K)=[]; wh = find(sum(tmp2) == min(sum(tmp2))); % next closest cluster
                            dif = (tmp2(:,wh) - kd(:,K)) ./ max([tmp2(:,wh) kd(:,K)]')';   % silhouette indices
                            avgsil(K) = mean(dif);
                        end
                           
                   else, K = K + 1;, 
                   end
   
             end
             
         else
             if verbose, disp('Data appears to be normally distributed')
             end
         end
            
         
         % adjust
         xout = tmp;
         if exist('IDX') == 1
            for c = 1:max(IDX), 
                xout(IDX==c) = tmp(IDX==c) - mean(tmp(IDX==c));
            end
                 
         else
            IDX = ones(size(tmp));
            xout = tmp;
         end
         
         
         % output plot
         
             if verbose
                 
                 figure; 
                 subplot(1,3,1);plot(tmp); hold on;plot(xout,'r'); legend({'Original' 'Declassed'})
                 subplot(1,3,2);hist(tmp,min(30,length(tmp)));
                 [hi,x] = hist(tmp,min(30,length(tmp)));
                 hold on; cols = {'r' 'g' 'y' 'c' 'm'};
                 if max(IDX) > 1
                     for c = 1:max(IDX), hh = hist(tmp(IDX==c),x); bar(x,hh,cols{c});,end
                 
                    subplot(1,3,3); plot(2:length(avgsil)+1,avgsil,'LineWidth',3),xlabel('K');,ylabel('Silhouette avg.')
                end
             end
             
             
             return
             
             