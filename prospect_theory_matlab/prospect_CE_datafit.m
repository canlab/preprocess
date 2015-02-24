%%
% assumes data is a mat file with matrix called data, from imported excel
% file
%
% excel file is exported from e-prime, and first column and first row are
% removed and saved as .xls, then imported into matlab and save as .mat
% file
%
% batch mode only; if doing single subject, read instructions below in this
% cell and start from cell 3
%
d = dir('*ce.mat');

% if you don't want to use batch mode, comment out the outside 'for' loop
% in the 3rd cell (for i=1:length(d))) and it's 'end' line; comment out the
% 3 lines following the for statement (clear... fname = ...); load the
% desired .mat file in the workspace and assign a filename (e.g. filename =
% '1009');
%
% if you're doing one subject at a time, also comment out lines that assign values to variables starting with
% "group" such as those in cell 2.
%%

group_prices = nan(50,length(d));   % certainty equivalences for each subject
group_fitu = nan(50,length(d));     % fitted utilities
group_a = nan(length(d),1);         % scaling parameter of weighting function
group_b = nan(length(d),1);         % exponent of weighting function
group_m = nan(length(d),1);         % exponent of value function
group_lo = nan(length(d),1);        % loss aversion parameter 
group_ev = nan(50,length(d));       % expected value of each gamble 


%%

% for individual plots
plotp = [.001:.001:.999];

% script designed for batch mode to print out a formatted table of
% variables where each row is a subject
nms = {'subid' 'a' 'b' 'm' 'lo' 'Stupid Answers' 'Fit' 'Risky picks total' 'Risky high prob' 'Risky low prob' 'Risky Gain' 'Risky Loss' 'Risky Mixed' 'Gain Hi' 'Gain Lo' 'Loss Hi' 'Loss Lo' 'Mix Hi' 'Mix Lo' 'TestGamble R' 'TestGamble P' 'Mean Pricedevs' 'SD PriceDevs'};
fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');

for i = 8:length(d),
    clear data best_params fit_u value1 y p cert_equiv name fname;
    load(d(i).name);
    fname = d(i).name(1:4);
    
% gets data (column numbers are based on imported data from excel)
value1 = data(:,66);    % becomes 'x' in equations; highest absolute value of two choices in loss or gain trials
y = data(:,67);         % alternate choice in gamble
p = data(:,46);         % probability of getting choice x (value 1)

% searches for trials on which no choice was made and removes them
wh= find(data(:,21)==0);    
data(wh,:) = [];
value1(wh) = [];
y(wh) = [];
p(wh) = [];

trialtype = NaN(length(value1),1); 

% in experiment, values listed under x refer to value presented 1st, but
% not necessarily highest absolute value, so below, we re-order the columns
% to make sure that the highest absolute value is in the x column
for j=1:length(value1),
    tmp = [value1(j) y(j)]; abstmp = [abs(value1(j)) abs(y(j))];
    if tmp(1).*tmp(2) <0,       % if negative product, it must be a mixed gamble
        whmax = find(tmp == max(tmp));  % for mixed gambles, order doens't matter really
        trialtype(j) = 0;
    else
        %for gain trial, x > y, and for loss trial x < y, so making sure
        %the highest absolute value is in the x position
        whmax = find(abstmp == max(abstmp)); 
    end
    
      if length(whmax) > 1, whmax = whmax(1); end
      if whmax ==1, whmin = 2; elseif whmax ==2, whmin =1; end
    newx(j,1) = tmp(whmax);
    newy(j,1) = tmp(whmin);
    
       
    if value1(j)==newx(j)   
        newp(j,1) = p(j); %if positions didn't change, keep p the same
    else newp(j,1) = 1-p(j); %if positions switched, then change p to 1-p
    end
    
      if value1(j) + y(j) > 0,
        trialtype(j) = 1;
    elseif value1(j) + y(j) < 0,
        trialtype(j) = -1;
    end
end

% reassign values 
value1 = newx;
y = newy;
p = newp;

%index for trialtype
wh_gain = find(trialtype == 1);
wh_loss = find(trialtype == -1);
wh_mix = find(trialtype == 0);

ev = (value1.*p)+(y.*(1-p));  % expected value of each gamble
%group_ev(1:length(ev),i) = ev; 

cert_equiv_index = zeros(size(data,1),1);
gamble_price = zeros(size(data,1),1);
adjusted_price = zeros(size(data,1),1);
cert_equiv = zeros(size(data,1),1);
% these are columns specific to the way data is imported from excel; 
% if you don't export the exact same way, these numbers will be incorrect
index = [35,38:1:45,36,37]; 

for j=1:size(data,1),
    cert_equiv_index(j) = data(j,21); %index is number 1-11, indicating which box they chose
    gamble_price(j) = data(j,index(cert_equiv_index(j))); %finds the value in that box for that trial
    
    if gamble_price(j) > max([value1(j,1) y(j,1)]),
        adjusted_price(j) = max([value1(j,1) y(j,1)]);
        cert_equiv(j) = adjusted_price(j);
    elseif gamble_price(j) < min([value1(j,1) y(j,1)]),
        adjusted_price(j) = min([value1(j,1) y(j,1)]);
        cert_equiv(j) = adjusted_price(j);
    else adjusted_price(j) = gamble_price(j);
        if index(cert_equiv_index(j)) < 11,
            cert_equiv(j) = (adjusted_price(j) + data(j,index(cert_equiv_index(j)+1)))./2; 
        else cert_equiv(j) = adjusted_price(j);
        end
    end
end


% testing for montonicity in pricing same gambles with diff probability
wh = find(abs(value1)==800 & abs(y) == 800);
tmpx = value1(wh); tmpy = y(wh); tmpp = p(wh);tmpce = cert_equiv(wh);
for j = 1:length(tmpx),
    if tmpx(j) < 0,
        p(j) = 1-p(j);
    end
end
[ordprob,ordind] = sort(p(wh));
%scatter(ordprob,tmpce(ordind));
[r,pval] = corr(ordprob,tmpce(ordind));
%plot it
scatter(ordprob,tmpce(ordind));
xlabel('Probability of X');
ylabel('Pricing of Gamble');
name = [fname '_trier_ce_plots'];
saveas(gcf,name, 'fig'); close(gcf);

%pricing deviation from expected values
pricedev_mean = nanmean(cert_equiv - ev);
pricedev_sd = nanstd(cert_equiv - ev);

stupid_answers = gamble_price ~= adjusted_price;

prob_hi = p > .5;
prob_lo = p <= .5;

risk_total = nanmean(cert_equiv > ev); % number of choices on which the subject subjectively valued gamble more than the objective EV
risk_hi = nanmean(cert_equiv(prob_hi)>ev(prob_hi));
risk_low = nanmean(cert_equiv(prob_lo)>ev(prob_lo));

risk_gain = nanmean(cert_equiv(wh_gain)>ev(wh_gain));
  risk_gain_hi = nanmean(cert_equiv(wh_gain(prob_hi(wh_gain)==1))>ev(wh_gain(prob_hi(wh_gain)==1)));
  risk_gain_lo = nanmean(cert_equiv(wh_gain(prob_lo(wh_gain)==1))>ev(wh_gain(prob_lo(wh_gain)==1)));
risk_loss = nanmean(cert_equiv(wh_loss)>ev(wh_loss));
  risk_loss_hi = nanmean(cert_equiv(wh_loss(prob_hi(wh_loss)==1))>ev(wh_loss(prob_hi(wh_loss)==1)));
  risk_loss_lo = nanmean(cert_equiv(wh_loss(prob_lo(wh_loss)==1))>ev(wh_loss(prob_lo(wh_loss)==1)));
risk_mix = nanmean(cert_equiv(wh_mix)>ev(wh_mix));
  risk_mix_hi = nanmean(cert_equiv(wh_mix(prob_hi(wh_mix)==1))>ev(wh_mix(prob_hi(wh_mix)==1)));
  risk_mix_lo = nanmean(cert_equiv(wh_mix(prob_lo(wh_mix)==1))>ev(wh_mix(prob_lo(wh_mix)==1)));


  
%group_prices(1:length(cert_equiv),i) = cert_equiv;

warning off
%runs fitting function
%make last input argument "1" if you want to include plots from
%prospect_fit_data.m
[best_params,fit_u,fval] = prospect_fit_data(value1,y,p,cert_equiv,0);

% group_a(i) = best_params(1);
% group_b(i) = best_params(2);
% group_m(i) = best_params(3);
% group_lo(i) = best_params(4);
% group_fitu(1:length(fit_u),i) = fit_u;

warning on

    % a = best_params(1);   b = best_params(2);
    % a : scaling. increase = convex, < 0 = hyperbolic, 1 = linear
    % a > 0 & a < 1, averse to gambling (probs are lower than nominal)
    % a > 1, seeks gambling; probs are higher than nominal
    %
    % b : exponent.  increase: sigmoid shape, < 1, reverse slope
    % b > 0 & b < 1, overweight low prob, underweight high
    % b > 1, underweight low prob, overweight high
    
        % a = scale, weighting function
    % b = exp., weighting function
    % m = exp., value function
    % lo = loss aversion

% % individual plots
%     xvals = [-100:100]; len = length(xvals); xvals = xvals';
%     u = prospect_utility(best_params,xvals,zeros(1,len),.5 * ones(1,len));
% 
%     
%     subplot(1,3,1);
%     
%     plot(xvals,u,'LineWidth',2); hold on;
%     plot_vertical_line(0);
%     title('Utility Function (at p = .5)')
%     xlabel('Nominal value');
%     ylabel('Est. Subjective Utility')
% 
% fp = best_params(1) .* (plotp.^best_params(2));
% wp = fp ./ (fp + (1-plotp).^best_params(2));
% 
% subplot(1,3,2);
% plot(plotp,wp,'LineWidth',2)
% title('Probability Weighting Function')
% xlabel('Probability');
% ylabel('Subjective Probability')
% 
% % scatterplot of pricing vs. gamble
% subplot(1,3,3);
% scatter(ev,cert_equiv);
% hold on; plot([min(ylim) max(ylim) ],[min(ylim) max(ylim) ],'k-');
% xlabel('Expected Value of Gamble');
% ylabel('Pricing of Gamble');
% 
% name = [fname '_trier_ce_plots'];
% saveas(gcf,name, 'fig'); %close(gcf);

%name = [fname '_ce_analysis2'];
%save(name)

%nms = {'a' 'b' 'm' 'lo' 'Stupid Ans'};
%fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
fprintf(1,'%s\t',d(i).name(1:4));
fprintf(1,'%3.4f\t',best_params);
fprintf(1,'%1.0f\t',nansum(stupid_answers));
fprintf(1,'%1.0f\t',fval);
fprintf(1,'%3.2f\t',risk_total); 
fprintf(1,'%3.2f\t',risk_hi); 
fprintf(1,'%3.2f\t',risk_low); 
fprintf(1,'%3.2f\t',risk_gain); 
fprintf(1,'%3.2f\t',risk_loss); 
fprintf(1,'%3.2f\t',risk_mix); 
fprintf(1,'%3.2f\t',risk_gain_hi); 
fprintf(1,'%3.2f\t',risk_gain_lo); 
fprintf(1,'%3.2f\t',risk_loss_hi); 
fprintf(1,'%3.2f\t',risk_loss_lo); 
fprintf(1,'%3.2f\t',risk_mix_hi); 
fprintf(1,'%3.2f\t',risk_mix_lo); 
fprintf(1,'%3.2f\t',r); 
fprintf(1,'%3.4f\t',pval); 
fprintf(1,'%3.4f\t',pricedev_mean);
fprintf(1,'%3.4f\t',pricedev_sd); 
fprintf(1,'\n'); 

end


%%
%
        