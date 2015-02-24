%allchosen = zeros(100,74);
%allchosensmooth = zeros(100,74);
%neglikelihood = zeros(74,1);
for i = 1:length(d)
    fprintf(1,'%3.0f ',i);
    clear data, load(d(i).name);
    [pchoice,LMneg,pchosen,actual_choices] = gambling_winstay(data);
    winstay_pchoice(:,i) = pchosen;
    winstay_preds(i) = mean(winstay_pchoice(51:end,i));
    %allchosensmooth(:,i) = moving_average('gaussian',pchosen,20);
    %neglikelihood(i,1) = LMneg; 
end

% figure; m = mean(allchosensmooth'); sterr = ste(allchosensmooth'); plot(m,'LineWidth',3);
% fill_around_line(m,sterr,'b');