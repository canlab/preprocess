function LMneg = gambling_negLM(phats,data,actual_choices,good,bad,fmode,st,en)
        %LMneg = gambling_negLM(phats,data,actual_choices,good,bad,fmode,st,en)
switch fmode
    case 'expval'
        phats(1) = max([0 phats(1)]);   % phat(1) = w, and we want to constrain it to be between 0 and 1
        phats(1) = min([1 phats(1)]);
        phats(2) = max([0 phats(2)]);
        phats(2) = min([1 phats(2)]);
        

        pchoice = predict_choices(phats,data,st,en);
        
    case 'decayEV'
        [pchoice,EVd] = decayed_expectedvalue(phats,data,st,en);
        
    case 'heuristic'
        phats(3) = max([0 phats(3)]);   % phat(1) = w, and we want to constrain it to be between 0 and 1
        phats(3) = min([1 phats(3)]);
        [pchoice,S,c] = gambling_heuristic_choice(phats,data,good,bad);
        
    case 'baseline'
        phats(1) = max([0 phats(1)]); phats(2) = max([0 phats(2)]); phats(3) = max([0 phats(3)]); phats(4) = max([0 phats(4)]);
        phats = phats./sum(phats);
        pchoice = [phats(1), phats(2), phats(3), phats(4)];
        pchoice = repmat(pchoice,100,1);
        
    case 'softmax'
        pchoice = gambling_softmax(phats,data,actual_choices);
      
end
%data = pchoice(actual_choices >0);
%what's the model's prediction of choosing the chosen deck?
pchosen = sum(pchoice .* actual_choices(st:en,:),2); 
%changes 0 to small number --> log of 0 is -INF
wh = find(pchosen<=0);
pchosen(wh) = .0001;
% sum of log likelihoods
LM = nansum(log(pchosen));
LMneg = -1*LM;
%LMneg = sum((pchosen(20:end)-data{4}(20:end)).^2);

return
