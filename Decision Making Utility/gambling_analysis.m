function gambling_analysis(data,fmode,st,en);
%[pchoice,w,a,c,LMneg,reorder_decks] = gambling_analysis(data,fmode,st,en);
%[pchoice,actual_choices,goodpicks,smgoodpicks,w,a,c,LMneg,smooth_choices,smooth_pchoices] = gambling_analysis(data,fmode,st,en);
% choices = vector representing the observed choices by subject
% pchoice = predicted choice probabilities, calculated from
% predict_choices.m


tmp = (st:1:en); ntrials = length(tmp);  %allows for variable start and end positions to analyze separate portions (e.g. reversal 101-180)
actual_choices = zeros(ntrials,4);      
tmpchoice = data{2}(st:en,1);
wh1 = find(data{2}(st:en,1)==1);
wh2 = find(data{2}(st:en,1)==2);
wh3 = find(data{2}(st:en,1)==3);
wh4 = find(data{2}(st:en,1)==4);

actual_choices(wh1,1) = 1;
actual_choices(wh2,2) = 1;
actual_choices(wh3,3) = 1;
actual_choices(wh4,4) = 1;

% good = strmatch('Advan',data{1});   %nominally good decks
% bad = strmatch('Disadvan',data{1}); %nominally bad decks

%put all of the decks in the same order for every subject; these labels
%"reorder_decks" indicates the index by which to sort other data if
%collapsing across subjects
% reorder_decks(1) = strmatch('Advantageous small',data{1}); 
% reorder_decks(2) = strmatch('Advantageous large',data{1});
% reorder_decks(3) = strmatch('Disadvantageous small',data{1});
% reorder_decks(4) = strmatch('Disadvantageous large',data{1});

switch fmode
    case 'expval'   % for expectancy valence model (Busemeyer & Stout)
        %startvalues = [.8, .2, 1.5];
        start_w = [0 .15 .25 .5 .75 1]; % w = weighting factor for wins v. losses
        start_a = [.01 .1 .5 .75 1 1.5];    % a = learning rate
        start_c = [0 .25 .5 .75 1 2];   % c = sensitivity parameter

        w = zeros(5,5,5);
        a = zeros(5,5,5);
        c = zeros(5,5,5);

        t = 100; % number of trials

        % for i = 1:length(t)
        %    tmp_lm(i) = actual_choices(i,1).*ln(pchoice(i,1)) + actual_choices(i,2).* ln(pchoice(i,2)) + ...
        %    actual_choices(i,3).*ln(pchoice(i,3)) + actual_choices(i,4).* ln(pchoice(i,4));
        % end
        % LM = sum(tmp_lm);
        %
        %LM = sum(ln(pchoice(actual_choices > 0))

        %minimize(-LM), - log likelood of model

        %f= @(phats) sum(nansum((predict_choices(phats,data)-actual_choices).^2));
        
       

     
 
        
%         f = @(phats) sum(nansum((predict_choices(phats,data)-probchoice).^2));
%         warning off
%         for i = 1:length(start_w),
%             for j = 1:length(start_a),
%                 for k = 1:length(start_c),
%                     [phats,fval(i,j,k)] = fminsearch(f,[start_w(i);start_a(j);start_c(k)]);
%                     w(i,j,k) = phats(1); a(i,j,k) = phats(2); c(i,j,k) = phats(3);
%                 end
%             end
%         end
%         warning on
warning off
        f = @(phats) gambling_negLM(phats,data,actual_choices,[],[],'expval',st,en) % minimizes negative log likelihood of expectancy valence model
        
        for i = 1:length(start_w),
            for j = 1:length(start_a),
                for k = 1:length(start_c),
                    %[phats,fval(i,j,k),exitflag] = fminsearchbnd3(f,[start_w(i);start_a(j);start_c(k)],[0 0 0], [1 1 Inf]);
                    [phats,fval(i,j,k),exitflag] = fminsearchbnd3(f,[start_w(i);start_a(j);start_c(k)]);
                    w(i,j,k) = phats(1); a(i,j,k) = phats(2); c(i,j,k) = phats(3);
                end
            end
        end
    
        %[phats,fval] = fminsearch(f,startvalues);
%         wh = find(w>0 & w<1 & a>0 & c>0);
%         w = w(wh); a = a(wh); c = c(wh); fval = fval(wh);
       
     
        [x,y,z] = matrix_min(fval);
        w = w(x,y,z);
        a = a(x,y,z);
        c = c(x,y,z);
        %SSE = fval(x,y,z);
        
        LMneg = fval(x,y,z);
        %varargout = {v,ev,strength,theta};
       [pchoice,v,ev,strength,theta] = predict_choices([w,a,c],data,st,en); %re-runs the expectancy valence model with chosen parameters
%         LM = nansum(log(pchoice(actual_choices > 0)));

        pchosen = sum(pchoice .* actual_choices,2);
        modelprob = mean(pchosen(end-10:end));
        
        % determine which choices were from good/bad decks --> calculate
        % the proportion of picks from each deck
         choices = data{2}(st:en,1);
         goodpicks = zeros(100,2);
         for i = 1:ntrials,
             goodpicks(i,1) = ismember(choices(i),good); 
             goodpicks(i,2) = ismember(choices(i),bad);
         end
        
        % moving average for illustration of changing probability of
        % choosing from good decks
%         smgoodpicks = moving_average('gaussian',goodpicks(:,1),20);    
%        
%         % moving average of choices from each deck
%         smooth_choices = moving_average('gaussian',actual_choices,20);
%         % moving average of probability of choice from each deck according
%         % to the expectancy valence model
%         smooth_pchoices = moving_average('gaussian',pchoice,20);
        
%         ntrials = 20;
%         kern = normpdf(-3:6/ntrials:3); kern = kern./sum(kern);
%         tmp1 = conv(probchoice(:,1),kern); tmp1 = tmp1(1:100);
%         tmp2 = conv(probchoice(:,2),kern); tmp2 = tmp2(1:100);
%         tmp3 = conv(probchoice(:,3),kern); tmp3 = tmp3(1:100);
%         tmp4 = conv(probchoice(:,4),kern); tmp4 = tmp4(1:100);
%         deckprops = [tmp1,tmp2,tmp3,tmp4];
% 
%         tmp1 = conv(pchoice(:,1),kern); tmp1 = tmp1(1:100);
%         tmp2 = conv(pchoice(:,2),kern); tmp2 = tmp2(1:100);
%         tmp3 = conv(pchoice(:,3),kern); tmp3 = tmp3(1:100);
%         tmp4 = conv(pchoice(:,4),kern); tmp4 = tmp4(1:100);
% 
%         smpchoice = [tmp1,tmp2,tmp3,tmp4];
% 
%         t=1:100;
%         figure; plot(t,y1(:,1),'b-',t,y2(:,1),'b-.',t,y1(:,2),'g-',t,y2(:,2),'g-.',...
%             t,y1(:,3),'r-',t,y2(:,3),'r-.',t,y1(:,4),'c-',t,y2(:,4),'c-.');
%         %title(files(1).name(1:3));
%         legend(data{1}{1}, 'predicted', data{1}{2}, 'predicted', data{1}{3},'predicted',data{1}{4},'predicted','Location','NW');
        
        %figure; plot(t,smgoodpicks);
% 
         avert = mean(data{2}(st:en,2));
         averwd = mean(data{2}(st:en,3));
% 
         pick1 = length(find(data{2}(st:en,1)==1));
         pick2 = length(find(data{2}(st:en,1)==2));
         pick3 = length(find(data{2}(st:en,1)==3));
         pick4 = length(find(data{2}(st:en,1)==4));
         
         picks = [pick1,pick2,pick3,pick4];

% group_goodpicks(:,ind) = goodpicks(:,1); group_smgoodpicks(:,ind) = smgoodpicks(:,1);
% group_choices(:,:,ind) = actual_choices(:,:,1);
% %group_smchoices(:,:,ind) = smooth_choices(:,:,1);
% group_pchoice(:,:,ind) = pchoice(:,:,1);
% group_smpchoice(:,:,ind) = smooth_pchoices(:,:,1);


%avert = nanmean(data{2}(:,2));  % average RT to choose overall
%averwd = nanmean(data{2}(:,3)); % average reward over all trials
%pgoodpicks = mean(goodpicks(:,1));
pgoodpicks = mean(goodpicks(:,1));  % proportion of trials picked from "good" decks

        fprintf(1,'\n'); fprintf(1,'\n');
        nms = {'w', 'a', 'c', 'LMneg','totreward','Ave RT', 'Ave Rwd','Prop Good'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',w); fprintf(1,'%3.4e\t',a); fprintf(1,'%3.4f\t',c);
        fprintf(1,'%3.4f\t',LMneg); fprintf(1,'%3.4f\t',sum(data{2}(st:en,3)));
        fprintf(1,'%3.4f\t',avert); fprintf(1,'%3.4f\t',averwd);
        fprintf(1,'%3.4f\t',pgoodpicks);
        fprintf(1,'\n');
        fprintf(1,'\n');
        fprintf(1,'%s\t',[data{1}{1}, data{1}{2}, data{1}{3}, data{1}{4}]); fprintf(1,'\n');
        fprintf(1,'%3.1f\t',picks); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',modelprob);
        fprintf(1,'\n');
        fprintf(1,'\n');
        
        
        %save([files(1).name(1:3) '_expval']);
    
    case 'softmax'
        warning off
        
        start_B = [-5 -1 -.5 0 .5 1 5];
        f = @(phats) gambling_negLM(phats,data,actual_choices,[],[],'softmax',st,en) % minimizes negative log likelihood of expectancy valence model
        
        for i = 1:length(start_B),
                    [phat_explore(i,1),fval(i,1)] = fminsearch(f,start_B(i));            
        end

        whreal = find(~isinf(fval));
        wh = find(fval(whreal)==min(fval(whreal))); wh = wh(1);
        explore_beta = phat_explore(whreal(wh));
        LMneg = fval(whreal(wh));
        pchoice = gambling_softmax(explore_beta,data,actual_choices);
        
        pchosen = sum(pchoice .* actual_choices,2);
        
        fprintf(1,'\n'); fprintf(1,'\n');
        nms = {'LMneg','Beta', 'blank'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',LMneg); fprintf(1,'%3.4f\t',explore_beta);fprintf(1,'%3.4f\t',100.00)
        fprintf(1,'\n');
        
       

    case 'decayEV'  % assumes that expected value decays over time by an exponential weighting function
        actual_choices = zeros(100,4);      
        tmpchoice = data{2}(:,1);
        wh1 = find(data{2}(:,1)==1);
        wh2 = find(data{2}(:,1)==2);
        wh3 = find(data{2}(:,1)==3);
        wh4 = find(data{2}(:,1)==4);

        % decks are entered in original order; varying by subject
        actual_choices(wh1,1) = 1;
        actual_choices(wh2,2) = 1;
        actual_choices(wh3,3) = 1;
        actual_choices(wh4,4) = 1;

        wltmp = data{2}(:,3); %wltmp = wltmp./1250;
        winloss = zeros(length(wltmp),4);
        winloss(wh1,1) = wltmp(wh1);
        winloss(wh2,2) = wltmp(wh2);
        winloss(wh3,3) = wltmp(wh3);
        winloss(wh4,4) = wltmp(wh4);

        EVd = zeros(100,4);     % decayed expected values
        choices = data{2}(:,1);



        warning off
        f = @(phat) gambling_negLM(phat,data,actual_choices,[],[],'decayEV',st,en)
        
        start_a = [0, .05, .5, 1 2];  %decay parameter
        for j = 1:length(start_a),
            [phat(j),fval(j)] = fminsearch(f,start_a(j));
        end
        phat(j)
        whmin = find(fval==min(fval)); whmin = whmin(1);
        a = phat(whmin); fval = fval(whmin);
        
        [pchoice,EVd,l,h,decay] = decayed_expectedvalue(a,data,st,en);    %re-run model with chosen params
        
        %reorder output to be consistent deck order across subjects (saved
        %in same order)
        pchoice = pchoice; EVd = EVd;
        l = l; h = h; decay = decay;
        
        fprintf(1,'\n'); fprintf(1,'\n');
        nms = {'LMneg','a'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',fval); fprintf(1,'%3.4f\t',a);
        fprintf(1,'\n');
        
        save decayEVdata.mat 
  
   
    case 'heuristic'
        
        data{1}{:}
        good = input('good:    ');
        bad = input('bad:    ');
        
        actual_choices = zeros(100,2);      %
        tmpchoice = data{2}(:,1);
        wh1 = find(data{2}(:,1)==good(1) | data{2}(:,1)==good(2));
        wh2 = find(data{2}(:,1)==bad(1) | data{2}(:,1)==bad(2));

        actual_choices(wh1,1) = tmpchoice(wh1);
        actual_choices(wh2,2) = tmpchoice(wh2);
        
        start_a = [.1 .25 .5 1 1.5 2];
        start_b = [.1 .25 .5 1 1.25 1.5];
        start_p1 = [.1 .15 .25 .5 .75 1];

        a = zeros(5,5,5);
        b = zeros(5,5,5);
        p1 = zeros(5,5,5);

        f = @(phats) gambling_negLM(phats,data,actual_choices,good,bad,'heuristic')
        warning off
        for i = 1:length(start_a),
            for j = 1:length(start_b),
                for k = 1:length(start_p1),
                    [phats,fval(i,j,k)] = fminsearch(f,[start_a(i);start_b(j);start_p1(k)]);
                    a(i,j,k) = phats(1); b(i,j,k) = phats(2); p1(i,j,k) = phats(3);
                end
            end
        end
        warning on
        
        [x,y,z] = matrix_min(fval);
        a = a(x,y,z);
        b = b(x,y,z);
        p1 = p1(x,y,z);
        LMneg = fval(x,y,z);
        
        [pchoice,S,c] = gambling_heuristic_choice([a,b,p1],data,good,bad);
        
        deckchoice = zeros(100,2);
        for i = 1:100,
            if data{2}(i,1)==good(1),
                deckchoice(i,1) = 1;
            elseif data{2}(i,1)==good(2),
                deckchoice(i,1) = 1;
            elseif data{2}(i,1)==bad(1),
                deckchoice(i,2) =1;
            elseif data{2}(i,1)==bad(2),
                deckchoice(i,2) = 1;
            end
        end

        probchoice = zeros(100,2);
        for i=1:100,
            probchoice(i,:) = sum(deckchoice(1:i,:))./i;
        end

        ntrials = 20;
        kern = normpdf(-3:6/ntrials:3); kern = kern./sum(kern);
        tmp1 = conv(probchoice(:,1),kern); tmp1 = tmp1(1:100);
        tmp2 = conv(probchoice(:,2),kern); tmp2 = tmp2(1:100);
        deckprops = [tmp1,tmp2];

        tmp1 = conv(pchoice(:,1),kern); tmp1 = tmp1(1:100);
        tmp2 = conv(pchoice(:,2),kern); tmp2 = tmp2(1:100);
        smpchoice = [tmp1,tmp2];

        t=1:100;
        figure; plot(t,deckprops(:,1),'b-',t,smpchoice(:,1),'b-.',t,deckprops(:,2),'g-',t,smpchoice(:,2),'g-.');
        
        legend('good decks', 'predicted good', 'bad decks', 'predicted bad','Location','NW');
        LMneg
        
        fprintf(1,'\n'); fprintf(1,'\n');
        nms = {'a', 'b', 'p1', 'LMneg'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',a); fprintf(1,'%3.4f\t',b); fprintf(1,'%3.4f\t',p1);
        fprintf(1,'%3.4f\t',LMneg);fprintf(1,'\n');fprintf(1,'\n');
        
    case 'baseline'
        start_p1 = .25;
        start_p2 = .25;
        start_p3 = .25;

        p1 = zeros(4,4,4);
        p2 = zeros(4,4,4);
        p3 = zeros(4,4,4);
        p4 = zeros(4,4,4);

        f = @(phats) gambling_negLM(phats,data,actual_choices,[],[],'baseline')
        warning off
        for i = 1:length(start_p1),
            for j = 1:length(start_p2),
                for k = 1:length(start_p3),
                    start_p4 = 1-(start_p1(i)+start_p2(j)+start_p3(k));
                    phats = [start_p1(i); start_p2(j); start_p3(k); start_p4];
                    phats = phats./sum(phats);
                    [phats,fval(i,j,k)] = fminsearch(f,phats);
                    p1(i,j,k) = phats(1); p2(i,j,k) = phats(2); p3(i,j,k) = phats(3); p4 = phats(4);
                end
            end
        end
        warning on
        
        [x,y,z] = matrix_min(fval);
        p1 = p1(x,y,z);
        p2 = p2(x,y,z);
        p3 = p3(x,y,z);
        p4 = 1-(p1+p2+p3);
        LMneg = fval(x,y,z);
        
        fprintf(1,'\n'); fprintf(1,'\n');
        nms = {'p1', 'p2', 'p3','p4', 'LMneg'};
        fprintf(1,'%s\t',nms{:}); fprintf(1,'\n');
        fprintf(1,'%3.4f\t',p1); fprintf(1,'%3.4f\t',p2); fprintf(1,'%3.4f\t',p3); fprintf(1,'%3.4f\t',p4);
        fprintf(1,'%3.4f\t',LMneg);fprintf(1,'\n');fprintf(1,'\n');
        
        
end



