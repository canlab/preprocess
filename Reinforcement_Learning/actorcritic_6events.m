%actor critic

function [p_choice, v_hist, delta_hist, pa_hist, pb_hist]=actorcritic_6events(choices,reward,epsilon,gamma,beta,con);
% function []=actorcritic(choices,reward,epsilon,gamma,beta);
% function []=actorcritic(choices,reward,epsilon,gamma,beta);
% function [pb_hist]=actorcritic(choices,reward,epsilon,gamma,beta);

v_hist=[];
% v1_hist=[];
% va_hist=[];
% vb_hist=[];
% vstart_hist=[];
% va1_hist=[];
% va2_hist=[];
% 
% vb1_hist=[];
% vb2_hist=[];
% 
% delta_start_hist=[];
% delta1_hist=[];
% delta2_hist=[];

delta_hist=[];
pa_hist=[];
pb_hist=[];
% qa_hist=[];
% qb_hist=[];
% Aa_hist=[];
% Ab_hist=[];
p_choice=[];

vstart=0;
va=[0 0 0 0 0];
vb=[0 0 0 0 0];

qa=0;
qb=0;

pa=0.5;
pb=0.5;

for(i=1:length(choices))

   switch choices(i);

   case 1;
       %compute value predictions
       delta(5)=gamma*va(5)-va(4)+reward(i);
       delta(4)=gamma*va(4)-va(3);
       delta(3)=gamma*va(3)-va(2);
       delta(2)=gamma*va(2)-va(1);
       delta(1)=gamma*va(1)-vstart;

       v_hist = [v_hist vstart va(1) va(2) va(3) va(4)];

       va(1:4)=va(1:4)+epsilon*delta(2:5);
       vstart=pa*va(1)+pb*vb(1);
       
        if (con == 1 | con == 3)
            qa=qa+epsilon*(1-pa)*delta(1);
            qb=qb+epsilon*(-pa)*delta(1);
        else                                    % necessary only for correct probability of choice
            qa=qa+epsilon*(-pa)*delta(1);
            qb=qb+epsilon*(1-pa)*delta(1);
        end
            
       p_choice = [p_choice pa];

   case 2;
       %compute value predictions
       delta(5)=gamma*vb(5)-vb(4)+reward(i);
       delta(4)=gamma*vb(4)-vb(3);
       delta(3)=gamma*vb(3)-vb(2);
       delta(2)=gamma*vb(2)-vb(1);
       delta(1)=gamma*vb(1)-vstart;

       v_hist = [v_hist vstart vb(1) vb(2) vb(3) vb(4)];

       vb(1:4)=vb(1:4)+epsilon*delta(2:5);
       vstart=pa*va(1)+pb*vb(1);
       
        if (con == 1 | con == 3)
            qa=qa+epsilon*(-pa)*delta(1);
            qb=qb+epsilon*(1-pa)*delta(1);
        else                                    % necessary only for correct probability of choice
            qa=qa+epsilon*(1-pa)*delta(1);
            qb=qb+epsilon*(-pa)*delta(1);
        end
            
       p_choice = [p_choice pb];

    end;
    
      %softmax
       pa=exp(qa/beta)/(exp(qa/beta)+exp(qb/beta));
       pb=exp(qb/beta)/(exp(qa/beta)+exp(qb/beta));

%        delta_hist=[delta_hist delta_start delta(1:3)];
       pa_hist=[pa_hist pa];
       pb_hist=[pb_hist pb];
%        va_hist=[va_hist vstart va(1) va(2) va(3)];
%        vb_hist=[vb_hist vstart vb(1) vb(2) vb(3)];

%         delta_hist = [delta_hist delta_start delta(2)];
        delta_hist = [delta_hist delta(1) delta(5)];
%         delta_hist = [delta_hist delta(5)];
        
%         delta_start_hist = [delta_start_hist delta(1)];        
%         delta1_hist = [delta1_hist delta(1)];
%         delta2_hist = [delta2_hist delta(2)];
%         vstart_hist = [vstart_hist vstart];
%         va1_hist = [va1_hist va(1)];
%         va2_hist = [va2_hist va(2)];
%         vb1_hist = [vb1_hist vb(1)];
%         vb2_hist = [vb2_hist vb(2)];
%         qa_hist = [qa_hist qa];
%         qb_hist = [qb_hist qb]; 
%         Aa_hist = [Aa_hist qa-(mean([va(1) vb(1)]))];
%         Ab_hist = [Ab_hist qb-(mean([va(1) vb(1)]))];
                
end;