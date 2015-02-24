function plot_physio(physio, subj, st, en)
% function plot_physio(physio, subj, [st, en])

hr_samprate = 100;

if(~exist('st') || isempty(st))
    st = 1;
end
if(~exist('en') || isempty(en))
    en = length(physio.(subj).scanner_pulse);
end

secs = (1:length(physio.(subj).scanner_pulse))./100;
secs = secs(st:en);

if(isfield(physio.(subj), 'hr'))
    plot_hr = 1;
    num_rows = 8;
else
    plot_hr = 0;
    num_rows = 6;
end

figure('Color',[0.9 0.9 0.9]);

subplot(num_rows,1,1);
plot(secs,physio.(subj).scanner_pulse(st:en), 'k');
ylabel({'Scanner';'Pulse'});
title(subj, 'FontSize', 20);

subplot(num_rows,1,2);
plot(secs,physio.(subj).r1(st:en),'g');
ylabel({'Resp.';'Chest'});
subplot(num_rows,1,3);
plot(secs,physio.(subj).r2(st:en),'g');
ylabel({'Resp.';'Diaphragm'});

subplot(num_rows,1,4);
plot(secs,physio.(subj).rawgsr(st:en),'c');
ylabel({'Raw';'GSR'});
subplot(num_rows,1,5);
plot(secs,physio.(subj).gsr(st:en),'b');
ylabel({'Filt.';'GSR'});

subplot(num_rows,1,6);
plot(secs,physio.(subj).pulse(st:en),'r');
ylabel('Pulse');

if(plot_hr)
    hr_st = round(st / 100);
    hr_en = round(en / 100);
    hr_secs = hr_st:hr_en;
    
    subplot(num_rows,1,7);
    plot(hr_secs, physio.(subj).hr(hr_st:hr_en),'r');
    ylabel('HR');
    subplot(num_rows,1,8);
    plot(hr_secs,physio.(subj).hrs(hr_st:hr_en),'r');
    ylabel({'Smoothed';'HR'});
end

end