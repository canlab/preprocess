clear events
ss = [1 2 4 6 8 9 10 11 12];
drive = [1 1 1 1 2 2 4 4 4];

load /data1/intext/model8times
for i = 1:size(ss,2)
   myevents = [];
   for j = 1:8
      eval(['myevents = [myevents; sub' num2str(ss(i)) 'model8.offset{j}''];']);
   end
   for j = 1:6
      events{i}{j} = myevents(:,j);
   end
   events_subjects(i) = ss(i);
   
   mytaskdirs{i} = fullfile(['/data' num2str(drive(i))],'intext',['sub' num2str(ss(i))],'task')
   myresultsdirs{i} = fullfile(['/data' num2str(drive(i))],'intext','RESULTS','model13',['sub' num2str(ss(i))])

end
