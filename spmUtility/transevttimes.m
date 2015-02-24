% transforms event times into c cell array, with onsets for each condition.
% for use with the output of tor's scripts for getting event times from behavioral data.

clear c
b = [];

for i = 1:numscans

	% adjust by padding with -1's if num rows differs.
  	if size(b,1) < size(a{i}',1)                 
		b = [b; -1*ones(size(a{i}',1)-size(b,1),size(b,2))]; 
	elseif size(b,1) > size(a{i}',1)  
		a{i} = [a{i} -1*ones(size(b,1)-size(a{i},2),size(a{i},1))'];
	end

   	b = [b a{i}'];
end

for i = 1:size(b,2)
   c{i} = b(b(:,i) > -.5,i);
end
clear a
