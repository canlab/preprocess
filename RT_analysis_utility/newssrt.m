function ssrt = newssrt(tabledat,stopperc,all_stoptimes)

 
nblocks = length(tabledat);
startblk = 1

for i = startblk:nblocks
  rt = tabledat{i}{3}(tabledat{i}{2} ==1);
  x = prctile(rt,100*stopperc(i,2));
  new_SSRT(i) = x - all_stoptimes(i);
end
    


nms ={'SSRT1' 'SSRT2' 'SSRT3'};
fprintf(1,'%s\t',nms{:});
fprintf(1,'\n');
fprintf(1,'%3.4f\t',new_SSRT);
    
