function E = f99_exist(b,f)
%-----------------------------
% b - zoek_directory
% f - dir of file te zoeken
%-----------------------------
[fi , di] = spm_list_files(b,'*');
for i = 1:size(fi,1)
  %fi(i,:)
  if strcmp(deblank(fi(i,:)),f)
    E=1; return;
  end
end
for i = 1:size(di,1)
  %di(i,:)
  if strcmp(deblank(di(i,:)),f)
    E=2;return;
  end
end
E=0;
return;
