function LogFile = f99_log(LogFile,LogText)
%==========================================
%-Open LogFile
%----------------------------------------------------------------------------
global fmriTEST
LogText = strrep(LogText,'\','\\');

if ~fmriTEST
  [fid,message] = fopen(LogFile,'a');
  if fid==-1
    fprintf('%cfmri96_log error : tNo logging!\n\t%s\n\n',7,message);
    return
  end
  
  
  %-Write log
  %----------------------------------------------------------------------------
  if ~isempty(LogText) 
     fprintf([LogText,'\n']);
     fprintf(fid,[LogText,'\n']);
  end
  
  %-Close log
  %----------------------------------------------------------------------------
  status = fclose(fid);
  if fid==-1, fprintf('%cfmri96_log : error closing file, fid = %d\n',7,fid); end
  
else
   if ~isempty(LogText)
    fprintf([LogText,'\n']);
  end
end
