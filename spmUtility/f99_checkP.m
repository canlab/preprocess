function f99_checkP(DIR,EXP,mat)
%--------------------------------------------------------------------
% gaat na of de files in de P matrix effectief bestaan op schijf
% er wordt enkel of de files gecheckt als flag = 1
%--------------------------------------------------------------------
global fmriTEST;


if fmriTEST
  for i=1:size(DIR,1)
    if f99_exist(deblank(DIR(i,:)),'.') == 0
      msg = ['    *** WARNING files ' DIR(i,:) filesep EXP(i,:) ' do not (YET) exist'];
      disp(msg)
    else
      [fi , di] = spm_list_files(DIR(i,:),EXP(i,:));
      if fi > 0
        msg = ['    ' num2str(size(fi,1)) ' ' EXP(i,:) ' files do exist'];
      else
        msg = ['    ' EXP(i,:) ' files do NOT (YET) exist'];
      end
      disp(msg);
      if mat
        EXP2 = [deblank(spm_str_manip(EXP(i,:),'r')) '.mat'];
        [fi , di] = spm_list_files(DIR(i,:),EXP2);
        if fi > 0
          msg = ['    ' num2str(size(fi,1)) ' ' EXP2 ' files do also exist'];  
        else
          msg = ['    ' EXP2 ' files do NOT (YET) exist'];
        end
        disp(msg);
      end
    end
  end
end

