function copyimgs2scanTESTdirs(fmriDIR)
   
   str = ['cd ' fmriDIR '/task'];
   eval(str)
   
   for KK = 1:8
      str = ['!mkdir scanTEST' num2str(KK)'];
      eval(str)
   end
   
   !mv *vol*_00[0-9][0-9]* scanTEST1/
   !mv *vol*_01[0-4][0-9]* scanTEST1/
   !mv *vol*_0150* scanTEST1/
   
   !mv *vol*_01[5-9][0-9]* scanTEST2/
   !mv *vol*_02[0-9][0-9]* scanTEST2/
   !mv *vol*_0300* scanTEST2/
   
   !mv *vol*_03[0-9][0-9]* scanTEST3/
   !mv *vol*_04[0-4][0-9]* scanTEST3/
   !mv *vol*_0450* scanTEST3/
   
   !mv *vol*_04[5-9][0-9]* scanTEST4/
   !mv *vol*_05[0-9][0-9]* scanTEST4/
   !mv *vol*_0600* scanTEST4/
   
   !mv *vol*_06[0-9][0-9]* scanTEST5/
   !mv *vol*_07[0-4][0-9]* scanTEST5/
   !mv *vol*_0750* scanTEST5/
   
   !mv *vol*_07[5-9][0-9]* scanTEST6/
   !mv *vol*_08[0-9][0-9]* scanTEST6/
   !mv *vol*_0900* scanTEST6/
   
   !mv *vol*_09[0-9][0-9]* scanTEST7/
   !mv *vol*_10[0-4][0-9]* scanTEST7/
   !mv *vol*_1050* scanTEST7/
   
   !mv *vol*_10[5-9][0-9]* scanTEST8/
   !mv *vol*_11[0-9][0-9]* scanTEST8/
   !mv *vol*_1200* scanTEST8/   
  return 
