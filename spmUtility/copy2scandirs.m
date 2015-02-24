% move all task imgs to their own directory

subs = [2 4 5 6 8 9 10 11 12];
drive = [1 1 1 1 2 2 4 4 4];


for JJ = subs
   snum = num2str(JJ);
   dnum = num2str(drive(subs == JJ));
   
   fmriDIR = ['/data' dnum '/intext'];
   str = ['cd ' fmriDIR '/sub' snum '/task'];
   eval(str)
   
   for KK = 1:8
      str = ['!mkdir scan' num2str(KK)'];
      eval(str)
   end
   
   !mv sra*_00[0-9][0-9]* scan1/
   !mv sra*_01[0-4][0-9]* scan1/
   !mv sra*_0150* scan1/
   
   !mv sra*_01[5-9][0-9]* scan2/
   !mv sra*_02[0-9][0-9]* scan2/
   !mv sra*_0300* scan2/
   
   !mv sra*_03[0-9][0-9]* scan3/
   !mv sra*_04[0-4][0-9]* scan3/
   !mv sra*_0450* scan3/
   
   !mv sra*_04[5-9][0-9]* scan4/
   !mv sra*_05[0-9][0-9]* scan4/
   !mv sra*_0600* scan4/
   
   !mv sra*_06[0-9][0-9]* scan5/
   !mv sra*_07[0-4][0-9]* scan5/
   !mv sra*_0750* scan5/
   
   !mv sra*_07[5-9][0-9]* scan6/
   !mv sra*_08[0-9][0-9]* scan6/
   !mv sra*_0900* scan6/
   
   !mv sra*_09[0-9][0-9]* scan7/
   !mv sra*_10[0-4][0-9]* scan7/
   !mv sra*_1050* scan7/
   
   !mv sra*_10[5-9][0-9]* scan8/
   !mv sra*_11[0-9][0-9]* scan8/
   !mv sra*_1200* scan8/   
   
   
   
   
end
