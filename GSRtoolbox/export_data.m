function export_data(mainfig)
   plotdata = get(mainfig, 'UserData');
   fs = plotdata.fs;
   
   filename = 'myfile.txt';
   
   fid = fopen(filename, 'w');
   dl = '\t';
   
   fprintf(fid,['subject' dl 'session' dl 'trial' dl 'trialLabel' dl 'responseOnset' dl 'peakOnset' dl 'responseDuration' dl 'responseArea' dl 'responseHeight' '\n']);
   fclose(fid);
   outlines = [];
   for subj = 1:length(plotdata.trials)
       for sess = 1:length(plotdata.trials{subj})
           for trial = 1:size(plotdata.trials{subj}{sess},1)
               triallabel = plotdata.trials{subj}{sess}(trial,3);
               gsrnum = 0;
               while gsrnum >= 0
               
                   gsrnum = search_for_GSR(plotdata.SCR.location{subj}{sess}, plotdata.trials{subj}{sess}(trial,1),plotdata.trials{subj}{sess}(trial,2),gsrnum);
                   
                   if gsrnum >= 0
                        
                        onset = (plotdata.SCR.location{subj}{sess}(gsrnum, 1) - plotdata.trials{subj}{sess}(trial,1))/fs;
                        peak = (plotdata.SCR.location{subj}{sess}(gsrnum, 2) - plotdata.trials{subj}{sess}(trial,1))/fs;
                        len = (plotdata.SCR.location{subj}{sess}(gsrnum, 3) - plotdata.trials{subj}{sess}(trial,1))/fs;
                        startheight = plotdata.signal{subj}{sess}(plotdata.SCR.location{subj}{sess}(gsrnum, 1));
                        gsrarea = sum(plotdata.signal{subj}{sess}(plotdata.SCR.location{subj}{sess}(gsrnum, 1):plotdata.SCR.location{subj}{sess}(gsrnum, 1))-startheight)/fs;
                        gsrheight = (plotdata.SCR.height{subj}{sess}(gsrnum));
                        
                        outlines = [outlines; subj sess trial triallabel onset peak len gsrarea gsrheight];
                        
                       
                   end;
                   
                   
               end;
               
               
           end;
       end;
   end;
       
   dlmwrite(filename, outlines, '-append', 'delimiter', dl, 'precision', 10);
   
   
end





function gsrfound = search_for_GSR(rectpos, gstart, gend, glast)
    flag = 0;
    for i = (glast+1):length(rectpos)
        if ((rectpos(i,1)>gstart) && (rectpos(i,1) < gend)) || ((rectpos(i,3)>gstart) && (rectpos(i,3) < gend)) || ((rectpos(i,1) < gstart) && (rectpos(i,3) > gend))
            flag = 1;
            break;
        end;
    end;
    if flag
        gsrfound = i;
    else
        gsrfound = -1;
    end;
end
    
    