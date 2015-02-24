
nsi = 1;
clear T, clear P, clear leg

for numsearch = 1:1000:10000

        dfi = 1;
    
        for df = 5:5:30
    
            [T(nsi,dfi), P(nsi,dfi)] = tor_bonf(.05,numsearch,df);
              
            
            if nsi == 1
                   leg{dfi} = [num2str(df) ' df'];
            end
            
            dfi = dfi + 1;
               
        end
        
        xlab(nsi) = numsearch;
        nsi = nsi + 1;
       
        
 end
    
 
 plot(T)
 set(gcf,'Color','w')
 xlabel('Number of voxels searched','FontSize',14)
 ylabel('Critical t value','FontSize',14)
 set(gca,'XTick',1:length(xlab));
 set(gca,'XTickLabel',xlab)
 legend(leg)