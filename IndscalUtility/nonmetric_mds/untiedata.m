function [ newseries , ntiedblocks , whichtiedblock ] = untiedata( ns , series );

whichtiedblock  = zeros( 1,ns );
ntiedblocks     = 0;
newseries       = 1;

for i=1:ns
   if (i>1)
      now = series( i );
      
      if (last==now)
         if (newseries==1)
            ntiedblocks = ntiedblocks + 1;
            newseries = 0;
         end   
         
         whichtiedblock( i   ) = ntiedblocks;
         whichtiedblock( i-1 ) = ntiedblocks;
      else
         newseries = 1;
      end
      
   end
   
   last = series( i );
end

order     = 1:ns;
newseries = 1:ns;

for i=1:ntiedblocks
   ii = find( whichtiedblock==i );
   
   av = mean(  order( ii ) );
   
   newseries( ii ) = av;
end
