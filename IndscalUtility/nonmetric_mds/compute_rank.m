function RST = compute_rank( ns , tt1 , tt2 );

[ s1 , index1 ] = sort( tt1 );
[ s2 , index2 ] = sort( tt2 );
  
[ order1 , indexx1 ] = sort( index1 );
[ order2 , indexx2 ] = sort( index2 );

untied1 = untiedata( ns , s1 );
untied2 = untiedata( ns , s2 );

newindexx1 = untied1( indexx1 );
newindexx2 = untied2( indexx2 );

diff = (newindexx1 - newindexx2).^2;
sumd = sum( diff );
RST = 1 - 6*sumd / (ns*(ns^2-1));
