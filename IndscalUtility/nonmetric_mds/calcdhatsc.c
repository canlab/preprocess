#include "d:\matlab\extern\include\mex.h"
#include <math.h>
 
calcdhatsc( double *DH,
            int    ns,
			
            double *DO,
			
			double *inblock,
			double *ave,
			double *sum,
			double *nmem     )
{
  int     i,j,nblocks,blockactive,ups,downs,upa,reference,nmem1,nmem2,whatblock;
  double  value,av,sum1,sum2;

  /* Kruskal's upward-downward algorithm in C */

  /* start with the finest partition possible; every point is in one block */
  nblocks = ns;
  for (i=0; i<ns; i++)
  {
	 /* the inblock has a different interpretation from the MATLAB implementation:
	    every entry i in this array is filled with the block number that it is a member of */

	 *(inblock + i) = (double) i;
     value          = *( DO + i);
     /* square the distances */
	 
	 /*
	 *(ave	   + i) = value * value;
	 *(sum     + i) = value * value;
	 */

	 /* DO NOT square the distances */ 
	 *(ave	   + i) = value;
	 *(sum     + i) = value;

     *(nmem    + i) = (double) 1;
  }

  /* do the iterative process until the highest-block is up and down satisfied */
  blockactive  = 0;	/* in c, start with zero index */

  ups          = 0;
  downs        = 1;
  upa          = 1;

  while ((blockactive+1 < nblocks) || (ups==0) || (downs==0))
  { 
    /* is the current block up or downactive? */
    if (upa==1)
	{
       /* is the current block upsatisfied ? */
       ups = 1; /* assume it is satisfied */
       if (blockactive+1 < nblocks)
       {  
		  if ( *( ave+blockactive+1 ) < *( ave+blockactive ) )
          {
			 ups = 0;
          }
       }
      
       /* if the block is upsatisfied then the active block becomes downactive */
       if (ups == 1)
	   {
          upa = 0;
       }
		  else
	   {
	      /* if not, the active block is joined with the next higher block  */
          
          /* MATLAB code:
              members1 = B(blockactive).mem;
		      members2 = B(blockactive+1).mem;
              B(blockactive).mem  = union( members1 , members2 );

              the blocks with numbers "blockactive" and "blockactive+1" are joined
              in the block with number "blockactive";
              this means that all references to "blockactive+1 should now be changed to 
              blockactive
		  */

		  for (i=0; i<ns; i++)
          {
			  reference = (int) *( inblock+i );

              /* decrement all references higher than blockactive+1 */ 
              if (reference>=blockactive+1)
              {
				 *( inblock+i ) = (double) reference - 1;
              }
          }

          sum1     = *( sum+blockactive  );
          nmem1    = *( nmem+blockactive );
     
          sum2     = *( sum+blockactive+1);
          nmem2    = *( nmem+blockactive+1);
      
          *(sum+blockactive)  = sum1 + sum2;
          *(nmem+blockactive) = nmem1 + nmem2;
          *(ave+blockactive)  = (sum1 + sum2) / (nmem1 + nmem2);
      
          /* and now move all the blocks one downward */
          if (blockactive<=nblocks-3)
          {
             for (j=blockactive+1; j<=nblocks-2; j++)
             {
			   *( sum+j  ) = *( sum+j+1  );
			   *( ave+j  ) = *( ave+j+1  );
               *( nmem+j ) = *( nmem+j+1 );
             }
          }
         
          nblocks--;
          
          /* and the new, larger block becomes downactive */
          upa = 0;
		}
    }   
       else /* current block is downactive */
    {
       /* is the current block down satisfied ? */
       downs = 1; /* assume it is satisfied */
       if (blockactive > 0)
	   {
          if ( *( ave+blockactive-1 ) > *( ave+blockactive ) )
          {
             downs = 0;
          }
       }
      
       /* if the block is downsatisfied then the active block becomes upactive */
       if ( downs == 1 )
       {
          upa = 1;
       }
          else
       {
          /* if not, the active block is joined with the next lower block */
          
		  /* MATLAB code:
              members1 = B(blockactive).mem;
              members2 = B(blockactive-1).mem;
		      B(blockactive-1).mem  = union( members1 , members2 );

              decrement all references >= blockactive 
		  */

          for (i=0; i<ns; i++)
          {
			 reference = (int) *( inblock+i );

             if (reference >= blockactive)
             {
                *( inblock+i) = (double) reference-1;
             }
		  }

          sum1     = *(sum+blockactive);
          nmem1    = *(nmem+blockactive);
          
          sum2     = *(sum+blockactive-1);
          nmem2    = *(nmem+blockactive-1);
          
          *(sum+blockactive-1)  = sum1 + sum2;
          *(nmem+blockactive-1) = nmem1 + nmem2;
          *(ave+blockactive-1)  = (sum1 + sum2) / (nmem1 + nmem2);
      
          /* and now move all the blocks, including the active block one downward */
          if (blockactive<=nblocks-2)
          {
             for (j=blockactive; j<=nblocks-2; j++)
             {
                *( sum+j  ) = *( sum+j+1  );
 			    *( ave+j  ) = *( ave+j+1  );
                *( nmem+j ) = *( nmem+j+1 );
          	 }
          }
         
          blockactive--;
          nblocks--;
          
          /* and the new, larger block becomes upactive */
          upa = 1;
       }
    } 
   
    /* check whether the active block is simultaneously up and down satisfied */
    ups = 1; /* assume it is satisfied */
    if (blockactive+1 < nblocks)
    {
       if ( *( ave+blockactive+1 ) < *( ave+blockactive ) )
       {
          ups = 0;
       }
    }
   
    downs = 1; /* assume it is satisfied */
    if (blockactive > 0)
    {
       if ( *(ave+blockactive-1) > *( ave+blockactive ) )
       {
          downs = 0;
       }
    }
   
    if ((ups==1) && (downs==1))
    {
       /* both up and down-satisfied */
       if (blockactive+1 < nblocks)
       {
          /* transfer to next block */
          blockactive++;
          upa = 1;
       }
    }
  }   


  /* now, fill the DH-array with the block averages */
  for (i=0; i<ns; i++)
  {
    whatblock = (int) *( inblock+i );
    av = *( ave+whatblock );
    *( DH+i ) = av; 
  }

  /* and take the square root of the distances */
  /*
  for (i=0; i<ns; i++)
  {
    *( DH+i ) = sqrt( *( DH+i ) );
  }
  */
}


void mexFunction(  int nlhs, mxArray *plhs[],
				   int nrhs, const mxArray *prhs[] )

{
   double *DH;
   int    ns;
			
   double *DO;
			
   double *inblock;
   double *ave;
   double *sum;
   double *nmem;
   
   int    m,n;
   
   /* the input arguments are:
   0 = ns
   1 = DO
   
   the output arguments are
   0:  DH
   1:  inblock
   2:  ave
   3:  sum
   4:  nmem
   */

   if (nrhs != 2)
   {
      mexErrMsgTxt( "Only two input argument allowed." );
   }
      else if (nlhs != 5) 
   {
      mexErrMsgTxt( "Only five output argument allowed." );
   }
   
   /* the input ns must be a scalar */
   m = mxGetM( prhs[ 0 ] );
   n = mxGetN( prhs[ 0 ] );   
   if (!mxIsDouble( prhs[0]) || mxIsComplex(prhs[0]) || !(m ==1 && n==1)) mexErrMsgTxt("Input ns must be a scalar." );

   ns = mxGetScalar( prhs[ 0 ] );

   /* check the dimensions of the DO array */
   m = mxGetM( prhs[ 1 ] );
   n = mxGetN( prhs[ 1 ] );   
 
   if ((m != ns) || (n != 1)) mexErrMsgTxt("The DO Array has the wrong dimensions." );
   DO = mxGetPr( prhs[ 1 ] );

   /* the output arguments are
   0:  DH
   1:  inblock
   2:  ave
   3:  sum
   4:  nmem

   /* Create matrixes for the two return arguments. */   
   plhs[ 0 ] = mxCreateDoubleMatrix( ns,1, mxREAL );
   plhs[ 1 ] = mxCreateDoubleMatrix( ns,1, mxREAL );
   plhs[ 2 ] = mxCreateDoubleMatrix( ns,1, mxREAL );
   plhs[ 3 ] = mxCreateDoubleMatrix( ns,1, mxREAL );
   plhs[ 4 ] = mxCreateDoubleMatrix( ns,1, mxREAL );
   
   DH      = mxGetPr( plhs[ 0 ] );
   inblock = mxGetPr( plhs[ 1 ] );
   ave     = mxGetPr( plhs[ 2 ] );
   sum     = mxGetPr( plhs[ 3 ] );
   nmem    = mxGetPr( plhs[ 4 ] );
   
   /* call the subroutine */
   calcdhatsc( DH,ns,DO,inblock,ave,sum,nmem );
}

