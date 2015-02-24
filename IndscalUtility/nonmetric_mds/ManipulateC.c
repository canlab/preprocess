#include "d:\matlab\extern\include\mex.h"
#include <math.h>
  
/*
	manipulate( c1,c2,no , DistsSub );
*/

void ManipulateC(
 
		double *NewSub,
 		double *DistsSub,
		double *c1,
		double *c2,
		int    no
						)

{
	double minimum,value,cc1,extra,maxextra;
	long index;
	int i,j,l;

	/* manipulate the subjects' rating so that the triangular inequality and
	positivity holds

    see if positivity holds; first, find the minimum value
	*/

	minimum = (double) 1000000000;
	index   = 0;
	for (i=0; i<no; i++)
	{
		for (j=0; j<no; j++)
		{
			*( NewSub + index ) = *( DistsSub + index );

			if (i!=j)
			{
				value = *( NewSub + index );
				if (value < minimum) minimum = value;
			}

			index++;
		}
	}

	/* now add this value if the minimum is smaller than zero */
	cc1 = (double) 0;
	if (minimum < 0)
	{
		index = 0;
		cc1 = -minimum;
		for (i=0; i<no; i++)
		{
			for (j=0; j<no; j++)
			{
				if (i!=j) *( NewSub + index ) = *( NewSub + index ) - minimum;
				index++;
			}
		}
	}

	/*	assume that all diagonal entries are zero and that the matrix is symmetric
		find the entry that violates the triangular inequality the most

	*/

	maxextra = (double) 0;
	for (i=0; i<no; i++)
	{
		for (j=0; j<no; j++)
		{
			for (l=0; l<no; l++)
			{
				if ((i!=j) && (j!=l) && (i!=l))
            	{
					/* extra = DistsSub( i,l ) - ( DistsSub(i,j) + DistsSub( j,l) ); 
					*/

					extra = *( NewSub + i + l*no ) - 
							( *( NewSub + i + j*no ) + *( NewSub + j + l*no ));

					if (extra > maxextra) maxextra = extra;
				}
			}
		}
	}


	if (maxextra > 0)
	{
		index = 0;
		for (i=0; i<no; i++)
		{
			for (j=0; j<no; j++)
			{
				if (i!=j) *( NewSub + index ) = *( NewSub + index ) + maxextra;
				index++;  
			}
		}
	}

	*c1 = cc1;
    *c2 = maxextra;
}	

/*
	manipulate( c1,c2, no , DistsSub );
*/

void mexFunction(  int nlhs, mxArray *plhs[],
				   int nrhs, const mxArray *prhs[] )

{
	double *DistsSub;
	double *NewSub;
	double *c1;
	double *c2;
	int    no;		

   if (nrhs != 2)
   {
      mexErrMsgTxt( "Only two input argument allowed." );
   }
      else if (nlhs != 3) 
   {
      mexErrMsgTxt( "Only three output arguments allowed." );
   }
   	
	 
	DistsSub	=			mxGetPr(		prhs[ 0 ] );
	no			= (int)		mxGetScalar(	prhs[ 1 ] ); 

	plhs[ 0 ]   = mxCreateDoubleMatrix( no,no, mxREAL );
	NewSub		=			mxGetPr(		plhs[ 0 ] );

	plhs[ 1 ]   = mxCreateDoubleMatrix( 1,1, mxREAL );
	c1			=			mxGetPr(		plhs[ 1 ] );

	plhs[ 2 ]   = mxCreateDoubleMatrix( 1,1, mxREAL );
	c2			=			mxGetPr(		plhs[ 2 ] );

   /* call the subroutine */	
	ManipulateC( NewSub,DistsSub,c1,c2,no ); 
}

