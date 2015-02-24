#include "d:\matlab\extern\include\mex.h"
#include <math.h>

/*
[ oldsum1,oldsum2 ] = NotEffectedC2( noteffecteds , DS , DHS , l , ns , no , minoption );
*/   
  
void NotEffectedC2(
 
 		double *oldsum1,
		double *oldsum2,
		double *noteffecteds,
		double *DS,
		double *DHS,
		int    l,
		int    ns,
		int    no,
		int    minoption )

{
	int i,noteffected_index;
	double sum1,sum2,DS_value,DHS_value,temp;
	
	/*
	noteffecteds = zeros( no , ns - no + 1 );

    noteffectedlist = noteffecteds( l,:);

    oldsum1 = sum( ( DS(  noteffectedlist ) - ...
                     DHS( noteffectedlist )       ).^2 );         
    oldsum2 = sum(   DS(  noteffectedlist ).^2 );
	*/

	sum1 = 0;
	sum2 = 0;

	for (i=0; i<ns-no+1; i++)
	{
		noteffected_index = (int) *( noteffecteds + (l-1) + i*no );

		DS_value	= *( DS	 + (noteffected_index-1) );
		DHS_value	= *( DHS + (noteffected_index-1) );
		
		temp = DS_value - DHS_value;

		sum1 += temp * temp;		 
		sum2 += DS_value * DS_value;
	}	

	*oldsum1 = sum1;
	*oldsum2 = sum2;
}

void mexFunction(  int nlhs, mxArray *plhs[],
				   int nrhs, const mxArray *prhs[] )

{
	double *oldsum1;
	double *oldsum2;
	double *noteffecteds;
	double *DS;
	double *DHS;
	int    l;
	int    ns;
	int    no;
	int    minoption;

   if (nrhs != 7)
   {
      mexErrMsgTxt( "Only 7 input argument allowed." );
   }
      else if (nlhs != 2) 
   {
      mexErrMsgTxt( "Two output arguments allowed." );
   }
   	
	noteffecteds=			mxGetPr(		prhs[ 0 ] );
   	DS			=			mxGetPr(		prhs[ 1 ] );
	DHS			=			mxGetPr(		prhs[ 2 ] );

	l			= (int)		mxGetScalar(	prhs[ 3 ] );
	ns			= (int)		mxGetScalar(	prhs[ 4 ] );
	no			= (int)		mxGetScalar(	prhs[ 5 ] );
	minoption	= (int)		mxGetScalar(	prhs[ 6 ] );

	plhs[ 0 ]   = mxCreateDoubleMatrix( 1,1, mxREAL );
	oldsum1		=			mxGetPr(		plhs[ 0 ] );

	plhs[ 1 ]   = mxCreateDoubleMatrix( 1,1, mxREAL );
	oldsum2		=			mxGetPr(		plhs[ 1 ] );

   /* call the subroutine */
	NotEffectedC2( oldsum1,oldsum2,noteffecteds,DS,DHS,l,ns,no );
}
