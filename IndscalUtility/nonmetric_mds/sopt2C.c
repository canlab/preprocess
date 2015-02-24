#include "d:\matlab\extern\include\mex.h"
#include <math.h>

/*
f = sopt2C( Config , l , x , no , nt , ns , referindex , r , DO , ...
   			DS , DHS , effected , effecteds , oldsum1 , oldsum2 , minoption );
*/

void sopt2C( 
 		double *f,
        double *Config,
        int    l,
        double *x,
		int    no,
        int    nt,
		int    ns,
		double *referindex,
		double r,
		double *DO,
		double *DS,
		double *DHS,
		double *effected,
		double *effecteds,
		double oldsum1,
		double oldsum2,
		int    minoption    )
      
{
	int i,j,k,ind,effected_index1,effected_index2;
	double d1,value,value2,temp,newsum1,newsum2,newsum3;

	/* update the configuration with new positions of object l */
	for (i=0; i<nt; i++)
	{
		*( Config + (l-1) + i * no )= *( x + i ); 
	}

	/* update the distances where object l is involved
	 these are the distances that change due to optimization */
    
    if (r==2)
	{  
      if (l<no)
	  {
		i=l;
		for (j=i+1; j<=no; j++)
		{
			ind = (int) *( referindex + (i-1) + (j-1) * no );
     
			d1 = (double) 0;
			for (k=0; k<nt; k++)
			{
				value =  fabs( *( Config + (i-1) + k * no ) - *( Config + (j-1) + k * no ) );
			    d1 += value * value; /* assuming r=2 !!!!!!!!!!!!!!!!!!!!!!!!!! */
			}
			
			d1 = sqrt( d1 );
			*( DO + (ind-1) ) = d1;
		}
	  }

	  if (l>1)
	  {
		  j=l;
		  for (i=1; i<=j-1; i++)
		  {
			ind = (int) *( referindex + (i-1) + (j-1) * no );
			
			d1 = (double) 0;
			for (k=0; k<nt; k++)
			{
				value =  fabs( *( Config + (i-1) + k * no ) - *( Config + (j-1) + k * no ) );
			    d1 += value * value; /* assuming r=2 !!!!!!!!!!!!!!!!!!!!!!!!!! */
			}
			
			d1 = sqrt( d1 );
			*( DO + (ind-1) ) = d1;
		  }   
	  }
	} else
	{
	  /* other R values */
	  if (l<no)
	  {
		i=l;
		for (j=i+1; j<=no; j++)
		{
			ind = (int) *( referindex + (i-1) + (j-1) * no );
     
			d1 = (double) 0;
			for (k=0; k<nt; k++)
			{
				value =  fabs( *( Config + (i-1) + k * no ) - *( Config + (j-1) + k * no ) );
			    d1 += pow( value , r );
			}
			

			d1 = pow( d1 , 1/r );
			*( DO + (ind-1) ) = d1;
		}
	  }

	  if (l>1)
	  {
		  j=l;
		  for (i=1; i<=j-1; i++)
		  {
			ind = (int) *( referindex + (i-1) + (j-1) * no );
			
			d1 = (double) 0;
			for (k=0; k<nt; k++)
			{
				value =  fabs( *( Config + (i-1) + k * no ) - *( Config + (j-1) + k * no ) );
			    d1 += pow( value , r );
			}
			
			d1 = pow( d1 , 1/r );
			*( DO + (ind-1) ) = d1;
		  }   
	  }
    }

      
	if (minoption==1)
    {
	  newsum1 = 0;
	  newsum2 = 0;

	  for (i=0; i<(no-1); i++)
	  {
	    effected_index1 = (int) *( effecteds + (l-1) + i * no ); /* effected list */
		effected_index2 = (int) *( effected  + (l-1) + i * no ); /* effected */
		
		value = *( DO + (effected_index2-1) );
		*( DS + (effected_index1-1) ) = value;
		
		value2 = *( DHS + (effected_index1-1) );

		temp = value - value2;
		
		newsum1 += temp * temp; 
		newsum2 += value * value;	
	  } 
	
	  /*
	  effectedlist = effecteds( l,:);
	  DS( effectedlist ) = DO( effected( l,:) );
	

	  newsum1 = sum( ( DS(  effectedlist ) - ...
		               DHS( effectedlist )       ).^2 );
              
	  newsum2 = sum(   DS(  effectedlist ).^2 );
	
	  f = sqrt( (newsum1+oldsum1)/(newsum2+oldsum2) );     
	  */

	  *f = sqrt( (newsum1+oldsum1)/(newsum2+oldsum2) );
	} 
	
	   else /* minoption=2 */
   
	{
    
	  newsum3 = 0;
      /* calculate the mean distance DO in newsum3 */
	  for (i=0; i<ns; i++)
	  {
         value    = *( DO + i );
		 newsum3 += value;
	  }
	  newsum3 /= (double) ns;

	  newsum2 = 0;
      /* calculate the sum of squared deviations from distances DO */
	  for (i=0; i<ns; i++)
	  {
         value    = *( DO + i );
		 temp     = ( value - newsum3 );
		 newsum2 += temp * temp;
	  }
      
      newsum1 = 0;
	  for (i=0; i<(no-1); i++)
	  {
	    effected_index1 = (int) *( effecteds + (l-1) + i * no ); /* effected list */
		effected_index2 = (int) *( effected  + (l-1) + i * no ); /* effected */
		
		value = *( DO + (effected_index2-1) );
		*( DS + (effected_index1-1) ) = value;
		
		value2 = *( DHS + (effected_index1-1) );

		temp = value - value2;
		newsum1 += temp * temp; 
	  } 
	
	  /*
	  effectedlist = effecteds( l,:);
	  DS( effectedlist ) = DO( effected( l,:) );
	

	  newsum1 = sum( ( DS(  effectedlist ) - ...
		               DHS( effectedlist )       ).^2 );
              
	  newsum2 = sum(   ( DS( : ) - mean( DS( : ) ) .^2 );
	
	  f = sqrt( (newsum1+oldsum1)/(newsum2) );     
	  */

	  *f = sqrt( (newsum1+oldsum1)/(newsum2) );

	}

}

void mexFunction(  int nlhs, mxArray *plhs[],
				   int nrhs, const mxArray *prhs[] )

{
	double *f;
	double *Config;
	int    l;
	double *x;
	int    no;
    int    nt;
	int    ns;
	int    minoption;
	double *referindex;
	double r;
	double *DO;
	double *DS;
	double *DHS;
	double *effected;
	double *effecteds;
	double oldsum1;
	double oldsum2;

   if (nrhs != 16)
   {
      mexErrMsgTxt( "Only sixteen input argument allowed." );
   }
      else if (nlhs != 1) 
   {
      mexErrMsgTxt( "One output arguments allowed." );
   }
   	
	Config		=			mxGetPr(		prhs[ 0 ] );
   	l			= (int)		mxGetScalar(	prhs[ 1 ] );
   	x			=			mxGetPr(		prhs[ 2 ] );
	no			= (int)		mxGetScalar(	prhs[ 3 ] );
	nt			= (int)		mxGetScalar(	prhs[ 4 ] );
	ns			= (int)		mxGetScalar(	prhs[ 5 ] );
	referindex	=			mxGetPr(		prhs[ 6 ] );
	r			= (double)	mxGetScalar(	prhs[ 7 ] );
	DO			=			mxGetPr(		prhs[ 8 ] );
	DS			=			mxGetPr(		prhs[ 9 ] );
	DHS			=			mxGetPr(		prhs[ 10] );
	effected	=			mxGetPr(		prhs[ 11] );
	effecteds	=			mxGetPr(		prhs[ 12] );
	oldsum1		= (double)	mxGetScalar(	prhs[ 13] );
	oldsum2		= (double)	mxGetScalar(	prhs[ 14] );
    minoption	= (int)		mxGetScalar(	prhs[ 15] );

	plhs[ 0 ]   = mxCreateDoubleMatrix( 1,1, mxREAL );
	f			=			mxGetPr(		plhs[ 0 ] );
   
   /* call the subroutine */
   sopt2C( f,Config,l,x,no,nt,ns,referindex,r,DO,DS,DHS,effected,effecteds,oldsum1,oldsum2,minoption );

}
