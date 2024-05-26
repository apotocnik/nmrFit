#include <math.h>
//#include <stdio.h>
#include "mex.h"

/*
 * nmrpowder.c - simulation of the powder average for chemical shift anisotropies
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix 
 *
 * This is a MEX-file for MATLAB.
 * Copyright 13.03.2010 - 16.03.2010   Anton Potocnik @ IJS F5  
 */

/* $Revision: 2.0 $ */

void nmrpowder(double *spc, double *f, mwSize N, double DW, mwSize TD, double A, double LB, double Kxx, double Kyy, double Kzz)
{
  int i,j,k;
  double Fs, ni_ref, fc, fi, costh, sinth, cosfi, sinfi, maxi=-100000;
  
	LB = LB*1e3;
	DW = DW*1e-6;
	Fs = 1/DW;
	ni_ref = 10000000;
	Kxx = Kxx*1e3 + ni_ref;
	Kyy = Kyy*1e3 + ni_ref;
	Kzz = Kzz*1e3 + ni_ref;

	/*N = (int)(sqrt(N));*/
	

	for (i=0;i<TD;i++)
	{
		*(f+i) = -Fs/2 + Fs/(TD-1)*i - Fs/TD/2;
		*(spc+i) = 0;
	}
	
	for (i=0; i<N; i++) {   /* fi */
		for (j=0; j<N; j++) {   /* th */
			fi = 6.283185*i/(N-1);
			cosfi = cos(fi);
			sinfi = sin(fi);
			costh = 2.0*j/(N-1) - 1;
			sinth = sin(acos(costh));
			
			fc = sqrt(Kxx*sinth*cosfi*Kxx*sinth*cosfi + Kyy*sinth*sinfi*Kyy*sinth*sinfi + Kzz*costh*Kzz*costh) -ni_ref;

			for (k=0; k<TD; k++)
				*(spc+k) += 0.6366198*LB/(4*(*(f+k)-fc)*(*(f+k)-fc) + LB*LB);
		}
	}

	for (k=0; k<TD; k++)
		if (*(spc+k) > maxi) maxi = *(spc+k);

	for (k=0; k<TD; k++) *(spc+k) /= maxi;
  
}





/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *spc,*f;
  double DW, A, LB, Kxx, Kyy, Kzz;
  mwSize mrows, ncols, N, TD, m, k;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=8) 
    mexErrMsgTxt("Eight inputs required.");
  if(nlhs!=2) 
    mexErrMsgTxt("Two output required.");
  
  /* check to make sure the first input argument is a vector 
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetN(prhs[0])*mxGetM(prhs[0])==1 ) {
    mexErrMsgTxt("Input x must be a vector.");
  }*/

  /*  create a pointer to the input matrix X 
  X = mxGetPr(prhs[0]);
  /  get the dimensions of the matrix input X /
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  if (mrows>=ncols) 
	m = mrows;
  else
	m = ncols;*/
	
  /*  get the scalars */
  N = mxGetScalar(prhs[0]); 
  DW = mxGetScalar(prhs[1]);
  TD = mxGetScalar(prhs[2]);
  A = mxGetScalar(prhs[3]);
  LB = mxGetScalar(prhs[4]);
  Kxx = mxGetScalar(prhs[5]);
  Kyy = mxGetScalar(prhs[6]);
  Kzz = mxGetScalar(prhs[7]);

  mrows = TD;
  ncols = 1;

  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(mrows,ncols, mxREAL); 
  /*  create a C pointer to a copy of the output matrix */
  spc = mxGetPr(plhs[0]);
  f = mxGetPr(plhs[1]);
  
  /*  call the C subroutine*/ 
	nmrpowder(spc, f, N, DW, TD, A, LB, Kxx, Kyy, Kzz);
}
