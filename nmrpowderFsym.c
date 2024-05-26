#include <math.h>
#include <stdio.h>
#include "mex.h"

/*
 * nmrpowderF.c - Fast simulation of the powder spectrum for chemical/knight shift anisotropy
 * implementing ZCW powder averaging method
 *
 * This is a MEX-file for MATLAB.
 * Copyright 13.03.2010 - 26.04.2012   Anton Potocnik @ IJS F5  
 */

/* $Revision: 2.0 $ */
double mod(double a, double b) {
    int r;
    r = (int)(a/b);
    return a - r;
}

//********* ZCW ROUTINES ***********
int getNumberZCW(int M)
//returns the number of ZCW angles for the given integer M2,3,4,. . .
{
    int j,gM=5,gMminus1=3;
    int sum=5;
    for(j=0;j<=M;j++) {
        sum = gM + gMminus1;
        gMminus1=gM;
        gM=sum;
    }
    return sum;
}

void nmrpowder(double *spc, double *f, mwSize M, double fmax, mwSize NOP, double A, double LB, double Kiso, double Kax, double Kasym)
{
  int i,j,g2,N;
  double Fs, fc, fi, costh, sinth, cos2fi, maxi=-100000, UNI, c[3];

    for (i=0;i<NOP;i++)
	{
		*(f+i) = -fmax + (double)i/(NOP-1)*2.0*fmax;
		*(spc+i) = 0;
	}
	
    //c[0]=1.;c[1]=2.;c[2]=1.; // sphereType = “full”
    //c[0]=-1.;c[1]=1.;c[2]=1.; // sphereType = “hemi”
    c[0]=-1.;c[1]=1.;c[2]=4.; // sphereType = “oct”

    N = getNumberZCW(M); //total number of angles
    g2 = getNumberZCW(M-2);

	for (i=0; i<N; i++) {
            UNI = (double)i/N;
            fi = 2*3.141592*mod(UNI*g2,1)/c[2];
            costh = c[0]*(c[1]*mod(UNI,1)-1);
			cos2fi = cos(2*fi);
			sinth = sin(acos(costh));
			
			fc = Kiso + Kax*(3*costh*costh-1)/2 + Kasym*sinth*sinth*cos2fi;
			for (j=0; j<NOP; j++)
				//*(spc+j) += 0.6366198*LB/(4*(*(f+j)-fc)*(*(f+j)-fc) + LB*LB);
                *(spc+j) += 1/LB*exp(-(*(f+j)-fc)*(*(f+j)-fc)/2/LB/LB);
	}

    // Normalize
	for (i=0; i<NOP; i++)
		if (*(spc+i) > maxi) maxi = *(spc+i);

	for (i=0; i<NOP; i++) *(spc+i) /= maxi/A;
  
}





/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *spc,*f;
  double fmax, A, LB, Kiso, Kax, Kasym;
  mwSize mrows, ncols, N, NOP, m, k;
  
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
  fmax = mxGetScalar(prhs[1]);
  NOP = mxGetScalar(prhs[2]);
  A = mxGetScalar(prhs[3]);
  LB = mxGetScalar(prhs[4]);
  Kiso = mxGetScalar(prhs[5]);
  Kax = mxGetScalar(prhs[6]);
  Kasym = mxGetScalar(prhs[7]);

  mrows = NOP;
  ncols = 1;

  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(mrows,ncols, mxREAL); 
  /*  create a C pointer to a copy of the output matrix */
  spc = mxGetPr(plhs[0]);
  f = mxGetPr(plhs[1]);
  
  /*  call the C subroutine*/ 
	nmrpowder(spc, f, N, fmax, NOP, A, LB, Kiso, Kax, Kasym);
}
