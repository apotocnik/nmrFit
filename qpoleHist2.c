#include <math.h>
#include <stdio.h>
#include "mex.h"
//#include "matrix.h"

#define PI 3.141592653589793

/*
 * qpoleHist2 - fast simulation of the powder histogram 
 * to use with nmrFit, qpole.m wrapper is needed!
 *
 * In second version also second order corrections are taken into account 
 * for satellites!
 *
 * Hamiltonian contains: Knight shift tensor, Quadrupole tensor
 *   tensor orientations of Knight shift and quadrupole int coincide!
 *
 * All units are in ppm!!!
 *
 * If Spin < 0 only central transition will be calculated
 *
 * from quadrupole.c
 *
 * This is a MEX-file for MATLAB.
 * Copyright 13.03.2010 - 06.03.2014   Anton Potocnik @ IJS F5  
 */
/* $Revision: 1.0 $ */



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


//********* First order Quadrupole Levels Splitting ***********
double quad_1stOrder(double Spin, double m, double niQ, double eta, double costh, double cos2fi)
{
    double mu2, smu2;
    mu2 = costh*costh;
    smu2 = 1-mu2; // sinth^2
    
    // Units: ppm
    return niQ/4.0*(m*m - Spin*(Spin+1)/3.0)*(3.0*mu2 - 1 + eta*smu2*cos2fi);
}

//********* Second order Quadrupole Levels Splitting ***********
double quad_2ndOrder(double Spin, double m, double niQ, double eta, double costh, double cos2fi, double niL)
{
    double a, mu2, smu2, B, T1, T2, T3;
    a = Spin*(Spin + 1);
    mu2 = costh*costh;
    B = cos2fi;
    smu2 = 1-mu2; // sinth^2
    
    T1 = -1./5.*(a - 3.*m*m)*(3.+eta*eta);
    T2 = 1./28.*(8.*a - 12.*m*m - 3.)*((eta*eta-3.)*(3.*mu2-1.) + 6.*eta*smu2*B);
    T3 = 1./8.*(18.*a - 34.*m*m - 5.)*(1./140*(18.+eta*eta)*(35.*mu2*mu2-30.*mu2+3.) + 3./7.*eta*smu2*(7.*mu2-1.)*B + 1./4.*eta*eta*smu2*smu2*(2.*B*B-1.));
    
    return 1./18*niQ*niQ/niL*m*(T1 + T2 + T3);   
}


void nmrpowder(double *f, double *Ws, mwSize N, mwSize g2, double Spin, double Kiso, double Kani, double Kasy, double niQ, double eta, double DniQ)
{
	int i,j,k,comp;
	double fi, costh, sinth, cos2fi, maxi=-100000, UNI, c[3], niL, ra;
    double U1, U2, Em1, Em2, fc, niQQ;
    double *Im, *W;
    
    niL = 1e6; //ppm
    
    

    // Calculate Quadrupole energy levels
	if (Spin < 0) {
        comp = 2;  
        Im = mxMalloc((comp+1)*sizeof(double*));
        W = mxMalloc((comp+1)*sizeof(double*));
        W[0] = 1; W[1] = 1;
        Im[0] = 0.5; Im[1] = -0.5;
        //printf("Only central line is calculated.\n");
    } else {
        comp = (int)(2*Spin+1);  
        Im = mxMalloc((comp+1)*sizeof(double*));
        W = mxMalloc((comp+1)*sizeof(double*));
        
        for (i=0; i<comp; i++) {
            if (i==0) Im[0] = Spin;
            else Im[i] = Im[i-1]-1;
            W[i] = Spin*(Spin+1) - Im[i]*(Im[i]-1);
            printf("Im = %f\n", Im[i]);
        }
    }
   
    
    // Initialize powder distribution
    //c[0]=1.;c[1]=2.;c[2]=1.; // sphereType = “full”
    //c[0]=-1.;c[1]=1.;c[2]=1.; // sphereType = “hemi”
    c[0]=-1.;c[1]=1.;c[2]=4.; // sphereType = “oct”
    // N = getNumberZCW(M); //total number of angles
    // g2 = getNumberZCW(M-2);
    
    //printf("N(ZCW) = %d\n", N);
    
    // Initialize niQ distribution
    srand(123456);
    
    j = 0; // count powder averages
    // prepare for niQ distribution
    U1 = (rand() % 1000)/1000.0; // [0 1)  uniform distribution
    U2 = (rand() % 1000)/1000.0; // [0 1)  uniform distribution
    for (k=0; k<comp-1; k++) {
        for (i=0; i<N; i++) {

            // Powder distribution
            UNI = (double)i/N;
            fi = 2*PI*mod(UNI*g2,1)/c[2];
            costh = c[0]*(c[1]*mod(UNI,1)-1);
            cos2fi = cos(2*fi);
            sinth = sin(acos(costh));

            // niQ distribution
            if (DniQ == 0) {
                niQQ = niQ;
            } else {
                // uniform distribution
                //ra = (rand() % 1000)/1000.0; // [0 1)  
                //niQQ = niQ + (2*ra-1)*DniQ;

                // Box–Muller ... normal distribution
                if (U1 < 0.00001) U1 = 0.00001;
                ra = sqrt(-2*log(U1))*cos(2*PI*U2);
                niQQ = niQ + ra*DniQ;
                U1 = U2;
                U2 = (rand() % 1000)/1000.0; // [0 1)  uniform distribution
            }

            
            Em1 = quad_1stOrder(fabs(Spin), Im[k], niQQ, eta, costh, cos2fi);
            Em2 = quad_1stOrder(fabs(Spin), Im[k+1], niQQ, eta, costh, cos2fi);
            
            if (Spin > 0 || (Im[k] == 0.5 && Im[k+1] == -0.5)) { // centralka
                Em1 += quad_2ndOrder(fabs(Spin), Im[k], niQQ, eta, costh, cos2fi, niL);
                Em2 += quad_2ndOrder(fabs(Spin), Im[k+1], niQQ, eta, costh, cos2fi, niL);
            }
            //printf("Imk = %f Imk+1 = %f\n",fabs(Im[k]),fabs(Im[k+1]));
  
            fc = Em1 - Em2;
            
//             if (Im[k] == 0.5 && Im[k+1] == -0.5) { // centralka
//                 fc0  = (-27.0/8.               + 9./4.*eta*cos2fi - 3./8.*eta*eta*cos2fi*cos2fi)*costh*costh*costh*costh;
//                 fc0 += ( 30.0/8. - 1./2.*eta*eta - 2.*eta*cos2fi   + 3./4.*eta*eta*cos2fi*cos2fi)*costh*costh;
//                 fc0 += ( -3.0/8. + 1./3.*eta*eta - 1./4.*eta*cos2fi - 3./8.*eta*eta*cos2fi*cos2fi);
//                 fc0 *= -1.0/6.*niQQ*niQQ/niL*(fabs(Spin)*(fabs(Spin)+1) - 3.0/4.);
//                 if (fabs(fc-fc0) > 1e-11) printf("dfc = %e\n",fc-fc0);
//             }
            
            fc += Kiso + Kani*(3.0*costh*costh-1)/2 + Kasy*sinth*sinth*cos2fi;
            
            f[k,j] = fc;
            j++;
        }
        
        Ws[k] = W[k];
   }
    
    mxFree(Im);
    //printf("Im has been freed!\n");
    mxFree(W);
    //printf("W has been freed!\n");
    
}



/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *f, *W;
  double Spin, Kiso, Kani, Kasy, niQ, eta, DniQ;
  mwSize mrows, ncols, N, m, g2, comp;
  
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
  }
   * 
   * create a pointer to the input matrix X 
   * X = mxGetPr(prhs[0]);
   * mrows = mxGetM(prhs[0]);
   * ncols = mxGetN(prhs[0]); */
  
	
    /*  get the scalars */
    N = mxGetScalar(prhs[0]); 
    Spin = mxGetScalar(prhs[1]);
    Kiso = mxGetScalar(prhs[2]);
    Kani = mxGetScalar(prhs[3]);
    Kasy = mxGetScalar(prhs[4]);
    niQ = mxGetScalar(prhs[5]);
    eta = mxGetScalar(prhs[6]);
    DniQ = mxGetScalar(prhs[7]);
    
    // Determine number of Quadrupole energy levels
	if (Spin < 0) {
        comp = 2;  
    } else {
        comp = (int)(2*Spin+1);  
    }
    comp--; //number of resonances is one less then number of levels

    m = getNumberZCW(N); //total number of angles
    g2 = getNumberZCW(N-2);
    mrows = m;
    ncols = 1;

    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(mrows,comp, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(comp,ncols, mxREAL);
    /*  create a C pointer to a copy of the output matrix */
    f = mxGetPr(plhs[0]);
    W = mxGetPr(plhs[1]);
    /*  call the C subroutine*/ 
    
    nmrpowder(f, W, m, g2, Spin, Kiso, Kani, Kasy, niQ, eta, DniQ);
}
