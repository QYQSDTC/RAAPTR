//#define GSL_DLL
#include "maxphase.h"
#include "LLR_PSOav.h"
#include "mat.h"
#include "matrix.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>  // root finding functions
#include <gsl/gsl_sort.h>
#include <gsl/gsl_integration.h>
#include <stdio.h>
#include <stdlib.h>

double LLR_PSOav(gsl_vector *xVec, /*!< Standardized Particle Coordinates*/
               void  *inParamsPointer /*!< Fitness function parameter structure
				                           containing information for conversion
				                          of standardized to real coordinates*/
			   ){

	//Set to 0 if angular variables do not have a periodic boundary conditions
	size_t wrapAngles = 1;

	size_t validPt;
   /* Cast from void * to known structure pointer before any
	   of the fields can be accessed.
	 */
	struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;

	double fitFuncVal;

	validPt = llrpsochkcoord(wrapAngles, inParams->rmin, inParams->rangeVec, xVec, inParams->realCoord);


	if (validPt){
		inParams->fitEvalFlag = 1;
	    fitFuncVal = AvPhaseLLR(inParams);
    }
	else{
		inParams->fitEvalFlag = 0;
		fitFuncVal = GSL_POSINF;
	}
   return fitFuncVal;
}

size_t  llrpsochkcoord(const size_t wrapAngles,
                       const gsl_vector *rmin, const gsl_vector *rangeVec,
                       gsl_vector *xVec, gsl_vector *realCoord){

	//Default is to assume coordinates are valid
    size_t validPt = 1;

    /* Convert from standardized to real coordinates. Return real coordinates
	   through realCoord. */
	s2rvector(xVec, rmin, rangeVec, realCoord);

	/* if angle wrapping not needed, just check standardized coordinates only */
	if(!wrapAngles){
		validPt = chkstdsrchrng(xVec);
		return validPt;
	}

	//Indices of coordinates in realCoord
	/*
	   coord #0: range of alpha is [0, 2*pi]
	         #1: range of delta is [pi/2,-pi/2]
	   Coord #2: Frequency
	         #3: range of phi0 is [0,pi]
	   Coord #4: Amplitude
	         #5: range of iota is [0,pi]
	         #6: range of thetaN is [0,pi]
	*/
	size_t alphaIndx = 0;
	size_t deltaIndx = 1;
	size_t omgIndx = 2;
	size_t phi0Indx = 3;
	size_t ampIndx = 4;
	size_t iotaIndx = 5;
	size_t thetaNIndx = 6;

	/* First check non-angular parameters. Standardized coordinates are enough. */
	if (gsl_vector_get(xVec,omgIndx)<0 || gsl_vector_get(xVec,omgIndx)>1){
		validPt = 0;
		return validPt;
	}
	if (gsl_vector_get(xVec,ampIndx)<0 || gsl_vector_get(xVec,ampIndx)>1){
		validPt = 0;
		return validPt;
	}

	//Check alpha
	double twoPi = 2*M_PI;
	double mn_alpha = gsl_vector_get(rmin,alphaIndx);
	double mx_alpha = gsl_vector_get(rangeVec,alphaIndx) + mn_alpha;
	double rng_alpha = mx_alpha - mn_alpha;
	if(mn_alpha < 0 || mx_alpha > twoPi){
		printf("Warning: Check the limits on alpha\n");
	}
	double alpha = gsl_vector_get(realCoord, alphaIndx);
	if (alpha < mn_alpha || alpha > mx_alpha){
	    if (alpha < 0){
	    	alpha = fmod(alpha, -twoPi);
			alpha = twoPi + alpha;
	    }
		else if(alpha > twoPi){
			alpha = fmod(alpha, twoPi);
		}
		//Reset real coordinate values
		gsl_vector_set(realCoord,alphaIndx,alpha);
		//Reset standardized coordinate value
		gsl_vector_set(xVec,alphaIndx,(alpha - mn_alpha)/rng_alpha);
		//check with boundaries in case search reagion does not cover the whole sphere
		if (alpha < mn_alpha || alpha > mx_alpha){
			validPt = 0;
			return validPt;
		}
	}

	//check delta
	double mn_delta = gsl_vector_get(rmin, deltaIndx);
	double mx_delta = gsl_vector_get(rangeVec, deltaIndx) + mn_delta;
	if (mn_delta < -M_PI_2 || mx_delta > M_PI_2){
		printf("Warning: Check the limits on delta\n");
	}
	double delta = gsl_vector_get(realCoord, deltaIndx);
	if (delta < mn_delta || delta > mx_delta){
		double polTheta = M_PI_2 - delta;
		if (polTheta < 0){
			polTheta = fmod(polTheta, -twoPi);
			if (polTheta > -M_PI){
				alpha = alpha + M_PI;
				polTheta = -polTheta;
				alpha = fmod(alpha, twoPi);
				//reset alpha
				gsl_vector_set(realCoord,alphaIndx,alpha);
				//Reset standardized value of alpha
				gsl_vector_set(xVec,alphaIndx,(alpha - mn_alpha)/rng_alpha);
			}
			else if (polTheta <= -M_PI){
				polTheta = twoPi + polTheta;
			}
		}
		else if (polTheta > M_PI){
			polTheta = fmod(polTheta, twoPi);
			if (polTheta < M_PI) {
				//Do nothing
			}
			else if (polTheta >= M_PI){
				polTheta = twoPi - polTheta;
				alpha = alpha + M_PI;
				//Reset alpha
				gsl_vector_set(realCoord,alphaIndx,alpha);
				//Reset standardized value of alpha
				gsl_vector_set(xVec,alphaIndx,(alpha - mn_alpha)/rng_alpha);
			}
		}
		delta = M_PI_2 - polTheta;
		//Reset real coordinate
		gsl_vector_set(realCoord,deltaIndx,delta);
		//Reset standardized coordinate value
		gsl_vector_set(xVec,deltaIndx,(delta-mn_delta)/(mx_delta - mn_delta));
		//check with boundaries in case search reagion does not cover the whole sphere
		if (alpha < mn_alpha || alpha > mx_alpha || delta < mn_delta || delta > mx_delta){
			validPt = 0;
			return validPt;
		}
	}

	//Check phi0
	double mn_phi0 = gsl_vector_get(rmin, phi0Indx);
	double mx_phi0 = gsl_vector_get(rangeVec, phi0Indx) + mn_phi0;
	if (mn_phi0 < 0 || mx_phi0 > M_PI){
		printf("Warning: Check the limits for phi0\n");
	}
	double phi0 = gsl_vector_get(realCoord,phi0Indx);
	if (phi0 < mn_phi0 || phi0 > mx_phi0){
		phi0 = wraphalfcircangle(phi0);
		//Reset real coordinate
		gsl_vector_set(realCoord,phi0Indx,phi0);
		//Reset standardized coordinate
		gsl_vector_set(xVec,phi0Indx,(phi0-mn_phi0)/(mx_phi0 - mn_phi0));
		if (phi0 < mn_phi0 || phi0 > mx_phi0){
			validPt = 0;
			return validPt;
		}
	}

	//Check iota
	double mn_iota = gsl_vector_get(rmin, iotaIndx);
	double mx_iota = gsl_vector_get(rangeVec, iotaIndx) + mn_iota;
	if (mn_iota < 0 || mx_iota > M_PI){
		printf("Warning: Check the limits for iota\n");
		abort();
	}
	double iota = gsl_vector_get(realCoord,iotaIndx);
	if (iota < mn_iota || iota > mx_iota){
		iota = wraphalfcircangle(iota);
		//Reset real coordinate
		gsl_vector_set(realCoord,iotaIndx,iota);
		//Reset standardized coordinate
		gsl_vector_set(xVec,iotaIndx,(iota - mn_iota)/(mx_iota - mn_iota));
		if (iota < mn_iota || iota > mx_iota){
			validPt = 0;
			return validPt;
		}
	}

	//Check thetaN
	double mn_thetaN = gsl_vector_get(rmin, thetaNIndx);
	double mx_thetaN = gsl_vector_get(rangeVec, thetaNIndx) + mn_thetaN;
	if (mn_thetaN < 0 || mx_thetaN > M_PI){
		printf("Warning: Check the limits for thetaN\n");
	}
	double thetaN = gsl_vector_get(realCoord, thetaNIndx);
	if (thetaN < mn_thetaN || thetaN > mx_thetaN){
		thetaN = wraphalfcircangle(thetaN);
		//Reset real coordinate
		gsl_vector_set(realCoord,thetaNIndx, thetaN);
		//Reset standardized coordinate
		gsl_vector_set(xVec,thetaNIndx,(thetaN - mn_thetaN)/(mx_thetaN - mn_thetaN));
		if (thetaN < mn_thetaN || thetaN > mx_thetaN){
			validPt = 0;
			return validPt;
		}
	}

	//Wrap up: No violations were found, so validPt = 1.
	return validPt;
}

/*! Wrap angles into [0, pi].*/
double wraphalfcircangle(double phi){
	double twoPi = 2*M_PI;
    double twoPhi = 2*phi;
	if (twoPhi < 0){
		twoPhi = twoPi + fmod(twoPhi, -twoPi);
	}
	else if(twoPhi > twoPi){
		twoPhi = fmod(twoPhi, twoPi);
	}
	return twoPhi/2;
}


double AvPhaseLLR(struct fitFuncParams *inParams){

  struct llr_pso_params *splParams = (struct llr_pso_params *)inParams->splParams;

  struct cfunc_OUTPUT * output;
  output = (struct cfunc_OUTPUT *)malloc(1 * sizeof(struct cfunc_OUTPUT));
  output->c = (double *)malloc(4 * sizeof(double));
  output->v = (double *)malloc(9 * sizeof(double));

	struct lh_OUTPUT * lhoutput;
  lhoutput = (struct lh_OUTPUT *)malloc(1 * sizeof(struct lh_OUTPUT));
  lhoutput->phiI = (double *)malloc(1 * sizeof(double));
  lhoutput->lhI = (double *)malloc(1 * sizeof(double));

  double tmp;
  double res;
  unsigned int lpr, i, j, j1, jj;
//  size_t pp[6];   not used.
  const size_t kk = 1, stride = 1, nn = 6;  // gsl_sort_largest_index
  //double src[6] = {1.0, 3.0, 6.0, 4.0, 5.0, 2.0};  //
  // src = (double *)malloc(6 * sizeof(double));
  // *(src+0)=1.0;
  // *(src+1)=3.0;
  // *(src+2)=6.0;
  // *(src+3)=4.0;
  // *(src+4)=5.0;
  // *(src+5)=2.0;

  unsigned int Np = splParams->Np;  // number of pulsars in PTA
  unsigned int N = splParams->N;   // number of samples

  // transfer parameters from structure inParams
  //printf("MP5: Np = %d\n", Np);
  //printf("MP5: N = %d\n", N);
  double *yr;
  //yr = (double *)malloc(N * sizeof(double));
  double *sd;
  //sd = malloc(Np * sizeof(double));
  double *alphaP, *deltaP;
  //alphaP = (double *)malloc(Np * sizeof(double));
  //deltaP = (double *)malloc(Np * sizeof(double));
  sd = splParams->sd;
  alphaP = splParams->alphaP;
  deltaP = splParams->deltaP;
  yr = splParams->yr;

  double *Phi;
  Phi = (double *)malloc(N * sizeof(double));
  double theta;
  double alpha, delta, omega, phi0, Amp, iota, thetaN;
  double **phiItmp, **lh, *phiI;


  unsigned int nDim = inParams->nDim;   // dimension of fitness function
  //printf("MP5: nDim = %d\n", nDim);

  double b[5];  // quartic equation coefficients, closed form solution
  b[0]=1.0;
  b[1]=2.0;
  b[2]=1.2;
  b[3]=2.5;
  b[4]=0.8;

  //double z[8];  // 4 (real) + 4 (complex) solutions from gsl_poly_complex_solve

  unsigned int nr;  // number of EFFECTIVE roots (real number && abs(r)<1)

  double LLR = 0.0;  // log likelihood ratio
  double **s;
  double C = 0.0;

  gsl_vector * skyLocSrc = gsl_vector_calloc(3); // sky location of source in Cartesian coordinate
  gsl_vector * skyLocPulsar = gsl_vector_calloc(3); // sky location of pulsars in Cartesian coordinate
  //gsl_vector * realCoord = gsl_vector_calloc(nDim);
  //gsl_vector_memcpy(realCoord,inParams->realCoord);

  alpha = gsl_vector_get(inParams->realCoord,0);
  delta = gsl_vector_get(inParams->realCoord,1);
  omega = gsl_vector_get(inParams->realCoord,2);
  phi0 = gsl_vector_get(inParams->realCoord,3);
  Amp = gsl_vector_get(inParams->realCoord,4);
	Amp = pow(10,Amp);  // physical amplitude
  gsl_vector_set(inParams->realCoord,4,Amp);
  iota = gsl_vector_get(inParams->realCoord,5);
  thetaN = gsl_vector_get(inParams->realCoord,6);

  s = (double **)malloc(Np * sizeof(double));
  for (i = 0; i < Np; i++) {
    //*(s+i) = (double *)malloc(N * sizeof(double));  // not needed!
    s[i] = splParams->s[i];
  }

  phiI = malloc(Np * sizeof(double));
// printf("MP5: s[0][1] = %f\n", *(*(s+0)+1));
// printf("MP5: s[1][2] = %f\n", *(*(s+1)+5));

  double *norm, *norm1,*LRn, M, norm2;
  norm = malloc(4 * sizeof(double));
  norm1 = malloc(4 * sizeof(double));
  LRn = malloc(4 * sizeof(double));
  double sign [4][5] = {{1.0,1.0,1.0,1.0,1.0},{1.0, -1.0, 1.0, -1.0, 1.0},
  {1.0, -1.0, -1.0, 1.0, 1.0},{1.0, 1.0, -1.0, -1.0, 1.0}};
  double intup[4] = {M_PI_2,M_PI,3 * M_PI_2,2 * M_PI};
  double intlow[4] = {0,M_PI_2,M_PI,3 * M_PI_2};


  gsl_vector_set(skyLocSrc, 0, cos(delta)*cos(alpha));
  gsl_vector_set(skyLocSrc, 1, cos(delta)*sin(alpha));
  gsl_vector_set(skyLocSrc, 2, sin(delta));


  for (i = 0; i < N; i++) {
    Phi[i] = yr[i] * omega;
    //printf("MP5: *(Phi+i) = %e, *(yr+i) = %e\n", *(Phi+i), *(yr+i));
  }


  double bs[5], tmp0, NN, result, error ;
  //double alpha0 = 2.0;

  gsl_integration_workspace *w
     = gsl_integration_workspace_alloc(1000);

double bx[2];
bx[0]= 1.0;
bx[1]= 2.0;

      gsl_function F;
      F.function = &f;
      F.params = b;
      gsl_integration_qags(&F,0,1.57,0,1e-5,1000,w,&result,&error);
      printf("LLR_PSOav: result = %.18f\n",result );

 //  for (i = 0; i < Np; i++){
 //  //printf("MP5: i = %d\n", i);
 //  //printf("alphaP = %f, deltaP = %f \n", *(alphaP+i), *(deltaP+i));
 //    gsl_vector_set(skyLocPulsar, 0, cos(deltaP[i])*cos(alphaP[i]) );
 //    gsl_vector_set(skyLocPulsar, 1, cos(deltaP[i])*sin(alphaP[i]) );
 //    gsl_vector_set(skyLocPulsar, 2, sin(deltaP[i]) );
 //
 //    gsl_blas_ddot(skyLocSrc, skyLocPulsar, &res);
 //  //printf("MP5: kp*k = %f\n", res);
 //
 //    theta = acos(res);
 //  //printf("theta = %f\n", theta);
 //
 //    cfunc(N, alpha, delta, alphaP[i], deltaP[i], theta,
 //          Amp, omega, iota, thetaN, phi0, Phi, s[i], (sd+i), output);
 //    tmp = (*output).v[2] - 0.5 * ((*output).v[6]+(*output).v[8]);
 //      b[0] = tmp;
 //    tmp = (*output).v[0] - (*output).v[5];
 //      b[1] = tmp;
 //    tmp = (*output).v[1] - (*output).v[7];
 //      b[2] = tmp;
 //    tmp = -0.5 * (*output).v[4];
 //      b[3] = tmp;
 //    tmp = 0.5 * ((*output).v[6]-(*output).v[3]);
 //      b[4] = tmp;
 //
 //    for (j = 0; j < 4; j++) {
 //      for ( j1 = 0; j1 < 4; j1++) {
 //        bs[j1] = b[j1] * sign[j][j1];
 //      }
 //      tmp0 = 0;
 //      for ( jj = 1; jj < 5; jj++) {
 //        if (bs(jj) > 0) {
 //          tmp0 = tmp0 + bs(jj);
 //        }
 //      }
 //      norm(j) = tmp0;
 //  gsl_function F;
 //  F.function = &f;
 //  F.params = &alpha0;
 //  gsl_integration_qags(&F,intlow(j),intup(j),0,1e-5,1000,w,&result,&error);
 //  printf("result = %.18f\n",result );
 //  LRn(j) = result;
 //  if (isinf LRn(j)) {
 //    printf("AvPhaseLLR Inf for PSR i = %d\n",i );
 //    printf("LRn at j = %d\n",j );
 //  else{ isnan (LRn(j))
 //    printf("AvPhaseLLR NAN for PSR i = %d\n",i );
 //    printf("LRn at j = %d\n",j );
 //      }
 //    }
 //  norm1(j) = norm(j) + log(LRn(j) + b[1]);
 //  }
 //    NN = max(norm1);
 //    norm2 = norm1 - NN;
 //    M = sum(exp(norm2));
 //    LLR = LLR + NN + log(M);
 // }
 gsl_integration_workspace_free(w);
 free(output->c);
 free(output->v);
 free(output);
 free(lhoutput->phiI);
 free(lhoutput->lhI);
 free(lhoutput);
 free(Phi);
 free(s);
 for (i = 0; i < Np; i++) {
   free(phiItmp[i]);
   free(lh[i]);
 }
 free(phiItmp);
 free(lh);
 free(phiI);
 gsl_vector_free(skyLocSrc);
 gsl_vector_free(skyLocPulsar);
   LLR = -LLR;
   return LLR;
}


double f (double x, void * params) {
  //double alpha0 = *(double *) params[0];
  double *b = (double *)params;
  printf("f: b = %f\n", b[1];
  //printf("f: b = %f\n", *(b+2);
  //b[1]=bb[1];
  //b[2]=bb[2];
  //b[3]=bb[3];
  //b[4]=bb[4];
  double f = exp(2.0) * cos(x); // + b[2] * sin(x) + b[3] * sin(2.0 * x) +
              //  b[4] * pow(cos(x),2.0); // - norm;
  return f;
}


// double f (double x, void * params) {
//   double alpha0 = *(double *) params;
//   //double b = *(double *)bb;
//   // printf("b = %f\n", b);
//   //b[1]=bb[1];
//   //b[2]=bb[2];
//   //b[3]=bb[3];
//   //b[4]=bb[4];
//   double f = exp(alpha0) * cos(x); // + b[2] * sin(x) + b[3] * sin(2.0 * x) +
//               //  b[4] * pow(cos(x),2.0); // - norm;
//   return f;
// }


void cfunc(unsigned int N, double alpha, double delta, double alphaP, double deltaP, double theta,
           double Amp, double omega, double iota, double thetaN, double phi0, double * Phi,
           double * s, double * sd, struct cfunc_OUTPUT * varargout)
{
  unsigned int i;
  double alphatilde;
  double a, b, c, d, e, f;
  double Pp, Pc, Fp, Fc;
  double A, psi;
  double *x, *y, *z;
  double sx, sy, sz, xx, xy, xz, yy, yz, zz;

  //double res = 9.876;
  //(*varargout).fitVal = res;

  //printf("cfunc: alpha = %f, alphaP = %f, Phi[5] = %f, sd = %f\n", alpha, alphaP, *(Phi+4), *sd);

  //for (i = 0; i < 6; i++) {
  //  printf("cfunc: s[%d] = %f\n", i, *(s+i));
  //}

  x = (double *)malloc(N * sizeof(double));
  y = (double *)malloc(N * sizeof(double));
  z = (double *)malloc(N * sizeof(double));

  // varargout->c = (double *)malloc(4 * sizeof(double));
  // varargout->v = (double *)malloc(9 * sizeof(double));

  //printf("cfunc: varargout.fitVal = %f\n", (*varargout).fitVal);
  alphatilde = alpha - alphaP;
  a = cos(deltaP);
  b = sin(deltaP);
  c = cos(alphatilde);
  d = sin(alphatilde);
  e = cos(delta);
  f = sin(delta);

  Pp = -pow(a,2.0) * ( 1.0 - 2.0 * pow(c,2.0) + pow(c,2.0) * pow(e,2.0) ) +
        pow(b,2.0) * pow(e,2.0) - 0.5 * sin(2.0 * deltaP) * c * sin(2.0 * delta);

  Pc = 2.0 * a * d * (a * c * f - b * e);
//printf("cfunc: Pp = %f\n", Pp);
  Fp = Pp / (1.0 - cos(theta));
  Fc = Pc / (1.0 - cos(theta));
//printf("cfunc: Fp = %f\n", Fp);
//printf("cfunc: Amp = %f\n", Amp);
  A = 2.0*Amp*sqrt( pow(1.0+pow(cos(iota),2.0),2.0)*pow(Fp*cos(2.0*thetaN)-Fc*sin(2.0*thetaN),2.0)
    + 4.0*pow(cos(iota),2.0)*pow((Fp*sin(2.0*thetaN)+Fc*cos(2.0*thetaN)),2.0) );
//printf("cfunc: A = %f\n", A);
   // tmp=-2*cos(iota)/(1+cos(iota)^2)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))/(Fp*cos(2*thetaN)-Fc*sin(2*thetaN));
   // solve psi, atan or atan2 ?
   // psi=atan(tmp);

  // psi = rt_atan2d_snf( -2.0 * cos(iota) * (Fp * sin(2.0 * thetaN) + Fc * cos(2.0 * thetaN)),
  //       (1.0 + pow(cos(iota),2.0)) * (Fp * cos(2.0 * thetaN) - Fc * sin(2.0 * thetaN)) );
  //
  psi = atan2( -2.0 * cos(iota) * (Fp * sin(2.0 * thetaN) + Fc * cos(2.0 * thetaN)),
         (1.0 + pow(cos(iota),2.0)) * (Fp * cos(2.0 * thetaN) - Fc * sin(2.0 * thetaN)) );
//("cfunc: psi = %f\n", psi);

  for (i = 0; i < N; i++) {
     x[i] = 0.5 * A * cos( psi + Phi[i] );
     y[i] = -0.5 * A * sin( psi + Phi[i] );
     z[i] = -0.5 * A * cos( 2*phi0 + psi + Phi[i] );
  }
//printf("cfunc: x[0] = %f, x[1] = %f\n", *(x+0), *(x+1));

  // c are combination of inner weighted product of s with X,Y,Z
//printf("cfunc: N= %d, s[5] = %f\n", N, *(s+5) );
  sx = InnProduct(N, s, x, *sd);  //  scalar
  sy = InnProduct(N, s, y, *sd);
  sz = InnProduct(N, s, z, *sd);
  xx = InnProduct(N, x, x, *sd);
  xy = InnProduct(N, x, y, *sd);
  xz = InnProduct(N, x, z, *sd);
  yy = InnProduct(N, y, y, *sd);
  yz = InnProduct(N, y, z, *sd);
  zz = InnProduct(N, z, z, *sd);
//printf("cfunc: sx= %f, sy = %f, zz =%f\n", sx,sy,zz);

  (*varargout).c[0] = -sx + xz;
  (*varargout).c[1] = sy - yz;
  (*varargout).c[2] = 0.5 * (xx - yy);
  (*varargout).c[3] = -xy;

  (*varargout).v[0] = sx;
  (*varargout).v[1] = sy;
  (*varargout).v[2] = sz;
  (*varargout).v[3] = xx;
  (*varargout).v[4] = xy;
  (*varargout).v[5] = xz;
  (*varargout).v[6] = yy;
  (*varargout).v[7] = yz;
  (*varargout).v[8] = zz;

  free(x);
  free(y);
  free(z);

}
// EOF cfunc

/*! \brief Inner product function for MaxPhase codes. */
double InnProduct(unsigned int N, double * X, double * Y, double sd)
{
  unsigned int i;
  double result;
  double c = 0.0;

  //printf("InnProduct: sd = %f\n", sd);
  //printf("InnProduct: X[0] = %f\n", *(X+0));

  for (i = 0; i < N; i++) {
    c += ( X[i] * Y[i] );
  }
//printf("InnProduct: c = %f\n", c);

  result = c / (sd * sd);
//printf("InnProduct: result = %f\n", result);

  return result;
}
// EOF: InnProduct

/*! Load data from .mat file into the special parameter structure
for LLR_PSO fitness function. Once used, the output should be
destroyed using llrparam_free. */
struct  llr_pso_params * loadfile2llrparam(MATFile *inputFilePr){
	/* Variables to hold input file content */
	mxArray *simParams, *timingResiduals, *yr;
	mxArray *field_Np, *field_N, *field_sd, *field_alphaP, *field_deltaP;
	double **s;

	yr = matGetVariable(inputFilePr,"yr");
	double *yr_Pr = mxGetPr(yr);
	timingResiduals = matGetVariable(inputFilePr,"timingResiduals");
	double *trPr = mxGetPr(timingResiduals);
	/* Extract fields of simParams */
	simParams = matGetVariable(inputFilePr,"simParams");
	field_Np = mxGetField(simParams, 0, "Np");
	field_N = mxGetField(simParams, 0, "N");
	field_sd = mxGetField(simParams, 0, "sd");
	field_alphaP = mxGetField(simParams, 0, "alphaP");
	field_deltaP = mxGetField(simParams, 0, "deltaP");
	double *sd_Pr = mxGetPr(field_sd);
	double *alphaP_Pr = mxGetPr(field_alphaP);
    double *deltaP_Pr = mxGetPr(field_deltaP);


	/* Load fitness function parameter structure */
	size_t lpc1, lpc2;
	size_t Np, N;
	struct llr_pso_params *llp = (struct llr_pso_params *)malloc(sizeof(struct llr_pso_params));
	Np = (size_t)mxGetScalar(field_Np);
	N = (size_t)mxGetScalar(field_N);
	llp->Np = (unsigned int)Np;
	llp->N = (unsigned int)N;
	llp->sd = (double *)malloc(Np*sizeof(double));
	llp->alphaP = (double *)malloc(Np*sizeof(double));
	llp->deltaP = (double *)malloc(Np*sizeof(double));
	llp->phiI = (double *)malloc(Np*sizeof(double));
	for (lpc1 = 0; lpc1 < Np; lpc1++){
		llp->sd[lpc1] = sd_Pr[lpc1];
		llp->alphaP[lpc1] = alphaP_Pr[lpc1];
		llp->deltaP[lpc1] = deltaP_Pr[lpc1];
	}
	llp->yr = (double *) malloc(N*sizeof(double));
	for (lpc1 = 0; lpc1 < N; lpc1++){
		llp->yr[lpc1] = yr_Pr[lpc1];
	}

	/* load timing residuals into fitness function param struct */
	/* Allocate storage for timing residuals field 's' in fitness
	function parameter struct. This is needed because
	timingResiduals is a 2D array while it is actually stored as
	a 1 D array in the mxArray.
	*/
	s = (double **)malloc(Np*sizeof(double *));
	unsigned int countElements = 0;
	for (lpc1 = 0; lpc1 < Np; lpc1++){
		s[lpc1] = (double *) malloc(N*sizeof(double));
	}
	for (lpc1 = 0; lpc1 < N; lpc1++){
		for (lpc2 = 0; lpc2 < Np; lpc2++){
			s[lpc2][lpc1] = trPr[countElements];
			countElements++;
		}
	}
	llp->s = s;

	//Wrap up
	mxDestroyArray(timingResiduals);
	mxDestroyArray(yr);
	/* Destroying an mxArray holding a structure also
	   destroys its fields. So no need to deallocate
	   Np, N, sd, deltaP, alphaP separately.
	*/
	mxDestroyArray(simParams);

	return llp;
}

/*! Deallocate special parameter structure specific to LLR_PSO fitness function.
*/
void llrparam_free(struct llr_pso_params *llp){

	size_t lpc;
	size_t Np = llp->Np;

	for (lpc = 0; lpc < Np; lpc++){
		free(llp->s[lpc]);
	}
	free(llp->s);
	free(llp->sd);
	free(llp->alphaP);
	free(llp->deltaP);
	free(llp->yr);
	free(llp->phiI);
	free(llp);
}
