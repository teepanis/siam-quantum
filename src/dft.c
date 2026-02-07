/*
	Siam Quantum [SQ] is a program package that performs quantum
	mechanical calculations on atomic systems.

	Copyright (C) 2015 Teepanis Chachiyo <teepanisc@nu.ac.th>,
	The Institute for Fundamental Study (IF),
	Naresuan University, Thailand.

	This file is a part of SQ.
	                                                       
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2, or (at your option) 
	any later version.                                                  
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "option.h"
#include "dft.h"
#include "uhf.h"

#ifdef LIBXC
#include <xc.h>
#endif

//
// cleanBGrid: clean memory allocated by the structure
//
// May 9, 2017 - Teepanis Chachiyo
//   Add gradient
//
// May 11, 2016 - Teepanis Chachiyo
//   Initial implementation and testing
//
struct BGrid_t * cleanBGrid(struct BGrid_t *bg){
	if(bg->maxr2 != NULL) free(bg->maxr2);
	if(bg->val   != NULL) free(bg->val);
	if(bg->gx    != NULL) free(bg->gx);
	if(bg->gy    != NULL) free(bg->gy);
	if(bg->gz    != NULL) free(bg->gz);
	if(bg->x     != NULL) free(bg->x);
	if(bg->y     != NULL) free(bg->y);
	if(bg->z     != NULL) free(bg->z);

	if(bg != NULL) free(bg);
	return NULL;
}

//
// genBasisGrid: evaluate basis function value at the grid points.
// It returns 1D array whose element is bGrid[n*nBasis+i]
//
// May 9, 2017 - Teepanis Chachiyo
//    Add gradients
//
// May 8, 2016 - Teepanis Chachiyo
//    Truncate grid using cutoff 
//
// May 5, 2016 - Teepanis Chachiyo
//    Initial implementation and testing
//
struct BGrid_t * genBasisGrid(
const int nBasis,               // number of basis function
const struct GTOBasis_t *gto,   // basis function data
const struct MolGrid_t *grid,   // grid point data
const double cutoff){           // cutoff value

	int i;                      // basis function index
	int c;                      // contraction index
	int n;                      // grid point index
	double ir2;                 // radial distance squared
	struct BGrid_t  *bg;        // bGrid data
	double *ptr;                // pointer to basis grid data
	double *ptrx, *ptry, *ptrz; // pointer to basis gradient data

	// allocation
	bg = calloc(1, sizeof(struct BGrid_t));
	if(bg==NULL){
		printf("genBasisGrid - error cannot allocate memory\n");
		exit(-1);
	}
	bg->nBasis = nBasis;
	bg->maxr2  = calloc(nBasis, sizeof(double));
	bg->x      = calloc(nBasis, sizeof(double));
	bg->y      = calloc(nBasis, sizeof(double));
	bg->z      = calloc(nBasis, sizeof(double));
	if(bg->maxr2==NULL || bg->x==NULL || bg->y==NULL || bg->z==NULL){
		printf("genBasisGrid - error cannot allocate memory\n");
		exit(-1);
	}

	// load basis centers
	for(i=0; i < nBasis; i++){
		bg->x[i] = gto[i].x0;
		bg->y[i] = gto[i].y0;
		bg->z[i] = gto[i].z0;
	}

	// compute maximum r2 for basis functions
	for(i=0; i < nBasis; i++){

		// get the minimum exponent
		double minexp = gto[i].exp[0];
		for(c=0; c < gto[i].nContract; c++)
			if(minexp > gto[i].exp[c]) minexp = gto[i].exp[c];

		// compute maximum r2 
		bg->maxr2[i] = -log(cutoff)/minexp;
	}

	// count the number of element
	bg->nPoint = 0;
	for(n=0; n < grid->nPoint; n++)
	for(i=0; i < nBasis; i++){

		ir2 = (grid->x[n] - bg->x[i]) * (grid->x[n] - bg->x[i]) +
		      (grid->y[n] - bg->y[i]) * (grid->y[n] - bg->y[i]) +
		      (grid->z[n] - bg->z[i]) * (grid->z[n] - bg->z[i]);
		if(ir2 > bg->maxr2[i]) continue;
		bg->nPoint++;
	}

	// allocate value storage
	bg->val = calloc(bg->nPoint, sizeof(double));
	bg->gx  = calloc(bg->nPoint, sizeof(double));
	bg->gy  = calloc(bg->nPoint, sizeof(double));
	bg->gz  = calloc(bg->nPoint, sizeof(double));
	if(bg->val==NULL || bg->gx==NULL || bg->gy==NULL || bg->gz==NULL){
		printf("genBasisGrid - error cannot allocate memory\n");
		exit(-1);
	}

	// main loop
	double rad;              // radial part
	double ang;              // angular part
	ptr  = bg->val;

	double radg;             // radial part of gradient
	double angx, angy, angz; // angular part in gradient
	ptrx = bg->gx;
	ptry = bg->gy;
	ptrz = bg->gz;
	for(n=0; n < grid->nPoint; n++)
	for(i=0; i < nBasis; i++){

		ir2 = (grid->x[n] - bg->x[i]) * (grid->x[n] - bg->x[i]) +
		      (grid->y[n] - bg->y[i]) * (grid->y[n] - bg->y[i]) +
		      (grid->z[n] - bg->z[i]) * (grid->z[n] - bg->z[i]);
		if(ir2 > bg->maxr2[i]) continue;

		// evaluate basis function at the point
		ang = pow_int(grid->x[n]-gto[i].x0, gto[i].l)
		    * pow_int(grid->y[n]-gto[i].y0, gto[i].m)
		    * pow_int(grid->z[n]-gto[i].z0, gto[i].n);
		for(rad=0.0,c=0; c < gto[i].nContract; c++)
			rad += gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );

		// evaluate gradient at the point
		angx = (gto[i].l<1)?0: gto[i].l * pow_int(grid->x[n]-gto[i].x0, gto[i].l-1)
		                                * pow_int(grid->y[n]-gto[i].y0, gto[i].m  )
		                                * pow_int(grid->z[n]-gto[i].z0, gto[i].n  );

		angy = (gto[i].m<1)?0: gto[i].m * pow_int(grid->x[n]-gto[i].x0, gto[i].l  )
		                                * pow_int(grid->y[n]-gto[i].y0, gto[i].m-1)
		                                * pow_int(grid->z[n]-gto[i].z0, gto[i].n  );

		angz = (gto[i].n<1)?0: gto[i].n * pow_int(grid->x[n]-gto[i].x0, gto[i].l  )
		                                * pow_int(grid->y[n]-gto[i].y0, gto[i].m  )
		                                * pow_int(grid->z[n]-gto[i].z0, gto[i].n-1);

		for(radg=0.0,c=0; c < gto[i].nContract; c++)
			radg -= 2.0 * gto[i].exp[c] * gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );

		// store value
		*ptr  = ang  * rad; ptr++;
		*ptrx = angx * rad + ang * (grid->x[n]-gto[i].x0) * radg; ptrx++;
		*ptry = angy * rad + ang * (grid->y[n]-gto[i].y0) * radg; ptry++;
		*ptrz = angz * rad + ang * (grid->z[n]-gto[i].z0) * radg; ptrz++;

	}

	return bg;
}


//
// getRho: computes electron density at the
// grid points using previously evaluated basis function grid.
//
// May 6, 2016 - Teepanis Chachiyo
//    Compute alpha and beta separately
//
// May 5, 2016 - Teepanis Chachiyo
//    Initial implementation and testing
//
void getRho(
const int nBasis,              // number of basis function
const double *PA,              // spin up density matrix
const double *PB,              // spin dn density matrix
const struct GTOBasis_t *gto,  // basis function info
const struct MolGrid_t *grid,  // grid point data
const double cutoff,           // basis function cutoff
      double *rhoa,            // returned spin up density 
      double *rhob){           // returned spin dn density

	int i,j;       // basis function index
	int n;         // grid point index
	int c;         // basis contraction index
	double *val;   // basis function values
	double *maxr2; // maximum radius

	// allocate basis function buffer and max radius
	val   = calloc(nBasis, sizeof(double));
	maxr2 = calloc(nBasis, sizeof(double));
	if(val==NULL || maxr2==NULL){
		printf("getRho - error cannot allocate memory\n");
		exit(-1);
	}

	// compute maximum r2 for basis functions
	for(i=0; i < nBasis; i++){

		// get the minimum exponent
		double minexp = gto[i].exp[0];
		for(c=0; c < gto[i].nContract; c++)
			if(minexp > gto[i].exp[c]) minexp = gto[i].exp[c];

		// compute maximum r2 
		maxr2[i] = -log(cutoff)/minexp;
	}

	// set zero
	for(n=0; n < grid->nPoint; n++){ 
		rhoa[n] = 0.0;
		rhob[n] = 0.0;
	}

	// main loop
	double ival, ir2;
	for(n=0; n < grid->nPoint; n++){

		// evaluate all basis function at this point
		for(i=0; i < nBasis; i++){

			ir2 = (grid->x[n] - gto[i].x0) * (grid->x[n] - gto[i].x0) +
			      (grid->y[n] - gto[i].y0) * (grid->y[n] - gto[i].y0) +
			      (grid->z[n] - gto[i].z0) * (grid->z[n] - gto[i].z0);
			if(ir2 > maxr2[i]){
				val[i] = 0.0;
				continue;
			}

			for(ival=0.0,c=0; c < gto[i].nContract; c++)
				ival += gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );
			ival = ival * pow_int(grid->x[n]-gto[i].x0, gto[i].l)
			            * pow_int(grid->y[n]-gto[i].y0, gto[i].m)
			            * pow_int(grid->z[n]-gto[i].z0, gto[i].n);

			val[i] = ival;
		}

		for(i=0; i < nBasis; i++){

			ival = val[i];
			if(ival==0.0) continue;

			for(j=0; j <= i; j++){
				// symmetric properties
				if(j==i){
					rhoa[n] +=     PA[i*nBasis+j]*ival*val[j];
					rhob[n] +=     PB[i*nBasis+j]*ival*val[j];
				}else{
					rhoa[n] += 2.0*PA[i*nBasis+j]*ival*val[j];
					rhob[n] += 2.0*PB[i*nBasis+j]*ival*val[j];
				}
			}
		}
	}

	// cut-off
	for(n=0; n < grid->nPoint; n++){
		if(rhoa[n] < -RHO_CUTOFF || rhob[n] < -RHO_CUTOFF){
			printf("getRho - error rho less than zero\n");
			exit(-1);
		}
		if(fabs(rhoa[n] < RHO_CUTOFF)) rhoa[n] = 0.0;
		if(fabs(rhob[n] < RHO_CUTOFF)) rhob[n] = 0.0;
	}

	// cleanup
	free(val);
	free(maxr2);
}


//
// getXC: computes exchange-correlation matrix,
// energy functional at grid points, and forces
//
// May 9, 2017 - Teepanis Chachiyo
//   Add gradient and The Chachiyo enhancement factor
//
// August 14, 2016 - Teepanis Chachiyo
//   Add Chachiyo correlation functional
//
// May 11, 2016 - Teepanis Chachiyo
//   Can now compute forces
//
// May 9, 2016 - Teepanis Chachiyo
//   Initial implementation and testing
// 
void getXC(
const int nBasis,              // number of basis function
const double *PA,              // spin up density matrix
const double *PB,              // spin dn density matrix
const struct MolGrid_t *grid,  // grid point data
const struct BGrid_t *bg,      // basis function grid
const struct option_t *opt,    // global options
const struct GTOBasis_t *gto,  // basis function info
const struct Molecule_t *mol,  // molecular info
const double cutoff,           // basis product cutoff
const int getMatrix,           // flag to compute XC matrix
const int getForce,            // flag to compute forces
      double *Exc,             // returned (incremental) energy grid
      double *XA,              // returned spin up matrix
      double *XB,              // returned spin dn matrix
      double *Fx,              // returned (incremental) force
      double *Fy,              // returned (incremental) force
      double *Fz){             // returned (incremental) force

	int i,j;                  // basis function index
	int c;                    // contraction index
	int n;                    // grid point index
	double ival;              // evaluated ith basis function value
	double ir2;               // radius from ith basis
	double *val;              // basis function values
	double *valx,*valy,*valz; // basis function gradient
	double dEdrhoa, dEdrhob;  // functional derivative
	double rhoa,  rhob;       // electron densities
	double rhoax, rhoay, rhoaz;  // gradient of spin-up electron densities
	double rhobx, rhoby, rhobz;  // gradient of spin-dn electron densities

	// forces on basis centers
	double *fx=NULL, *fy=NULL, *fz=NULL;
	if(getForce){
		fx = calloc(nBasis, sizeof(double));
		fy = calloc(nBasis, sizeof(double));
		fz = calloc(nBasis, sizeof(double));
		if(fx==NULL || fy==NULL || fz==NULL){
			printf("getXC - error cannot allocate memory\n");
			exit(-1);
		}
	}

	// allocation
	val  = calloc(nBasis, sizeof(double));
	valx = calloc(nBasis, sizeof(double));
	valy = calloc(nBasis, sizeof(double));
	valz = calloc(nBasis, sizeof(double));
	if(val==NULL || valx==NULL || valy==NULL || valz==NULL){
		printf("getXC - error cannot allocate memory\n");
		exit(-1);
	}

	// main loop
	register double suma, sumb;
	double maxval;
	double *ptr  = bg->val;
	double *ptrx = bg->gx;
	double *ptry = bg->gy;
	double *ptrz = bg->gz;
	for(n=0; n < grid->nPoint; n++){

		// load all basis function at this point
		maxval=0.0;
		for(i=0; i < nBasis; i++){

			ir2 = (grid->x[n] - bg->x[i]) * (grid->x[n] - bg->x[i]) +
			      (grid->y[n] - bg->y[i]) * (grid->y[n] - bg->y[i]) +
			      (grid->z[n] - bg->z[i]) * (grid->z[n] - bg->z[i]);
			if(ir2 > bg->maxr2[i]){ val [i] = 0.0;
			                        valx[i] = 0.0;
			                        valy[i] = 0.0;
			                        valz[i] = 0.0; continue; }

			val [i] = *ptr ;  ptr++;
			valx[i] = *ptrx; ptrx++;
			valy[i] = *ptry; ptry++;
			valz[i] = *ptrz; ptrz++;

			if(maxval < fabs(val[i])) maxval = fabs(val[i]);
		}

		// compute electron densities at this point
		rhoa   = 0.0; rhob   = 0.0;
		rhoax  = 0.0; rhoay  = 0.0; rhoaz = 0.0;
		rhobx  = 0.0; rhoby  = 0.0; rhobz = 0.0;

/*
		for(i=0; i < nBasis; i++){

			if(fabs(val[i])*maxval < cutoff) continue;

			suma  = 0.0; sumb  = 0.0;
			for(j=0; j < nBasis; j++){
				suma  += PA[i*nBasis+j]*val[j];
				sumb  += PB[i*nBasis+j]*val[j];
			}
			rhoa  += val[i]*suma; rhob += val[i]*sumb;
		}
*/

		for(i=0; i < nBasis; i++){
			//if(fabs(val[i])*maxval < cutoff) continue;
			// Above screening caused error for testing the Chachiyo 
			// exchange for B,C,N,O atoms in the 3rd decimal place.
			// Putting it back should give the result in the preprint
			// archive.
			// Teepanis Chachiyo 6/7/2017
			//
			suma = 0.0; sumb = 0.0;
			for(j=0; j < i; j++){
				suma  += PA[i*nBasis+j]*val[j];
				sumb  += PB[i*nBasis+j]*val[j];
			}
			rhoa += 2.0*val[i]*suma; rhob += 2.0*val[i]*sumb;
			rhoa += val[i]*PA[i*nBasis+i]*val[i];
			rhob += val[i]*PB[i*nBasis+i]*val[i];

			for(j=0; j < nBasis; j++){
				rhoax += PA[i*nBasis+j] * (valx[i]*val[j]+valx[j]*val[i]);
				rhoay += PA[i*nBasis+j] * (valy[i]*val[j]+valy[j]*val[i]);
				rhoaz += PA[i*nBasis+j] * (valz[i]*val[j]+valz[j]*val[i]);
				rhobx += PB[i*nBasis+j] * (valx[i]*val[j]+valx[j]*val[i]);
				rhoby += PB[i*nBasis+j] * (valy[i]*val[j]+valy[j]*val[i]);
				rhobz += PB[i*nBasis+j] * (valz[i]*val[j]+valz[j]*val[i]);
			}
		}

		// enforce positive electron densities
		rhoa = fabs(rhoa);
		rhob = fabs(rhob);

		// compute energy functional and potential
		dEdrhoa = 0.0; dEdrhob = 0.0;
		switch(opt->method){

		case METHOD_SVWN5:
			getDFT_xSlater(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			getDFT_cVWN5  (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_HALF:
			getDFT_xSlater(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			dEdrhoa = 0.5*dEdrhoa;
			dEdrhob = 0.5*dEdrhob;
			Exc[n]  = 0.5*Exc[n];
			getDFT_cVWN5  (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_HF2xDFT:

			if((rhoa+rhob) > RHO_CUTOFF){
				// compute Dirac-Slater contribution for uniform electron gas
				double exa = -3.0/4.0*pow(3.0/M_PI,1.0/3.0)*pow(2.0*rhoa,4.0/3.0);
				double exb = -3.0/4.0*pow(3.0/M_PI,1.0/3.0)*pow(2.0*rhob,4.0/3.0);

				double grhoa = 2.0*sqrt(rhoax*rhoax + rhoay*rhoay + rhoaz*rhoaz);
				double grhob = 2.0*sqrt(rhobx*rhobx + rhoby*rhoby + rhobz*rhobz);
				double sa = grhoa/pow(2.0*rhoa,4.0/3.0)/2.0/pow(3.0*M_PI*M_PI,1.0/3.0);
				double sb = grhob/pow(2.0*rhob,4.0/3.0)/2.0/pow(3.0*M_PI*M_PI,1.0/3.0);
				double Fa=1.0; // spin-up enhancement factor
				double Fb=1.0; // spin-dn enhancement factor

				switch(opt->whichHF2xDFT){
				case HF2xDFT_xSLATER:
#define F(s) 1.0
					Fa = F(sa);
					Fb = F(sb);
#undef F
				break;

				case HF2xDFT_xPY86:
#define F(s) pow(1.0+1.296*s*s+14.0*s*s*s*s+0.2*s*s*s*s*s*s,1.0/15.0)
					Fa = F(sa);
					Fb = F(sb);
#undef F
				break;

				case HF2xDFT_xBECKE88:
#define x(s) (s*2.0*pow(3.0*M_PI*M_PI,1.0/3.0))
#define Ax (3.0/4.0*pow(3.0/M_PI,1.0/3.0))
#define beta (0.0053)
#define F(s) 1.0 + beta/Ax*x(s)*x(s)/(1.0+6*x(s)*beta*asinh(x(s)*pow(2.0,1.0/3.0)))
					Fa = F(sa);
					Fb = F(sb);
#undef F
#undef x
#undef Ax
#undef beta
				break;

				case HF2xDFT_xPBE:
#define kappa (0.804)
#define mu (0.21951)
#define F(s) 1.0+kappa-kappa/(1.0+mu/kappa*s*s)
					Fa = F(sa);
					Fb = F(sb);
#undef F
#undef mu
#undef kappa
				break;

				case HF2xDFT_xCHACHIYO:
#define F(x) (3*x*x + M_PI*M_PI*log(x+1.0))/(3*x+M_PI*M_PI)/log(x+1.0)
					Fa = F(sa*4.0*M_PI/9.0);
					Fb = F(sb*4.0*M_PI/9.0);
#undef F

				break;
				}

				dEdrhoa = 0.0; // do not use DFT potential
				dEdrhob = 0.0; // do not use DFT potential
				if(rhoa > RHO_CUTOFF) Exc[n]  = exa*Fa/2.0;
				if(rhob > RHO_CUTOFF) Exc[n] += exb*Fb/2.0;

			}
		break;

		case METHOD_SCHACHIYO:
			getDFT_xSlater   (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			getDFT_cChachiyo (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		default:
			printf("getXC - error unknown DFT METHOD\n");
			exit(-1);
		break;
		}

		// coupled with grid weight
		dEdrhoa *= grid->w[n];
		dEdrhob *= grid->w[n];

		// evaluate matrix element
		if(getMatrix)
		for(i=0; i < nBasis; i++){

			// load basis function value
			ival = val[i];
			if(fabs(ival*maxval) < cutoff) continue;

			// quadrature sum
			for(j=0; j <= i; j++){
				XA[i*nBasis+j] += dEdrhoa*ival*val[j];
				XB[i*nBasis+j] += dEdrhob*ival*val[j];
			}
		}

		// evaluate molecular forces
		if(getForce){

			// compute gradient of basis function
			for(i=0; i < nBasis; i++){

				ir2 = (grid->x[n] - bg->x[i]) * (grid->x[n] - bg->x[i]) +
				      (grid->y[n] - bg->y[i]) * (grid->y[n] - bg->y[i]) +
				      (grid->z[n] - bg->z[i]) * (grid->z[n] - bg->z[i]);

				// screening
				if(ir2 > bg->maxr2[i]){ valx[i] = valy[i] = valz[i] = 0.0; continue; }

				suma=sumb=0.0;
				for(c=0; c < gto[i].nContract; c++){
					double radial = gto[i].norm[c]*gto[i].coef[c]*exp(-gto[i].exp[c]*ir2);
					suma += radial;
					sumb -= gto[i].exp[c]*radial;
				}
				double ax = pow_int(grid->x[n]-gto[i].x0, gto[i].l);
				double ay = pow_int(grid->y[n]-gto[i].y0, gto[i].m);
				double az = pow_int(grid->z[n]-gto[i].z0, gto[i].n);

				valx[i] = 2.*(grid->x[n] - gto[i].x0)*ax*ay*az*sumb;
				valy[i] = 2.*(grid->y[n] - gto[i].y0)*ax*ay*az*sumb;
				valz[i] = 2.*(grid->z[n] - gto[i].z0)*ax*ay*az*sumb;

				if(gto[i].l>0) valx[i] += gto[i].l * pow_int(grid->x[n]-gto[i].x0, gto[i].l-1)*ay*az*suma;
				if(gto[i].m>0) valy[i] += gto[i].m * pow_int(grid->y[n]-gto[i].y0, gto[i].m-1)*ax*az*suma;
				if(gto[i].n>0) valz[i] += gto[i].n * pow_int(grid->z[n]-gto[i].z0, gto[i].n-1)*ax*ay*suma;
			}

			// evaluate forces
			for(i=0; i < nBasis; i++){

				if(valx[i]==0.0 && valy[i]==0.0 && valz[i]==0.0) continue;

				suma  = 0.0; sumb  = 0.0;
				for(j=0; j < nBasis; j++){
					suma  += PA[i*nBasis+j]*val[j];
					sumb  += PB[i*nBasis+j]*val[j];
				}

				// accumulate forces
				fx[i] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * valx[i];
				fy[i] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * valy[i];
				fz[i] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * valz[i];
			}

		}

	}

	// make use of the symmetric property
	if(getMatrix)
	for(i=0; i < nBasis; i++)
	for(j=0; j < i; j++){
		XA[j*nBasis+i] = XA[i*nBasis+j];
		XB[j*nBasis+i] = XB[i*nBasis+j];
	}

	// translate forces to atomic center
	if(getForce){
		for(j=0; j < nBasis; j++){
			for(i=0; i < mol->nAtom; i++)
				if(gto[j].x0==mol->x[i] && 
				   gto[j].y0==mol->y[i] && 
				   gto[j].z0==mol->z[i])
					break;
			Fx[i] += fx[j];
			Fy[i] += fy[j];
			Fz[i] += fz[j];
		}
	}

	// cleanup
	free(val);
	free(valx);
	free(valy);
	free(valz);

	if(fx!=NULL) free(fx);
	if(fy!=NULL) free(fy);
	if(fz!=NULL) free(fz);
}

struct gridBasisCenter_t * cleanGridBasisCenter(
int nBasis,
struct gridBasisCenter_t *gbCenter){

	if(gbCenter->gCenter!= NULL) free(gbCenter->gCenter);
	if(gbCenter->bR2Max != NULL) free(gbCenter->bR2Max);
	if(gbCenter->bR2    != NULL) free(gbCenter->bR2);
	if(gbCenter->bAtomID!= NULL) free(gbCenter->bAtomID);
	if(gbCenter         != NULL) free(gbCenter);
	return NULL;
}

struct gridBasisCenter_t * genGridBasisCenter(
int nBasis,
const struct GTOBasis_t *gto,
const struct Molecule_t *mol,
const struct MolGrid_t *grid){

	struct gridBasisCenter_t *gB;
	int i,j,n,c;
	double ival,ir2;

	// allocate main structure
	gB = calloc(1, sizeof(struct gridBasisCenter_t));
	if(gB==NULL){
		printf("genGridBasisCenter - error cannot allocate gB\n");
		exit(-1);
	}

	gB->gCenter = calloc(grid->nPoint,sizeof(int));
	gB->bR2Max  = calloc(nBasis, sizeof(double));
	gB->bR2     = calloc(nBasis*nBasis, sizeof(double));
	gB->bAtomID = calloc(nBasis, sizeof(int));
	if(   gB->gCenter==NULL || gB->bR2Max ==NULL ||
	      gB->bR2    ==NULL || gB->bAtomID==NULL){
			printf("genGridBasisCenter - error cannot allocate screening\n");
			exit(-1);
	}
	// compute distance matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
		gB->bR2[i*nBasis+j] = DIST2(gto[i].x0, gto[i].y0, gto[i].z0, gto[j].x0, gto[j].y0, gto[j].z0);
	}
	// find nearest center
	double avgR2=0.0;
	for(n=0; n < grid->nPoint; n++){
		gB->gCenter[n] = 0;
		ival = 1.0E30;
		for(i=0; i < nBasis; i++){
			ir2 = DIST2(grid->x[n],grid->y[n],grid->z[n], gto[i].x0, gto[i].y0, gto[i].z0);
			if(ir2 < ival){
				gB->gCenter[n] = i;
				ival = ir2;
			}
		}
		avgR2 += ival;
	}
	avgR2 = avgR2/grid->nPoint;
#undef DIST2
	// compute basis maximum radius
	for(i=0; i < nBasis; i++){

		// get the minimum exponent
		double minexp = gto[i].exp[0];
		double weight = gto[i].norm[0] * gto[i].coef[0];
		for(c=0; c < gto[i].nContract; c++)
			if(minexp > gto[i].exp[c]){
				minexp = gto[i].exp[c];
				weight = gto[i].norm[c] * gto[i].coef[c];
			}

		// compute maximum r2 
		gB->bR2Max[i] = pow(sqrt(-log(BGRID_CUTOFF/weight)/minexp) + sqrt(avgR2), 2.0);
	}

	// compute basis atomID
	for(i=0; i < nBasis; i++){
		gB->bAtomID[i] = -1;  // undetermined
		for(j=0; j < mol->nAtom; j++)
			if(gto[i].x0 == mol->x[j] &&
			   gto[i].y0 == mol->y[j] &&
			   gto[i].z0 == mol->z[j]){
				gB->bAtomID[i] = j;
				break;
			}
		if(gB->bAtomID[i] == -1){
			printf("genGridBasisCenter : error cannot find matching atom center\n");
			exit(-1);
		}
	}

	return gB;
}

//
// getXCDirect: computes exchange-correlation matrix,
// energy functional at grid points, and forces
//
// July 6, 2017 - Teepanis Chachiyo
//   Modified from getXC to compute the basis function
//   at a point instead of storing them in memory.
//   This turned out to be faster and more accurate
//   than the getXC counter part. From tests I suspect
//   the scaling is also linear, but need systematic 
//   test to confirm this observation.
//
//////// Copied from getXC LOG /////////
//
// May 9, 2017 - Teepanis Chachiyo
//   Add gradient and The Chachiyo enhancement factor
//
// August 14, 2016 - Teepanis Chachiyo
//   Add Chachiyo correlation functional
//
// May 11, 2016 - Teepanis Chachiyo
//   Can now compute forces
//
// May 9, 2016 - Teepanis Chachiyo
//   Initial implementation and testing
// 
void getXCDirect(
const int nBasis,             // number of basis function
const int nEA,                // number of electron for spin-up
const int nEB,                // number of electron for spin-dn
const double *CA,             // molecular orbitals for spin-up
const double *CB,             // molecular orbitals for spin-dn
const double *PA,             // spin up density matrix
const double *PB,             // spin dn density matrix
const struct MolGrid_t *grid, // grid point data
const struct gridBasisCenter_t *gB, // grid-basis-center screening
const struct option_t *opt,         // global options
const struct GTOBasis_t *gto,       // basis function info
const struct Molecule_t *mol,       // molecular info
const double cutoff,                // basis product cutoff
const int getMatrix,                // flag to compute XC matrix
const int getForce,                 // flag to compute forces
      double *Exc,            // returned (incremental) energy grid
      double *XA,             // returned spin up matrix
      double *XB,             // returned spin dn matrix
      double *Fx,             // returned (incremental) force
      double *Fy,             // returned (incremental) force
      double *Fz){            // returned (incremental) force

	int i,j;                  // basis function index
	int c;                    // contraction index
	int n;                    // grid point index
	double ir2;               // radius from ith basis
	double dEdrhoa, dEdrhob;  // functional derivative
	double rhoa,  rhob;       // electron densities
	double rhoax, rhoay, rhoaz;  // gradient of spin-up electron densities
	double rhobx, rhoby, rhobz;  // gradient of spin-dn electron densities
	double Ga, Gb;               // divergent term of exchange
	double Gc;                   // divergent term of correlation
	double taua, taub;           // meta-GGA tau for each spin type

	// truncated basis values
	double *VAL;                 // truncated set
	double *VALX, *VALY, *VALZ;  // truncated grad
	double *VALDa, *VALDb;       // truncated div term = G DOT grad(basis)
	int *I2i;                    // index parsing
	int I,J,nSet;                // loop and size

	// forces allocation
	double *fx=NULL, *fy=NULL, *fz=NULL; // forces on basis centers
	double *VALGG=NULL;                  // truncated second derivative
	if(getForce){
		VALGG = calloc(6*nBasis, sizeof(double));
		fx    = calloc(  nBasis, sizeof(double));
		fy    = calloc(  nBasis, sizeof(double));
		fz    = calloc(  nBasis, sizeof(double));
		if(fx==NULL || fy==NULL || fz==NULL || VALGG==NULL){
			printf("getXCDirect - error cannot allocate force array\n");
			exit(-1);
		}
	}

	// allocation
	VAL   = calloc(nBasis, sizeof(double));
	VALX  = calloc(nBasis, sizeof(double));
	VALY  = calloc(nBasis, sizeof(double));
	VALZ  = calloc(nBasis, sizeof(double));
	VALDa = calloc(nBasis, sizeof(double));
	VALDb = calloc(nBasis, sizeof(double));
	I2i   = calloc(nBasis,sizeof(int));
	if(   VAL ==NULL || VALDa==NULL || VALDb==NULL || I2i==NULL
	   || VALX==NULL || VALY ==NULL || VALZ ==NULL){
		printf("getXCDirect - error cannot allocate memory\n");
		exit(-1);
	}

	// main loop
	register double suma, sumb;
	double rad;              // radial part
	double ang;              // angular part
	double radg;             // radial part of gradient
	double angx, angy, angz; // angular part in gradient
	double dx,dy,dz;         // grid point --> basis center

	for(n=0; n < grid->nPoint; n++){

		// compute basis function at this point and truncate set
		for(I=0,i=0; i < nBasis; i++){

			dx = grid->x[n] - gto[i].x0;
			dy = grid->y[n] - gto[i].y0;
			dz = grid->z[n] - gto[i].z0;
			ir2 = dx*dx + dy*dy + dz*dz;
			// screening
			if(ir2 > gB->bR2Max[i]) continue;

			for(rad=0.0, radg=0.0,c=0; c < gto[i].nContract; c++){
				suma  = gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );
				rad  += suma;
				radg -= 2.0 * gto[i].exp [c] * suma;
			}

			ang = pow_int(dx, gto[i].l) *
			      pow_int(dy, gto[i].m) *
			      pow_int(dz, gto[i].n);

			// evaluate gradient at the point
			angx = (gto[i].l<1)?0: gto[i].l * pow_int(dx, gto[i].l-1)
			                                * pow_int(dy, gto[i].m  )
			                                * pow_int(dz, gto[i].n  );

			angy = (gto[i].m<1)?0: gto[i].m * pow_int(dx, gto[i].l  )
			                                * pow_int(dy, gto[i].m-1)
			                                * pow_int(dz, gto[i].n  );

			angz = (gto[i].n<1)?0: gto[i].n * pow_int(dx, gto[i].l  )
			                                * pow_int(dy, gto[i].m  )
			                                * pow_int(dz, gto[i].n-1);

			VAL [I] = ang  * rad;
			VALX[I] = angx * rad + ang * dx * radg;
			VALY[I] = angy * rad + ang * dy * radg;
			VALZ[I] = angz * rad + ang * dz * radg;

			// compute second derivative as needed
			if(getForce){
				double radGG;

				for(radGG=0.0,c=0; c < gto[i].nContract; c++)
					radGG += 4.0 * gto[i].exp [c] * gto[i].exp [c] 
					             * gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );

#define angxx  ( (gto[i].l<2)?0: gto[i].l * (gto[i].l-1) * pow_int(dx, gto[i].l-2)  \
                                                         * pow_int(dy, gto[i].m  )  \
                                                         * pow_int(dz, gto[i].n  )  )

#define angyy  ( (gto[i].m<2)?0: gto[i].m * (gto[i].m-1) * pow_int(dx, gto[i].l  )  \
                                                         * pow_int(dy, gto[i].m-2)  \
                                                         * pow_int(dz, gto[i].n  )  )

#define angzz  ( (gto[i].n<2)?0: gto[i].n * (gto[i].n-1) * pow_int(dx, gto[i].l  )  \
                                                         * pow_int(dy, gto[i].m  )  \
                                                         * pow_int(dz, gto[i].n-2)  )

#define angxy  ( (gto[i].l<1 || gto[i].m<1)?0: gto[i].l * gto[i].m * pow_int(dx, gto[i].l-1)  \
                                                                   * pow_int(dy, gto[i].m-1)  \
                                                                   * pow_int(dz, gto[i].n  )  )

#define angxz  ( (gto[i].l<1 || gto[i].n<1)?0: gto[i].l * gto[i].n * pow_int(dx, gto[i].l-1)  \
                                                                   * pow_int(dy, gto[i].m  )  \
                                                                   * pow_int(dz, gto[i].n-1)  )

#define angyz  ( (gto[i].m<1 || gto[i].n<1)?0: gto[i].m * gto[i].n * pow_int(dx, gto[i].l  )  \
                                                                   * pow_int(dy, gto[i].m-1)  \
                                                                   * pow_int(dz, gto[i].n-1)  )

				VALGG[6*I+0] = dx*dx*ang*radGG + rad*angxx + (ang+2.0*dx*angx)*radg;
				VALGG[6*I+1] = dy*dy*ang*radGG + rad*angyy + (ang+2.0*dy*angy)*radg;
				VALGG[6*I+2] = dz*dz*ang*radGG + rad*angzz + (ang+2.0*dz*angz)*radg;
				VALGG[6*I+3] = dx*dy*ang*radGG + rad*angxy + (dx*angy+dy*angx)*radg;
				VALGG[6*I+4] = dx*dz*ang*radGG + rad*angxz + (dx*angz+dz*angx)*radg;
				VALGG[6*I+5] = dy*dz*ang*radGG + rad*angyz + (dy*angz+dz*angy)*radg;
			}

			// truncate
			I2i [I] = i;

			I++;

		}
		nSet=I;

		// compute electron densities at this point
		rhoa   = 0.0; rhob   = 0.0;
		rhoax  = 0.0; rhoay  = 0.0; rhoaz = 0.0;
		rhobx  = 0.0; rhoby  = 0.0; rhobz = 0.0;

		// compute rho and its gradient
		for(I=0; I < nSet; I++){

			suma = 0.0; sumb = 0.0;
			for(J=0; J < nSet; J++){
				suma  += PA[I2i[I]*nBasis+I2i[J]]*VAL[J];
				sumb  += PB[I2i[I]*nBasis+I2i[J]]*VAL[J];
			}
			rhoa  +=      VAL[I]*suma; rhob  +=      VAL[I]*sumb;
			rhoax += 2.0*VALX[I]*suma; rhobx += 2.0*VALX[I]*sumb;
			rhoay += 2.0*VALY[I]*suma; rhoby += 2.0*VALY[I]*sumb;
			rhoaz += 2.0*VALZ[I]*suma; rhobz += 2.0*VALZ[I]*sumb;

		}

		// enforce positive electron densities
		rhoa = fabs(rhoa); rhob = fabs(rhob);

		// compute energy functional and potential
		dEdrhoa = 0.0; dEdrhob = 0.0;
		Ga      = 0.0; Gb      = 0.0;
		Gc      = 0.0;
		switch(opt->method){

		case METHOD_SVWN5:
			getDFT_xSlater(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			getDFT_cVWN5  (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_HALF:
			getDFT_xSlater(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			dEdrhoa = 0.5*dEdrhoa;
			dEdrhob = 0.5*dEdrhob;
			Exc[n]  = 0.5*Exc[n];
			getDFT_cVWN5(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_HF2xDFT:

			// compute tau for meta-GGA
			taua = 0.0;
			for(i=0; i < nEA; i++){
				double gx=0.0,gy=0.0,gz=0.0;
				for(I=0; I < nSet; I++){
					gx += CA[i*nBasis+I2i[I]]*VALX[I];
					gy += CA[i*nBasis+I2i[I]]*VALY[I];
					gz += CA[i*nBasis+I2i[I]]*VALZ[I];
				}
				taua += (gx*gx + gy*gy + gz*gz)/2.0;
			}
			taub = 0.0;
			for(i=0; i < nEB; i++){
				double gx=0.0,gy=0.0,gz=0.0;
				for(I=0; I < nSet; I++){
					gx += CB[i*nBasis+I2i[I]]*VALX[I];
					gy += CB[i*nBasis+I2i[I]]*VALY[I];
					gz += CB[i*nBasis+I2i[I]]*VALZ[I];
				}
				taub += (gx*gx + gy*gy + gz*gz)/2.0;
			}

			switch(opt->whichHF2xDFT){

			case HF2xDFT_xSLATER:
			getDFT_xGGA(GGAxSLATER,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xMVS:
			getDFT_xMVS(rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            taua, taub,
			            Exc+n);
			break;

			case HF2xDFT_xPY86:
			getDFT_xGGA(GGAxPY86,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xBECKE88:
			getDFT_xGGA(GGAxBECKE88,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xPBE:
			getDFT_xGGA(GGAxPBE,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xCHACHIYO:
			getDFT_xGGA(GGAxCHACHIYO,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;
			}

			// do not use DFT potential
			dEdrhoa = 0.0; dEdrhob = 0.0;
			Ga      = 0.0; Gb      = 0.0;
			Gc      = 0.0;
		break;

		case METHOD_SCHACHIYO:
			getDFT_xSlater   (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			getDFT_cChachiyo (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_CHACHIYO:
			getDFT_xGGA_Chachiyo(rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			getDFT_cGGA_Chachiyo(
			            rhoa, rhob,
			            sqrt((rhoax+rhobx)*(rhoax+rhobx) +
			                 (rhoay+rhoby)*(rhoay+rhoby) +
			                 (rhoaz+rhobz)*(rhoaz+rhobz)),
			            Exc+n, &dEdrhoa, &dEdrhob, &Gc);
		break;

		case METHOD_xBECKE88:
			getDFT_xGGA(GGAxBECKE88,rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
		break;

		default:
			printf("getXCDirect - error unknown DFT METHOD\n");
			exit(-1);
		break;
		}

		// coupled with grid weight
		dEdrhoa *= grid->w[n]; dEdrhob *= grid->w[n];
		Ga      *= grid->w[n]; Gb      *= grid->w[n];
		Gc      *= grid->w[n];

		// prepared G DOT grad(Basis)
		for(I=0; I < nSet; I++){
			double gaDOTgchi = rhoax*VALX[I]+rhoay*VALY[I]+rhoaz*VALZ[I];
			double gbDOTgchi = rhobx*VALX[I]+rhoby*VALY[I]+rhobz*VALZ[I];

			VALDa[I] = Ga*gaDOTgchi + Gc*(gaDOTgchi+gbDOTgchi);
			VALDb[I] = Gb*gbDOTgchi + Gc*(gaDOTgchi+gbDOTgchi);

			//VALDa[I] = Ga*( rhoax*VALX[I]
			//               +rhoay*VALY[I]
			//               +rhoaz*VALZ[I]);
			//VALDb[I] = Gb*( rhobx*VALX[I]
			//               +rhoby*VALY[I]
			//               +rhobz*VALZ[I]);
		}

		// evaluate matrix element
		if(getMatrix)
		for(I=0; I < nSet; I++){

			//suma = dEdrhoa*VAL[I] + VALDa[I];
			//sumb = dEdrhob*VAL[I] + VALDb[I];

			for(J=0; J <= I; J++){
				//XA[I2i[I]*nBasis+I2i[J]] += suma*VAL[J] + VAL[I]*VALDa[J]; 
				//XB[I2i[I]*nBasis+I2i[J]] += sumb*VAL[J] + VAL[I]*VALDb[J];
				XA[I2i[I]*nBasis+I2i[J]] += VAL[I]*dEdrhoa*VAL[J] + VAL[I]*VALDa[J] + VAL[J]*VALDa[I]; 
				XB[I2i[I]*nBasis+I2i[J]] += VAL[I]*dEdrhob*VAL[J] + VAL[I]*VALDb[J] + VAL[J]*VALDb[I];

			}
		}

		// evaluate molecular forces
		if(getForce){
			for(I=0; I < nSet; I++){

				suma = 0.0; sumb = 0.0;
				for(J=0; J < nSet; J++){
					suma  += PA[I2i[I]*nBasis+I2i[J]]*VAL[J];
					sumb  += PB[I2i[I]*nBasis+I2i[J]]*VAL[J];
				}

				// accumulate forces
				fx[I2i[I]] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * VALX[I];
				fy[I2i[I]] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * VALY[I];
				fz[I2i[I]] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * VALZ[I];
				// double negative is hidden in the above expression
				// First  (-) comes from force definition which is negative of gradient
				// Second (-) comes from d(basis)/dX X >> atomic center position
				//            the derivative with respect to the center X is minus
				//            of that with respect to the coordinate x.
				//            i.e. VALX is d(basis)/dx where x is the coordinate
				// Teepanis Chachiyo 8/7/2017

				// divergence term
				double sumg=0.0;
				for(J=0; J < nSet; J++){
/*
					sumg += ( VALX[J]*Ga*rhoax+
				              VALY[J]*Ga*rhoay+
				              VALZ[J]*Ga*rhoaz )*PA[I2i[I]*nBasis+I2i[J]];

					sumg += ( VALX[J]*Gb*rhobx+
				              VALY[J]*Gb*rhoby+
				              VALZ[J]*Gb*rhobz )*PB[I2i[I]*nBasis+I2i[J]];
*/
					sumg += ( VALX[J]*(Ga*rhoax+Gc*rhoax+Gc*rhobx)+
				              VALY[J]*(Ga*rhoay+Gc*rhoay+Gc*rhoby)+
				              VALZ[J]*(Ga*rhoaz+Gc*rhoaz+Gc*rhobz) )*PA[I2i[I]*nBasis+I2i[J]];

					sumg += ( VALX[J]*(Gb*rhobx+Gc*rhoax+Gc*rhobx)+
				              VALY[J]*(Gb*rhoby+Gc*rhoay+Gc*rhoby)+
				              VALZ[J]*(Gb*rhobz+Gc*rhoaz+Gc*rhobz) )*PB[I2i[I]*nBasis+I2i[J]];
				}

#define xx 6*I+0
#define yy 6*I+1
#define zz 6*I+2
#define xy 6*I+3
#define xz 6*I+4
#define yz 6*I+5
/*
				fx[I2i[I]] += 2.0*(sumg*VALX[I] + (Ga*rhoax*suma+Gb*rhobx*sumb)*VALGG[xx] +
				                                  (Ga*rhoay*suma+Gb*rhoby*sumb)*VALGG[xy] +
				                                  (Ga*rhoaz*suma+Gb*rhobz*sumb)*VALGG[xz] );

				fy[I2i[I]] += 2.0*(sumg*VALY[I] + (Ga*rhoax*suma+Gb*rhobx*sumb)*VALGG[xy] +
				                                  (Ga*rhoay*suma+Gb*rhoby*sumb)*VALGG[yy] +
				                                  (Ga*rhoaz*suma+Gb*rhobz*sumb)*VALGG[yz] );

				fz[I2i[I]] += 2.0*(sumg*VALZ[I] + (Ga*rhoax*suma+Gb*rhobx*sumb)*VALGG[xz] +
				                                  (Ga*rhoay*suma+Gb*rhoby*sumb)*VALGG[yz] +
				                                  (Ga*rhoaz*suma+Gb*rhobz*sumb)*VALGG[zz] );
*/
				double x = (Ga*rhoax+Gc*rhoax+Gc*rhobx)*suma + (Gb*rhobx+Gc*rhoax+Gc*rhobx)*sumb;
				double y = (Ga*rhoay+Gc*rhoay+Gc*rhoby)*suma + (Gb*rhoby+Gc*rhoay+Gc*rhoby)*sumb;
				double z = (Ga*rhoaz+Gc*rhoaz+Gc*rhobz)*suma + (Gb*rhobz+Gc*rhoaz+Gc*rhobz)*sumb;

				fx[I2i[I]] += 2.0*(sumg*VALX[I] + x*VALGG[xx] + y*VALGG[xy] + z*VALGG[xz] );
				fy[I2i[I]] += 2.0*(sumg*VALY[I] + x*VALGG[xy] + y*VALGG[yy] + z*VALGG[yz] );
				fz[I2i[I]] += 2.0*(sumg*VALZ[I] + x*VALGG[xz] + y*VALGG[yz] + z*VALGG[zz] );
#undef xx
#undef yy
#undef zz
#undef xy
#undef xz
#undef yz

			}
		}


	}

	// make use of the symmetric property
	if(getMatrix)
	for(i=0; i < nBasis; i++)
	for(j=0; j < i; j++){
		XA[j*nBasis+i] = XA[i*nBasis+j];
		XB[j*nBasis+i] = XB[i*nBasis+j];
	}

	// translate forces to atomic center
	if(getForce){
		for(j=0; j < nBasis; j++){
			for(i=0; i < mol->nAtom; i++)
				if(gto[j].x0==mol->x[i] && 
				   gto[j].y0==mol->y[i] && 
				   gto[j].z0==mol->z[i])
					break;
			Fx[i] += fx[j];
			Fy[i] += fy[j];
			Fz[i] += fz[j];
		}
	}

	// cleanup
	free(VAL);
	free(VALX);
	free(VALY);
	free(VALZ);
	free(VALDa);
	free(VALDb);
	free(I2i);

	if(VALGG!=NULL) free(VALGG);
	if(fx!=NULL)    free(fx);
	if(fy!=NULL)    free(fy);
	if(fz!=NULL)    free(fz);
}



//
// getXCDirect_LibXCcompat: computes exchange-correlation matrix,
// energy functional at grid points, and forces
//
// Aug 17, 2019 - Teepanis Chachiyo
//  Forked from getXCDirect so the the calculation
//  involving the divergent terms is compatible with
//  the "sigma" from LibXC
//
//////// Copied from getXCDirect LOG /////////
//
// getXCDirect: computes exchange-correlation matrix,
// energy functional at grid points, and forces
//
// July 6, 2017 - Teepanis Chachiyo
//   Modified from getXC to compute the basis function
//   at a point instead of storing them in memory.
//   This turned out to be faster and more accurate
//   than the getXC counter part. From tests I suspect
//   the scaling is also linear, but need systematic 
//   test to confirm this observation.
//
//////// Copied from getXC LOG /////////
//
// May 9, 2017 - Teepanis Chachiyo
//   Add gradient and The Chachiyo enhancement factor
//
// August 14, 2016 - Teepanis Chachiyo
//   Add Chachiyo correlation functional
//
// May 11, 2016 - Teepanis Chachiyo
//   Can now compute forces
//
// May 9, 2016 - Teepanis Chachiyo
//   Initial implementation and testing
// 
void getXCDirect_LibXCcompat(
const int nBasis,             // number of basis function
const int nEA,                // number of electron for spin-up
const int nEB,                // number of electron for spin-dn
const double *CA,             // molecular orbitals for spin-up
const double *CB,             // molecular orbitals for spin-dn
const double *PA,             // spin up density matrix
const double *PB,             // spin dn density matrix
const struct MolGrid_t *grid, // grid point data
const struct gridBasisCenter_t *gB, // grid-basis-center screening
const struct option_t *opt,         // global options
const struct GTOBasis_t *gto,       // basis function info
const struct Molecule_t *mol,       // molecular info
const double cutoff,                // basis product cutoff
const int getMatrix,                // flag to compute XC matrix
const int getForce,                 // flag to compute forces
      double *Exc,            // returned (incremental) energy grid
      double *XA,             // returned spin up matrix
      double *XB,             // returned spin dn matrix
      double *Fx,             // returned (incremental) force
      double *Fy,             // returned (incremental) force
      double *Fz){            // returned (incremental) force

	int i,j;                  // basis function index
	int c;                    // contraction index
	int n;                    // grid point index
	double ir2;               // radius from ith basis
	double dEdrhoa, dEdrhob;  // functional derivative
	double rhoa,  rhob;       // electron densities
	double rhoax, rhoay, rhoaz;  // gradient of spin-up electron densities
	double rhobx, rhoby, rhobz;  // gradient of spin-dn electron densities
	double Ga, Gb;               // divergent term of exchange
	double Gc;                   // divergent term of correlation
	double taua, taub;           // meta-GGA tau for each spin type
	double vtaua, vtaub;         // meta-GGA potential

	// truncated basis values
	double *VAL;                 // truncated set
	double *VALX, *VALY, *VALZ;  // truncated grad
	double *VALDa, *VALDb;       // truncated div term = G DOT grad(basis)
	int *I2i;                    // index parsing
	int I,J,nSet;                // loop and size

#ifdef LIBXC
	xc_func_type func1;      // first functional
	xc_func_type func2;      // second functional

	double rho[2];           // input alpha and beta density to libxc
	double sigma[3];         // input sigma to libxc
	double exc[1];           // output energy density from libxc
	double vrho[2];          // output potential
	double vsigma[3];        // output potential the divergent term

	double tau[2];                // input tau to libxc
	double lapl[2] = {0.0, 0.0};  // laplacian is ignored
	double vlapl[2];              // output potential
	double vtau[2];               // output meta-GGA potential 

	if(opt->func1Id) xc_func_init(&func1, opt->func1Id, XC_POLARIZED);
	if(opt->func2Id) xc_func_init(&func2, opt->func2Id, XC_POLARIZED);
#endif

	// forces allocation
	double *fx=NULL, *fy=NULL, *fz=NULL; // forces on basis centers
	double *VALGG=NULL;                  // truncated second derivative
	if(getForce){
		VALGG = calloc(6*nBasis, sizeof(double));
		fx    = calloc(  nBasis, sizeof(double));
		fy    = calloc(  nBasis, sizeof(double));
		fz    = calloc(  nBasis, sizeof(double));
		if(fx==NULL || fy==NULL || fz==NULL || VALGG==NULL){
			printf("getXCDirect_LibXCcompat - error cannot allocate force array\n");
			exit(-1);
		}
	}

	// allocation
	VAL   = calloc(nBasis, sizeof(double));
	VALX  = calloc(nBasis, sizeof(double));
	VALY  = calloc(nBasis, sizeof(double));
	VALZ  = calloc(nBasis, sizeof(double));
	VALDa = calloc(nBasis, sizeof(double));
	VALDb = calloc(nBasis, sizeof(double));
	I2i   = calloc(nBasis,sizeof(int));
	if(   VAL ==NULL || VALDa==NULL || VALDb==NULL || I2i==NULL
	   || VALX==NULL || VALY ==NULL || VALZ ==NULL){
		printf("getXCDirect_LibXCcompat - error cannot allocate memory\n");
		exit(-1);
	}

	// main loop
	register double suma, sumb;
	double rad;              // radial part
	double ang;              // angular part
	double radg;             // radial part of gradient
	double angx, angy, angz; // angular part in gradient
	double dx,dy,dz;         // grid point --> basis center

	for(n=0; n < grid->nPoint; n++){

		// compute basis function at this point and truncate set
		for(I=0,i=0; i < nBasis; i++){

			dx = grid->x[n] - gto[i].x0;
			dy = grid->y[n] - gto[i].y0;
			dz = grid->z[n] - gto[i].z0;
			ir2 = dx*dx + dy*dy + dz*dz;
			// screening
			if(ir2 > gB->bR2Max[i]) continue;

			for(rad=0.0, radg=0.0,c=0; c < gto[i].nContract; c++){
				suma  = gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );
				rad  += suma;
				radg -= 2.0 * gto[i].exp [c] * suma;
			}

			ang = pow_int(dx, gto[i].l) *
			      pow_int(dy, gto[i].m) *
			      pow_int(dz, gto[i].n);

			// evaluate gradient at the point
			angx = (gto[i].l<1)?0: gto[i].l * pow_int(dx, gto[i].l-1)
			                                * pow_int(dy, gto[i].m  )
			                                * pow_int(dz, gto[i].n  );

			angy = (gto[i].m<1)?0: gto[i].m * pow_int(dx, gto[i].l  )
			                                * pow_int(dy, gto[i].m-1)
			                                * pow_int(dz, gto[i].n  );

			angz = (gto[i].n<1)?0: gto[i].n * pow_int(dx, gto[i].l  )
			                                * pow_int(dy, gto[i].m  )
			                                * pow_int(dz, gto[i].n-1);

			VAL [I] = ang  * rad;
			VALX[I] = angx * rad + ang * dx * radg;
			VALY[I] = angy * rad + ang * dy * radg;
			VALZ[I] = angz * rad + ang * dz * radg;

			// compute second derivative as needed
			if(getForce){
				double radGG;

				for(radGG=0.0,c=0; c < gto[i].nContract; c++)
					radGG += 4.0 * gto[i].exp [c] * gto[i].exp [c] 
					             * gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );

#define angxx  ( (gto[i].l<2)?0: gto[i].l * (gto[i].l-1) * pow_int(dx, gto[i].l-2)  \
                                                         * pow_int(dy, gto[i].m  )  \
                                                         * pow_int(dz, gto[i].n  )  )

#define angyy  ( (gto[i].m<2)?0: gto[i].m * (gto[i].m-1) * pow_int(dx, gto[i].l  )  \
                                                         * pow_int(dy, gto[i].m-2)  \
                                                         * pow_int(dz, gto[i].n  )  )

#define angzz  ( (gto[i].n<2)?0: gto[i].n * (gto[i].n-1) * pow_int(dx, gto[i].l  )  \
                                                         * pow_int(dy, gto[i].m  )  \
                                                         * pow_int(dz, gto[i].n-2)  )

#define angxy  ( (gto[i].l<1 || gto[i].m<1)?0: gto[i].l * gto[i].m * pow_int(dx, gto[i].l-1)  \
                                                                   * pow_int(dy, gto[i].m-1)  \
                                                                   * pow_int(dz, gto[i].n  )  )

#define angxz  ( (gto[i].l<1 || gto[i].n<1)?0: gto[i].l * gto[i].n * pow_int(dx, gto[i].l-1)  \
                                                                   * pow_int(dy, gto[i].m  )  \
                                                                   * pow_int(dz, gto[i].n-1)  )

#define angyz  ( (gto[i].m<1 || gto[i].n<1)?0: gto[i].m * gto[i].n * pow_int(dx, gto[i].l  )  \
                                                                   * pow_int(dy, gto[i].m-1)  \
                                                                   * pow_int(dz, gto[i].n-1)  )

				VALGG[6*I+0] = dx*dx*ang*radGG + rad*angxx + (ang+2.0*dx*angx)*radg;
				VALGG[6*I+1] = dy*dy*ang*radGG + rad*angyy + (ang+2.0*dy*angy)*radg;
				VALGG[6*I+2] = dz*dz*ang*radGG + rad*angzz + (ang+2.0*dz*angz)*radg;
				VALGG[6*I+3] = dx*dy*ang*radGG + rad*angxy + (dx*angy+dy*angx)*radg;
				VALGG[6*I+4] = dx*dz*ang*radGG + rad*angxz + (dx*angz+dz*angx)*radg;
				VALGG[6*I+5] = dy*dz*ang*radGG + rad*angyz + (dy*angz+dz*angy)*radg;
			}

			// truncate
			I2i [I] = i;

			I++;

		}
		nSet=I;

		// compute electron densities at this point
		rhoa   = 0.0; rhob   = 0.0;
		rhoax  = 0.0; rhoay  = 0.0; rhoaz = 0.0;
		rhobx  = 0.0; rhoby  = 0.0; rhobz = 0.0;

		// compute rho and its gradient
		for(I=0; I < nSet; I++){

			suma = 0.0; sumb = 0.0;
			for(J=0; J < nSet; J++){
				suma  += PA[I2i[I]*nBasis+I2i[J]]*VAL[J];
				sumb  += PB[I2i[I]*nBasis+I2i[J]]*VAL[J];
			}
			rhoa  +=      VAL[I]*suma; rhob  +=      VAL[I]*sumb;
			rhoax += 2.0*VALX[I]*suma; rhobx += 2.0*VALX[I]*sumb;
			rhoay += 2.0*VALY[I]*suma; rhoby += 2.0*VALY[I]*sumb;
			rhoaz += 2.0*VALZ[I]*suma; rhobz += 2.0*VALZ[I]*sumb;

		}

		// enforce positive electron densities
		rhoa = fabs(rhoa); rhob = fabs(rhob);

		// compute energy functional and potential
		dEdrhoa = 0.0; dEdrhob = 0.0;
		Ga      = 0.0; Gb      = 0.0;
		Gc      = 0.0;
		
		// compute tau for meta-GGA
		/*
		taua = 0.0;
		for(i=0; i < nEA; i++){
			double gx=0.0,gy=0.0,gz=0.0;
			for(I=0; I < nSet; I++){
				gx += CA[i*nBasis+I2i[I]]*VALX[I];
				gy += CA[i*nBasis+I2i[I]]*VALY[I];
				gz += CA[i*nBasis+I2i[I]]*VALZ[I];
			}
			taua += (gx*gx + gy*gy + gz*gz)/2.0;
		}

		taub = 0.0;
		for(i=0; i < nEB; i++){
			double gx=0.0,gy=0.0,gz=0.0;
			for(I=0; I < nSet; I++){
				gx += CB[i*nBasis+I2i[I]]*VALX[I];
				gy += CB[i*nBasis+I2i[I]]*VALY[I];
				gz += CB[i*nBasis+I2i[I]]*VALZ[I];
			}
			taub += (gx*gx + gy*gy + gz*gz)/2.0;
		}
		*/

		taua = 0.0;
		taub = 0.0;
		for(I=0; I < nSet; I++)
		for(J=0; J < nSet; J++){
			double gDOTg = VALX[I]*VALX[J] + VALY[I]*VALY[J] + VALZ[I]*VALZ[J];
			taua += PA[I2i[I]*nBasis+I2i[J]]*gDOTg/2.0;
			taub += PB[I2i[I]*nBasis+I2i[J]]*gDOTg/2.0;
		}

		vtaua = 0.0; vtaub = 0.0;

		switch(opt->method){

		case METHOD_SVWN5:
			getDFT_xSlater(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			getDFT_cVWN5  (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_HALF:
			getDFT_xSlater(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			dEdrhoa = 0.5*dEdrhoa;
			dEdrhob = 0.5*dEdrhob;
			Exc[n]  = 0.5*Exc[n];
			getDFT_cVWN5(rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_HF2xDFT:

			switch(opt->whichHF2xDFT){

			case HF2xDFT_xSLATER:
			getDFT_xGGA(GGAxSLATER,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xMVS:
			getDFT_xMVS(rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            taua, taub,
			            Exc+n);
			break;

			case HF2xDFT_xPY86:
			getDFT_xGGA(GGAxPY86,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xBECKE88:
			getDFT_xGGA(GGAxBECKE88,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xPBE:
			getDFT_xGGA(GGAxPBE,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;

			case HF2xDFT_xCHACHIYO:
			getDFT_xGGA(GGAxCHACHIYO,
			            rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			break;
			}

			// do not use DFT potential
			dEdrhoa = 0.0; dEdrhob = 0.0;
			Ga      = 0.0; Gb      = 0.0;
			Gc      = 0.0;
		break;

		case METHOD_SCHACHIYO:
			getDFT_xSlater   (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
			getDFT_cChachiyo (rhoa, rhob, Exc+n, &dEdrhoa, &dEdrhob);
		break;

		case METHOD_CHACHIYO:
			getDFT_xGGA_Chachiyo(rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
			getDFT_cGGA_Chachiyo(
			            rhoa, rhob,
			            sqrt((rhoax+rhobx)*(rhoax+rhobx) +
			                 (rhoay+rhoby)*(rhoay+rhoby) +
			                 (rhoaz+rhobz)*(rhoaz+rhobz)),
			            Exc+n, &dEdrhoa, &dEdrhob, &Gc);
			Ga += Gc;
			Gb += Gc;
		break;

		case METHOD_xBECKE88:
			getDFT_xGGA(GGAxBECKE88,rhoa, rhob,
			            sqrt(rhoax*rhoax+rhoay*rhoay+rhoaz*rhoaz),
			            sqrt(rhobx*rhobx+rhoby*rhoby+rhobz*rhobz),
			            Exc+n, &dEdrhoa, &dEdrhob, &Ga, &Gb);
		break;

#ifdef LIBXC
		case METHOD_LIBXC:
			// prepare inputs for libxc calls
			Exc[n]   = 0.0;
			dEdrhoa  = 0.0;
			dEdrhob  = 0.0;
			Ga       = 0.0;
			Gb       = 0.0;
			Gc       = 0.0;
			rho[0]   = rhoa;
			rho[1]   = rhob;
			sigma[0] = rhoax*rhoax + rhoay*rhoay + rhoaz*rhoaz;
			sigma[2] = rhobx*rhobx + rhoby*rhoby + rhobz*rhobz;
			sigma[1] = rhoax*rhobx + rhoay*rhoby + rhoaz*rhobz;

			// compute functional 1
			if(opt->func1Id)
			switch(func1.info->family){
				case XC_FAMILY_LDA:
					xc_lda_exc_vxc(&func1, 1, rho, exc, vrho);
					Exc[n]  += exc[0]*(rhoa+rhob);
					dEdrhoa += vrho[0];
					dEdrhob += vrho[1];
					break;
				case XC_FAMILY_GGA:
				case XC_FAMILY_HYB_GGA:
					xc_gga_exc_vxc(&func1, 1, rho, sigma, exc, vrho, vsigma);
					Exc[n]  += exc[0]*(rhoa+rhob);
					dEdrhoa += vrho[0];
					dEdrhob += vrho[1];
					Ga += 2.0*vsigma[0];
					Gb += 2.0*vsigma[2];
					Gc += vsigma[1];
					break;
				case XC_FAMILY_MGGA:
				case XC_FAMILY_HYB_MGGA:
					tau[0] = taua; tau[1] = taub;
					xc_mgga_exc_vxc(&func1, 1, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);
					Exc[n]  += exc[0]*(rhoa+rhob);
					dEdrhoa += vrho[0];
					dEdrhob += vrho[1];
					Ga += 2.0*vsigma[0];
					Gb += 2.0*vsigma[2];
					Gc += vsigma[1];
					vtaua += vtau[0];
					vtaub += vtau[1];
					break;
			}

			// compute functional 2
			if(opt->func2Id)
			switch(func2.info->family){
				case XC_FAMILY_LDA:
					xc_lda_exc_vxc(&func2, 1, rho, exc, vrho);
					Exc[n]  += exc[0]*(rhoa+rhob);
					dEdrhoa += vrho[0];
					dEdrhob += vrho[1];
					break;
				case XC_FAMILY_GGA:
				case XC_FAMILY_HYB_GGA:
					xc_gga_exc_vxc(&func2, 1, rho, sigma, exc, vrho, vsigma);
					Exc[n]  += exc[0]*(rhoa+rhob);
					dEdrhoa += vrho[0];
					dEdrhob += vrho[1];
					Ga += 2.0*vsigma[0];
					Gb += 2.0*vsigma[2];
					Gc += vsigma[1];
					break;
				case XC_FAMILY_MGGA:
				case XC_FAMILY_HYB_MGGA:
					tau[0] = taua; tau[1] = taub;
					xc_mgga_exc_vxc(&func2, 1, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);
					Exc[n]  += exc[0]*(rhoa+rhob);
					dEdrhoa += vrho[0];
					dEdrhob += vrho[1];
					Ga += 2.0*vsigma[0];
					Gb += 2.0*vsigma[2];
					Gc += vsigma[1];
					vtaua += vtau[0];
					vtaub += vtau[1];
					break;
			}
		break;
#endif

		default:
			printf("getXCDirect_LibXCcompat - error unknown DFT METHOD\n");
			exit(-1);
		break;
		}

		// coupled with grid weight
		dEdrhoa *= grid->w[n]; dEdrhob *= grid->w[n];
		Ga      *= grid->w[n]; Gb      *= grid->w[n];
		Gc      *= grid->w[n];
		vtaua   *= grid->w[n]; vtaub   *= grid->w[n];

		// prepared G DOT grad(Basis)
		for(I=0; I < nSet; I++){
			double gaDOTgchi = rhoax*VALX[I]+rhoay*VALY[I]+rhoaz*VALZ[I];
			double gbDOTgchi = rhobx*VALX[I]+rhoby*VALY[I]+rhobz*VALZ[I];

			VALDa[I] = Ga*gaDOTgchi + Gc*gbDOTgchi;
			VALDb[I] = Gb*gbDOTgchi + Gc*gaDOTgchi;

			//VALDa[I] = Ga*gaDOTgchi + Gc*(gaDOTgchi+gbDOTgchi);
			//VALDb[I] = Gb*gbDOTgchi + Gc*(gaDOTgchi+gbDOTgchi);

			//VALDa[I] = Ga*( rhoax*VALX[I]
			//               +rhoay*VALY[I]
			//               +rhoaz*VALZ[I]);
			//VALDb[I] = Gb*( rhobx*VALX[I]
			//               +rhoby*VALY[I]
			//               +rhobz*VALZ[I]);
		}

		// evaluate matrix element
		if(getMatrix)
		for(I=0; I < nSet; I++){

			//suma = dEdrhoa*VAL[I] + VALDa[I];
			//sumb = dEdrhob*VAL[I] + VALDb[I];

			for(J=0; J <= I; J++){
				//XA[I2i[I]*nBasis+I2i[J]] += suma*VAL[J] + VAL[I]*VALDa[J]; 
				//XB[I2i[I]*nBasis+I2i[J]] += sumb*VAL[J] + VAL[I]*VALDb[J];
				XA[I2i[I]*nBasis+I2i[J]] += VAL[I]*dEdrhoa*VAL[J] + VAL[I]*VALDa[J] + VAL[J]*VALDa[I]; 
				XB[I2i[I]*nBasis+I2i[J]] += VAL[I]*dEdrhob*VAL[J] + VAL[I]*VALDb[J] + VAL[J]*VALDb[I];

				// meta-gga contribution
				XA[I2i[I]*nBasis+I2i[J]] += vtaua*(VALX[I]*VALX[J] + VALY[I]*VALY[J] + VALZ[I]*VALZ[J])/2.0;
				XB[I2i[I]*nBasis+I2i[J]] += vtaub*(VALX[I]*VALX[J] + VALY[I]*VALY[J] + VALZ[I]*VALZ[J])/2.0;
			}	
		}

		// evaluate molecular forces
		if(getForce){
			for(I=0; I < nSet; I++){

				suma = 0.0; sumb = 0.0;
				for(J=0; J < nSet; J++){
					suma  += PA[I2i[I]*nBasis+I2i[J]]*VAL[J];
					sumb  += PB[I2i[I]*nBasis+I2i[J]]*VAL[J];
				}

				// accumulate forces
				fx[I2i[I]] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * VALX[I];
				fy[I2i[I]] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * VALY[I];
				fz[I2i[I]] += 2.0*(dEdrhoa*suma+dEdrhob*sumb) * VALZ[I];
				// double negative is hidden in the above expression
				// First  (-) comes from force definition which is negative of gradient
				// Second (-) comes from d(basis)/dX X >> atomic center position
				//            the derivative with respect to the center X is minus
				//            of that with respect to the coordinate x.
				//            i.e. VALX is d(basis)/dx where x is the coordinate
				// Teepanis Chachiyo 8/7/2017

				// divergence term
				double sumg=0.0;
				for(J=0; J < nSet; J++){
/*
					sumg += ( VALX[J]*Ga*rhoax+
				              VALY[J]*Ga*rhoay+
				              VALZ[J]*Ga*rhoaz )*PA[I2i[I]*nBasis+I2i[J]];

					sumg += ( VALX[J]*Gb*rhobx+
				              VALY[J]*Gb*rhoby+
				              VALZ[J]*Gb*rhobz )*PB[I2i[I]*nBasis+I2i[J]];

					sumg += ( VALX[J]*(Ga*rhoax+Gc*rhoax+Gc*rhobx)+
				              VALY[J]*(Ga*rhoay+Gc*rhoay+Gc*rhoby)+
				              VALZ[J]*(Ga*rhoaz+Gc*rhoaz+Gc*rhobz) )*PA[I2i[I]*nBasis+I2i[J]];

					sumg += ( VALX[J]*(Gb*rhobx+Gc*rhoax+Gc*rhobx)+
				              VALY[J]*(Gb*rhoby+Gc*rhoay+Gc*rhoby)+
				              VALZ[J]*(Gb*rhobz+Gc*rhoaz+Gc*rhobz) )*PB[I2i[I]*nBasis+I2i[J]];
*/
					sumg += ( VALX[J]*(Ga*rhoax+Gc*rhobx)+
				              VALY[J]*(Ga*rhoay+Gc*rhoby)+
				              VALZ[J]*(Ga*rhoaz+Gc*rhobz) )*PA[I2i[I]*nBasis+I2i[J]];

					sumg += ( VALX[J]*(Gb*rhobx+Gc*rhoax)+
				              VALY[J]*(Gb*rhoby+Gc*rhoay)+
				              VALZ[J]*(Gb*rhobz+Gc*rhoaz) )*PB[I2i[I]*nBasis+I2i[J]];
				}

#define xx 6*I+0
#define yy 6*I+1
#define zz 6*I+2
#define xy 6*I+3
#define xz 6*I+4
#define yz 6*I+5
/*
				fx[I2i[I]] += 2.0*(sumg*VALX[I] + (Ga*rhoax*suma+Gb*rhobx*sumb)*VALGG[xx] +
				                                  (Ga*rhoay*suma+Gb*rhoby*sumb)*VALGG[xy] +
				                                  (Ga*rhoaz*suma+Gb*rhobz*sumb)*VALGG[xz] );

				fy[I2i[I]] += 2.0*(sumg*VALY[I] + (Ga*rhoax*suma+Gb*rhobx*sumb)*VALGG[xy] +
				                                  (Ga*rhoay*suma+Gb*rhoby*sumb)*VALGG[yy] +
				                                  (Ga*rhoaz*suma+Gb*rhobz*sumb)*VALGG[yz] );

				fz[I2i[I]] += 2.0*(sumg*VALZ[I] + (Ga*rhoax*suma+Gb*rhobx*sumb)*VALGG[xz] +
				                                  (Ga*rhoay*suma+Gb*rhoby*sumb)*VALGG[yz] +
				                                  (Ga*rhoaz*suma+Gb*rhobz*sumb)*VALGG[zz] );
*/
				//double x = (Ga*rhoax+Gc*rhoax+Gc*rhobx)*suma + (Gb*rhobx+Gc*rhoax+Gc*rhobx)*sumb;
				//double y = (Ga*rhoay+Gc*rhoay+Gc*rhoby)*suma + (Gb*rhoby+Gc*rhoay+Gc*rhoby)*sumb;
				//double z = (Ga*rhoaz+Gc*rhoaz+Gc*rhobz)*suma + (Gb*rhobz+Gc*rhoaz+Gc*rhobz)*sumb;

				double x = (Ga*rhoax+Gc*rhobx)*suma + (Gb*rhobx+Gc*rhoax)*sumb;
				double y = (Ga*rhoay+Gc*rhoby)*suma + (Gb*rhoby+Gc*rhoay)*sumb;
				double z = (Ga*rhoaz+Gc*rhobz)*suma + (Gb*rhobz+Gc*rhoaz)*sumb;

				fx[I2i[I]] += 2.0*(sumg*VALX[I] + x*VALGG[xx] + y*VALGG[xy] + z*VALGG[xz] );
				fy[I2i[I]] += 2.0*(sumg*VALY[I] + x*VALGG[xy] + y*VALGG[yy] + z*VALGG[yz] );
				fz[I2i[I]] += 2.0*(sumg*VALZ[I] + x*VALGG[xz] + y*VALGG[yz] + z*VALGG[zz] );

				// meta-gga contribution
				/*
				for(J=0; J < nSet; J++){
					fx[I2i[I]] += vtaua*PA[I2i[I]*nBasis+I2i[J]]*(VALGG[xx]*VALX[J] + VALGG[xy]*VALY[J] + VALGG[xz]*VALZ[J]);
					fx[I2i[I]] += vtaub*PB[I2i[I]*nBasis+I2i[J]]*(VALGG[xx]*VALX[J] + VALGG[xy]*VALY[J] + VALGG[xz]*VALZ[J]);

					fy[I2i[I]] += vtaua*PA[I2i[I]*nBasis+I2i[J]]*(VALGG[xy]*VALX[J] + VALGG[yy]*VALY[J] + VALGG[yz]*VALZ[J]);
					fy[I2i[I]] += vtaub*PB[I2i[I]*nBasis+I2i[J]]*(VALGG[xy]*VALX[J] + VALGG[yy]*VALY[J] + VALGG[yz]*VALZ[J]);

					fz[I2i[I]] += vtaua*PA[I2i[I]*nBasis+I2i[J]]*(VALGG[xz]*VALX[J] + VALGG[yz]*VALY[J] + VALGG[zz]*VALZ[J]);
					fz[I2i[I]] += vtaub*PB[I2i[I]*nBasis+I2i[J]]*(VALGG[xz]*VALX[J] + VALGG[yz]*VALY[J] + VALGG[zz]*VALZ[J]);
				}
				*/

				double sumax = 0.0, sumay = 0.0, sumaz = 0.0;
				double sumbx = 0.0, sumby = 0.0, sumbz = 0.0;
				for(J=0; J < nSet; J++){
					sumax += PA[I2i[I]*nBasis+I2i[J]]*VALX[J];
					sumay += PA[I2i[I]*nBasis+I2i[J]]*VALY[J];
					sumaz += PA[I2i[I]*nBasis+I2i[J]]*VALZ[J];

					sumbx += PB[I2i[I]*nBasis+I2i[J]]*VALX[J];
					sumby += PB[I2i[I]*nBasis+I2i[J]]*VALY[J];
					sumbz += PB[I2i[I]*nBasis+I2i[J]]*VALZ[J];
				}
				fx[I2i[I]] += vtaua*(VALGG[xx]*sumax + VALGG[xy]*sumay + VALGG[xz]*sumaz);
				fx[I2i[I]] += vtaub*(VALGG[xx]*sumbx + VALGG[xy]*sumby + VALGG[xz]*sumbz);

				fy[I2i[I]] += vtaua*(VALGG[xy]*sumax + VALGG[yy]*sumay + VALGG[yz]*sumaz);
				fy[I2i[I]] += vtaub*(VALGG[xy]*sumbx + VALGG[yy]*sumby + VALGG[yz]*sumbz);

				fz[I2i[I]] += vtaua*(VALGG[xz]*sumax + VALGG[yz]*sumay + VALGG[zz]*sumaz);
				fz[I2i[I]] += vtaub*(VALGG[xz]*sumbx + VALGG[yz]*sumby + VALGG[zz]*sumbz);

#undef xx
#undef yy
#undef zz
#undef xy
#undef xz
#undef yz

			}
		}


	}

	// make use of the symmetric property
	if(getMatrix)
	for(i=0; i < nBasis; i++)
	for(j=0; j < i; j++){
		XA[j*nBasis+i] = XA[i*nBasis+j];
		XB[j*nBasis+i] = XB[i*nBasis+j];
	}

	// translate forces to atomic center
	if(getForce){
		for(j=0; j < nBasis; j++){
			for(i=0; i < mol->nAtom; i++)
				if(gto[j].x0==mol->x[i] && 
				   gto[j].y0==mol->y[i] && 
				   gto[j].z0==mol->z[i])
					break;
			Fx[i] += fx[j];
			Fy[i] += fy[j];
			Fz[i] += fz[j];
		}
	}

	// cleanup
	free(VAL);
	free(VALX);
	free(VALY);
	free(VALZ);
	free(VALDa);
	free(VALDb);
	free(I2i);

	if(VALGG!=NULL) free(VALGG);
	if(fx!=NULL)    free(fx);
	if(fy!=NULL)    free(fy);
	if(fz!=NULL)    free(fz);

#ifdef LIBXC
	if(opt->func1Id) xc_func_end(&func1);
	if(opt->func2Id) xc_func_end(&func2);
#endif
}

/////////////////////////////////////////////////////////
//                                                     //
// DFT Functional is based on the magnificent document //
// by Timothy J. Giese, "Density functional theory"   //
// (Version April 15, 2008) which can be downloaded    //
// on the Internet. The PDF version of the document    //
// is also available in this Siam Quantum package.     //
//                                                     //
/////////////////////////////////////////////////////////


//
// getDFT_xSlater: computes Slater exchange functional at grid point
//
// May 10, 2016 - Teepanis Chachiyo
//    Initial implementation and testing
//
void getDFT_xSlater(
double rhoa,     // spin up electron density
double rhob,     // spin dn electron density
double *Exc,     // returned (incremental) exchange-correlation energy
double *dEdrhoa, // returned (incremental) spin up potential
double *dEdrhob){// returned (incremental) spin dn potential

	// parameters
	double alpha = 2./3;
	double Cx = 3./4*pow(3./M_PI,1./3);

	double rho  = rhoa+rhob;
	if(rho > RHO_CUTOFF){

		double z = (rhoa-rhob)/rho;
		double f = (pow(1.+z,4./3) + pow(1.-z,4./3) - 2.)
		           /2/(pow(2.,1./3)-1.);
		double dfdz = 4./3*(pow(1.+z,1./3)-pow(1.-z,1./3))
		                /2/(pow(2.,1./3)-1.);
		double e0 = -3./2*alpha*Cx*pow(rho,1./3);
		double e1 = pow(2.,1./3)*e0;
		double e  = e0 + (e1-e0)*f;
		double delta = e1-e0;

		*dEdrhoa += e + 1./3*e + delta*dfdz * 2. * rhob / rho;
		*dEdrhob += e + 1./3*e - delta*dfdz * 2. * rhoa / rho;
		*Exc     += e * rho;
	}
}

//
// getDFT_cVWN5: computes VWN correlation functional at grid point
//
// May 12, 2016 - Teepanis Chachiyo
//    Initial implementation and testing
//

//parameters from Table IV (T.Giese, 2008)
double VWN5_A[3] = { 0.0621814/2, 0.0621814/4, -1./3/M_PI/M_PI/2};
double VWN5_x[3] = {-0.10498    ,-0.32500    , -0.0047584     };
double VWN5_b[3] = { 3.72744    , 7.06042    ,  1.13107       };
double VWN5_c[3] = {12.9352     ,18.0578     , 13.0045        };

void getDFT_cVWN5(
double rhoa,     // spin up electron density
double rhob,     // spin dn electron density
double *Exc,     // returned (incremental) exchange-correlation energy
double *dEdrhoa, // returned (incremental) spin-up potential
double *dEdrhob){// returned (incremental) spin-dn potential

	double g(double Ai,double xi,double bi,double ci,double rs){
		double x = sqrt(rs);
		double Qi = sqrt(4.*ci-bi*bi);
		double Xix  = x * x +  bi*x  + ci;
		double Xixi = xi* xi + bi*xi + ci;
		double atanQ = atan(Qi/(x+x+bi));
		return Ai*( log(x*x/Xix) + 2.*bi/Qi*atanQ
		            -bi*xi/Xixi*( log((x-xi)*(x-xi)/Xix)
		                        + 2.*(xi+xi+bi)/Qi*atanQ
		                      )
		          );
	}

	double dgdrs(double Ai,double xi,double bi,double ci,double rs){
		double x = sqrt(rs);
		double Qi = sqrt(4.*ci-bi*bi);
		double Xix  = x * x +  bi*x  + ci;
		double Xixi = xi* xi + bi*xi + ci;
		double b2Q2x = bi*bi+Qi*Qi+4.*x*(bi+x);
		double dgdx = 2.*Ai/x 
		              - 4.*Ai*bi/b2Q2x
		              - Ai*bi*xi/Xixi*( 2./(x-xi) 
		                                - 4.*(xi+xi+bi)/b2Q2x
		                              );
		double dgdX = Ai/Xix*(bi*xi/Xixi-1.);
		double dXdx = bi+2*x;
		double dxdrs = 1./2/x;

		return dgdx*dxdrs + dgdX*dXdx*dxdrs;
	}

	double rho=rhoa+rhob;
	if(rho>RHO_CUTOFF){
		double z = (rhoa-rhob)/rho;
		double f = (pow(1.+z,4./3) + pow(1.-z,4./3) - 2.)
		           /2/(pow(2.,1./3)-1.);
		double dfdz = 4./3*(pow(1.+z,1./3)-pow(1.-z,1./3))
		                /2/(pow(2.,1./3)-1.);
		double ddf0 = 4./9/(pow(2.,1./3)-1.);
		double rs = pow(3./4/M_PI/rho,1./3);
		double e0 = g(VWN5_A[0],VWN5_x[0],VWN5_b[0],VWN5_c[0],rs);
		double e1 = g(VWN5_A[1],VWN5_x[1],VWN5_b[1],VWN5_c[1],rs);
		double a  = g(VWN5_A[2],VWN5_x[2],VWN5_b[2],VWN5_c[2],rs);
		double e  = e0 + f*( a/ddf0*(1.-pow(z,4))
		                     +(e1-e0)*pow(z,4)
		                   );
		double dedz = 4.*pow(z,3)*f*(e1-e0-a/ddf0)
		              + dfdz*(a/ddf0*(1.-pow(z,4))
		                      + (e1-e0)*pow(z,4)
		                     );
		double de0drs = dgdrs(VWN5_A[0],VWN5_x[0],VWN5_b[0],VWN5_c[0],rs);
		double de1drs = dgdrs(VWN5_A[1],VWN5_x[1],VWN5_b[1],VWN5_c[1],rs);
		double dadrs  = dgdrs(VWN5_A[2],VWN5_x[2],VWN5_b[2],VWN5_c[2],rs);
		double dedrs  = de0drs
		                + f*( dadrs*(1.-pow(z,4))/ddf0
		                      +(de1drs-de0drs)*pow(z,4)
		                    );
		*dEdrhoa += e + dedrs*(-1./3*rs) + dedz*2.0*rhob/rho;
		*dEdrhob += e + dedrs*(-1./3*rs) - dedz*2.0*rhoa/rho;
		*Exc     += e * rho;
	}

}

void getDFT_cChachiyo(
//
// Teepanis Chachiyo (2016). "Simple and accurate uniform electron gas 
// correlation energy for the full range of densities". J. Chem. Phys. 
// 145 (2): 021101.
//
// T.Chachiyo and H.Chachiyo (2018) "Understanding electron correlation energy
// through density functional theory"
//
double rhoa,     // spin up   electron density
double rhob,     // spin down electron density
double *Ec,      // returned (incremental)  correlation energy
double *dEdrhoa, // returned (incremental) spin up   potential
double *dEdrhob){// returned (incremental) spin down potential

	double rho=rhoa+rhob;
	if(rho>RHO_CUTOFF){

		double z = (rhoa-rhob)/rho;
		double g = (pow(1.+z,2./3)+pow(1.-z,2./3))/2.0;
		double f = 2.0*(1.0-g*g*g);
		double dfdz;
		if(fabs(1.0-z*z)>RHO_CUTOFF) dfdz = -2.0*g*g*(pow(1.+z,-1./3)-pow(1.-z,-1.0/3));
		else                         dfdz =  0.0;

#define Chachiyo_formula(a,b,rs) a*log(1.+b/rs+b/rs/rs)
#define decdrs(a,b,rs)           a/(1. + b/rs + b/rs/rs)*b*(-1./rs/rs-2./rs/rs/rs);
		double rs = pow(3./4/M_PI/rho,1./3);
		double a  = (log(2.)-1.)/M_PI/M_PI/2.;
		double b  = 20.4562557;
		double e0 = Chachiyo_formula(a,b,rs);
		double de0drs = decdrs(a,b,rs);
		       a  = (log(2.)-1.)/M_PI/M_PI/4.;
		       b  = 27.4203609;
		double e1 = Chachiyo_formula(a,b,rs);
		double de1drs = decdrs(a,b,rs)
		double e  = e0 + (e1-e0)*f;

		*dEdrhoa += e - rs/3.0*(de0drs + (de1drs-de0drs)*f) + (e1-e0)*dfdz*2.0*rhob/rho;
		*dEdrhob += e - rs/3.0*(de0drs + (de1drs-de0drs)*f) - (e1-e0)*dfdz*2.0*rhoa/rho;

		*Ec      += e * rho;
	}

}


//
// getDFT_xGGA: computes the GGA exchange at grid point
//
// July 7, 2017 - Teepanis Chachiyo
//    Initial implementation and testing
//
void getDFT_xGGA(
int which,       // which exchange functional
double rhoa,     // spin up electron density
double rhob,     // spin dn electron density
double rhoag,    // magnitude of the gradient of spin up electron density
double rhobg,    // magnitude of the gradient of spin dn electron density
double *Ex,      // returned (incremental) exchange energy
double *dEdrhoa, // returned (incremental) spin up potential
double *dEdrhob, // returned (incremental) spin dn potential
double *Ga,      // returned (incremental) coef of spin up scaled gradient
double *Gb){     // returned (incremental) coef of spin dn scaled gradient


/////////////////////////////////////
////////// Slater Exchange //////////
/////////////////////////////////////

// Cx = 3./4*pow(3./M_PI,1./3);
#define Cx 0.738558766382
#define eLDA(rho) (-Cx*pow(rho,1.0/3))



///////////////////////////////////////////
////////// Becke-88 Exchange //////////////
#define coef_BECKE88 1.0

// 0.0042*2^1/3
#define beta 0.00529166841
// 2^(1/3)
#define TWO13 1.259921049895
#define BH(x) (6.0*beta*Cx*x*asinh(TWO13*x))
#define F_BECKE88(x)      1.0 + beta*x*x/(Cx+BH(x))
#define Fp_BECKE88(x,F)   beta*x/(Cx+BH(x))/(Cx+BH(x))                             \
                          *(2.0*(Cx+BH(x)) - BH(x)                                 \
                                - 6.0*beta*Cx*x*x*TWO13/sqrt(1.0+x*x*TWO13*TWO13)  )



/////////////////////////////////////////////////
///////////// Perdew-Yue 1986 Exchange //////////
/////////////////////////////////////////////////

// 1/2.0/pow(3.0*M_PI*M_PI,1.0/3.0)
#define coef_PY86 0.161620459674
#define F_PY86(x)      pow(1.0+1.296*x*x+14.0*x*x*x*x+0.2*x*x*x*x*x*x,1.0/15.0)
#define Fp_PY86(x,F)   0.0



//////////////////////////////////////////////
////////// Perdew-Burke-Ernzerhof ////////////
//////////////////////////////////////////////

// 1/2.0/pow(3.0*M_PI*M_PI,1.0/3.0)
#define coef_PBE 0.161620459674

#define kappa (0.804)
#define mu (0.21951)
#define F_PBE(x) 1.0+kappa-kappa/(1.0+mu/kappa*x*x)
#define Fp_PBE(x,F) 0.0



//////////////////////////////////////
////////// Chachiyo Exchange /////////
//////////////////////////////////////

// PI square
#define PI2 9.869604401089

// 2/9 (PI/3)^1/3
#define coef_CHACHIYO 0.225664732792

#define F_CHACHIYO(x)    (3.0*x*x+PI2*log(x+1.0))                                 \
                         /(3.0*x+PI2)                                             \
                         /log(x+1.0)

#define Fp_CHACHIYO(x,F) ( 6.0*x*(x+1.0) + PI2                                    \
                           - F*(3.0*x + PI2 + 3.0*(x+1.0)*log(x+1.0))   )         \
                         /(3.0*x+PI2)                                             \
                         /(x+1.0)                                                 \
                         /log(x+1.0)



	// parameters
	double coef,x,e,F,Fp;

	switch(which){
	case GGAxSLATER:   coef = 1.0;           break;
	case GGAxBECKE88:  coef = coef_BECKE88;  break;
	case GGAxPY86:     coef = coef_PY86;     break;
	case GGAxPBE:      coef = coef_PBE;      break;
	case GGAxCHACHIYO: coef = coef_CHACHIYO; break;
	default:
		printf("getDFT_xGGA : error unrecognized functional\n");
		exit(-1);
	}

#define GET_FACTOR                                                     \
	switch(which){                                                     \
	case GGAxSLATER:  F = 1.0;           Fp = 0.0;              break; \
	case GGAxBECKE88: F = F_BECKE88(x);  Fp = Fp_BECKE88(x,F);  break; \
	case GGAxPY86:    F = F_PY86(x);     Fp = Fp_PY86(x,F);     break; \
	case GGAxPBE:     F = F_PBE(x);      Fp = Fp_PBE(x,F);      break; \
	case GGAxCHACHIYO:F = F_CHACHIYO(x); Fp = Fp_CHACHIYO(x,F); break; \
	default:          F = 1.0;           Fp = 0.0;              break; \
	}

	if(rhoa > RHO_CUTOFF){
		e = eLDA(2.0*rhoa);
		x = coef*(2.0*rhoag)/pow(2.0*rhoa,4.0/3);
		if(x < RHO_CUTOFF){ F = 1.0; Fp = 0.0; } else{ GET_FACTOR; }
		*Ex      += rhoa  * e * F;
		*dEdrhoa += 4.0/3 * e * (F-x*Fp);
		if(rhoag > RHO_CUTOFF) *Ga += (-Cx)*coef*Fp/rhoag;
	}

	if(rhob > RHO_CUTOFF){
		e = eLDA(2.0*rhob);
		x = coef*(2.0*rhobg)/pow(2.0*rhob,4.0/3);
		if(x < RHO_CUTOFF){ F = 1.0; Fp = 0.0; } else{ GET_FACTOR; }
		*Ex      += rhob  * e * F;
		*dEdrhob += 4.0/3 * e * (F-x*Fp);
		if(rhobg > RHO_CUTOFF) *Gb += (-Cx)*coef*Fp/rhobg;
	}

#undef GET_FACTOR

#undef Cx
#undef eLDA

#undef F_BECKE88
#undef F_BECKE88
#undef coef_BECKE88
#undef BH
#undef TWO13
#undef beta

#undef F_PY86
#undef Fp_PY86
#undef coef_PY86

#undef F_PBE
#undef Fp_PBE
#undef coef_PBE
#undef mu
#undef kappa

#undef F_CHACHIYO
#undef Fp_CHACHIYO
#undef coef_CHACHIYO
#undef PI2
}

//
// getDFT_xMVS: computes the meta-GGA Made Very Simple exchange.
// The potential is not supported.
//
// September 20, 2017 - Teepanis Chachiyo
//    Initial implementation and testing
//
void getDFT_xMVS(
double rhoa,     // spin up electron density
double rhob,     // spin dn electron density
double rhoag,    // magnitude of the gradient of spin up electron density
double rhobg,    // magnitude of the gradient of spin dn electron density
double taua,     // meta-GGA tau for spin-up density
double taub,     // meta-GGA tau for spin-dn density
double *Ex){     // returned (incremental) exchange energy

// Dirac Exchange
// Cx = 3./4*pow(3./M_PI,1./3);
#define Cx 0.738558766382
#define eLDA(rho) (-Cx*pow(rho,1.0/3))

double e,F,s;
double f, alpha, alpha2, tauW, tauU;

// 1/2.0/pow(3.0*M_PI*M_PI,1.0/3.0)
#define coef 0.161620459674
#define k0 0.174
#define e1 (-1.6665)
#define c1 0.7438
#define b 0.0233

#define GET_FACTOR(rho,rhog,tau)                                       \
			tauW = (2.0*rhog)*(2.0*rhog)/8.0/(2.0*rho);                \
			tauU = 3.0/10*pow(3.0*M_PI*M_PI,2.0/3)*pow(2.0*rho,5.0/3); \
			alpha = (2.0*tau - tauW)/tauU;                             \
			alpha2 = alpha*alpha;                                      \
			f = (1.0-alpha)/pow(1.0                                    \
			                    +2.0*e1*alpha2                         \
			                    + e1*e1*alpha2*alpha2                  \
			                    +    c1*alpha2*alpha2,1.0/4);          \
			F = (1.0 + f*k0)/pow(1.0+b*s*s*s*s,1.0/8);                 \

	if(rhoa > RHO_CUTOFF){
		e = eLDA(2.0*rhoa);
		s = coef*(2.0*rhoag)/pow(2.0*rhoa,4.0/3);
		if(s < RHO_CUTOFF) F = 1.0; else{ GET_FACTOR(rhoa,rhoag,taua); }
		*Ex      += rhoa  * e * F;
	}

	if(rhob > RHO_CUTOFF){
		e = eLDA(2.0*rhob);
		s = coef*(2.0*rhobg)/pow(2.0*rhob,4.0/3);
		if(s < RHO_CUTOFF) F = 1.0; else{ GET_FACTOR(rhob,rhobg,taub); }
		*Ex      += rhob  * e * F;
	}

#undef coef
#undef k0
#undef e1
#undef c1
#undef b
#undef Cx
#undef eLDA
}

// T.Chachiyo and H.Chachiyo (2018) "Understanding electron correlation energy
// through density functional theory"
//
void getDFT_xGGA_Chachiyo(
double rhoa,     // spin up   electron density
double rhob,     // spin down electron density
double rhoag,    // spin up   density gradient 
double rhobg,    // spin down density gradient
double *Ex,      // returned (incremental) exchange energy
double *dEdrhoa, // returned (incremental) spin up   potential
double *dEdrhob, // returned (incremental) spin down potential
double *Ga,      // returned (incremental) coef. of alpha density gradient
double *Gb){     // returned (incremental) coef. of beta  density gradient

#define PI2  M_PI*M_PI
#define eDirac(rho) -3./4*pow(3./M_PI*rho,1./3)
#define GET_F    F    = (3.0*x*x+PI2*log(x+1.0))/(3.0*x+PI2)/log(x+1.0)
#define GET_dFdx dFdx = ( 6.0*x*(x+1.0) + PI2                         \
                          - F*(3.0*x + PI2 + 3.0*(x+1.0)*log(x+1.0)) )\
                        /(3.0*x+PI2)/(x+1.0)/log(x+1.0)

	double e,x,F,dFdx;

	// spin up contribution
	if(rhoa > RHO_CUTOFF){
		e = eDirac(2.0*rhoa);
		x = (2.0*rhoag)/pow(2.0*rhoa,4./3) * 2.0/9*pow(M_PI/3.,1./3); 
		if(x<RHO_CUTOFF){ F = 1.0; dFdx = 0.0; } else{ GET_F; GET_dFdx; }
		*Ex      += rhoa  * e * F;
		*dEdrhoa += 4.0/3 * e * (F-x*dFdx);
		if(rhoag > RHO_CUTOFF) *Ga += -dFdx/rhoag/6.0;
	}

	// spin down contribution
	if(rhob > RHO_CUTOFF){
		e = eDirac(2.0*rhob);
		x = (2.0*rhobg)/pow(2.0*rhob,4./3) * 2.0/9*pow(M_PI/3.,1./3);
		if(x<RHO_CUTOFF){ F = 1.0; dFdx = 0.0; } else{ GET_F; GET_dFdx; }
		*Ex      += rhob  * e * F;
		*dEdrhob += 4.0/3 * e * (F-x*dFdx);
		if(rhobg > RHO_CUTOFF) *Gb += -dFdx/rhobg/6.0;
	}
}

// T.Chachiyo and H.Chachiyo (2018) "Understanding electron correlation energy
// through density functional theory"
//
void getDFT_cGGA_Chachiyo(
double rhoa,     // spin up   electron density
double rhob,     // spin down electron density
double rhog,     // total gradient
double *Ec,      // returned (incremental) correlation energy
double *dEdrhoa, // returned (incremental) spin up   potential
double *dEdrhob, // returned (incremental) spin down potential
double *G){      // returned (incremental) coef. of total density gradient

	double rho  = rhoa  + rhob;
	if(rho > RHO_CUTOFF){

		// From the Chachiyo correlation functional
		double va=0., vb=0., e_unif=0.;
		getDFT_cChachiyo(rhoa, rhob, &e_unif, &va, &vb);
		e_unif = e_unif/rho;

		double t = rhog/pow(rho,7.0/6) * pow(M_PI/3.0, 1.0/6)/4.0;
		double h = 0.06672632;
		double S = pow(1.0+t*t,h/e_unif);

		*dEdrhoa += S*(va-7.0/3*h*t*t/(1.0+t*t)-h*log(1.0+t*t)*(va/e_unif-1.0));
		*dEdrhob += S*(vb-7.0/3*h*t*t/(1.0+t*t)-h*log(1.0+t*t)*(vb/e_unif-1.0));

		if(rhog*rhog > RHO_CUTOFF)
			*G += 2.0*h*S*t*t/(1.0+t*t)*rho/rhog/rhog;

		*Ec += rho * S * e_unif;
	}
}
