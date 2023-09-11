/*
	Siam Quantum [SQ] is a program pacakage that performs quantum
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basis.h"
#include "matrix.h"
#include "lin.h"
#include "int.h"
#include "fgamma.h"
#include "util.h"
#include "option.h"
#include "rhf.h"
#include "uhf.h"
#include "conv.h"
#include "check.h"
#include "rpc.h"
#include "grid.h"
#include "dft.h"

// normalizeC : normalized eigen vector 
//
// Jan 15, 2014 - Teepanis Chachiyo
//	Can select how many orbital to normalize
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static void normalizeC(int nBasis, int nOcc, double *S, double *C){
	double sum;
	int i,j,k;

	for(k=0; k < nOcc; k++){
		sum = 0.0;
		for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++)
				sum += C[k*nBasis+i]
				      *C[k*nBasis+j]
				      *S[i*nBasis+j];
		sum = 1.0/sqrt(sum);
		for(i=0; i < nBasis; i++)
			C[k*nBasis+i] = sum * C[k*nBasis+i];
	}
	return;
}

// uhf_getDMatrix : compute density matrix
//
// July 10, 2010 - Teepanis Chachiyo
//    Migrate to uhf scheme
//
// 2008 - Teepanis Chachiyo
//    Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//    Not using static anymore
//
void uhf_getDMatrix(int nBasis, int nOcc, double *C, double *P){
	int i,j,k;
	double sum;

	for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			sum = 0.0;
			for(k=0; k < nOcc; k++){
				sum += C[k*nBasis+i] * C[k*nBasis+j];
			}
			P[i*nBasis+j] = sum;
		}
	return;
}

// uhf_getGMatrix : compute G matrix element for both alpha and beta spin
//
// May 6, 2016 - Teepanis Chachiyo
//  Use JT,KA,KB instead of GA,GB
//
// Jan 16, 2014 - Teepanis Chachiyo
//	Not using static in the declaration anymore
//
// Mar 6, 2013 - Teepanis Chachiyo
//  Use parallel version
//
// Nov 19, 2012 - Teepanis Chachiyo
//  Passing cutoff value as an argument
//
// July 10, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculations
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
void uhf_getGMatrix(
	int nBasis,
	const struct GTOBasis_t *gto,
	const double *Schwarz,
	double cutoff,
	const double *PA, const double *PB,
	double *JT, double *KA, double *KB,
	struct option_t *opt){

	int i,j;
	
	// reset to zero
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		JT[i*nBasis+j] = 0.0;
		KA[i*nBasis+j] = 0.0;
		KB[i*nBasis+j] = 0.0;
	}

	/////////////////////
	// parallel version
	/////////////////////
	double *JTset, *KAset, *KBset; // set of JT,KA,KB for each cpu
	int *status;                   // status for each cpu
	int alldone;                   // all idle flag
	int n;                         // cpu counter

	// allocate memory
	status=calloc(opt->nCPU,sizeof(int));
	JTset =calloc(nBasis*nBasis*opt->nCPU,sizeof(double));
	KAset =calloc(nBasis*nBasis*opt->nCPU,sizeof(double));
	KBset =calloc(nBasis*nBasis*opt->nCPU,sizeof(double));
	if(status==NULL || JTset==NULL || KAset==NULL || KBset==NULL){
		printf("uhf_getGMatrix - error cannot allocate memory\n");
		exit(-1);
	}

	// reset status to idle
	for(n=(opt->nCPU-1);n>=0;n--) status[n] = RPC_IDLE;

	// loop thru all cpu and compute GA and GB
	do{
		// remote call to all cpu except childID=0
		for(n=(opt->nCPU-1);n>0;n--)
			if(status[n] != RPC_DONE)
			status[n] = rpc_GTO_JK_Matrix_Quartet_Parallel(status[n],n,nBasis,
			                                               PA, PB, gto, Schwarz,
			                                               cutoff,
			                                               JTset+nBasis*nBasis*n,
			                                               KAset+nBasis*nBasis*n,
			                                               KBset+nBasis*nBasis*n,
			                                               opt);

		// local call for childID=0
		if(status[0]==RPC_IDLE){
			GTO_JK_Matrix_Quartet_Parallel(0, nBasis, PA, PB, gto, Schwarz, cutoff,
			                               JTset, KAset, KBset, opt);
			status[0]=RPC_DONE;
		}

		// check if all done
		alldone=1;
		for(n=(opt->nCPU-1);n>=0;n--) if(status[n] != RPC_DONE) alldone=0;

	}while(!alldone);

	// accumulate GA and GB
	for(n=(opt->nCPU-1);n>=0;n--){
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			JT[i*nBasis+j] += JTset[nBasis*nBasis*n +i*nBasis+j];
			KA[i*nBasis+j] += KAset[nBasis*nBasis*n +i*nBasis+j];
			KB[i*nBasis+j] += KBset[nBasis*nBasis*n +i*nBasis+j];
		}
	}

	// clean memory
	free(status);
	free(JTset);
	free(KAset);
	free(KBset);

	///////////////////////////////////
	// perform G matrix computation  //
	///////////////////////////////////
	//GTO_JK_Matrix_Quartet(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB, opt);
	//GTO_JK_Matrix_ShellSet(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB, opt);
	//GTO_JK_Matrix_NoSymm(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB);
	//GTO_JK_Matrix_PrimeSpace(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB);
	//GTO_JK_Matrix_PrimeSpaceShellSet(nBasis, PA, PB, gto, Schwarz, cutoff, GA, GB);
}

// getEtotal : compute total energy this is equal to
// electronic part + nuclei part
//
// Aug 18, 2019 - Teepanis Chachiyo
//  Add support for LibXC: GGA and Hybrid
//
// May 6, 2016 - Teepanis Chachiyo
//  Use J,KA,KB instead
//
// March 16, 2016 - Teepanis Chachiyo
//  Add uniform electric field
//
// July 12, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculations
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static double uhf_getEtotal(
	const int nBasis,
	const struct Molecule_t *mol,
	const struct option_t *opt,
	const double *PA,
	const double *PB,
	const double *H,
	const double *JT,
	const double *KA,
	const double *KB,
	const struct MolGrid_t *grid,
	const double *Exc){

	double E=0.0;
	int i,j;

	// compute nuclei - uniform electric field energy
	for(i=0; i < mol->nAtom; i++){
		E -= mol->Z[i] * mol->x[i] * opt->Ex;
		E -= mol->Z[i] * mol->y[i] * opt->Ey;
		E -= mol->Z[i] * mol->z[i] * opt->Ez;
	}

	// compute nuclei repulsion energy
	E += nuclei_coulomb(mol);

	/*
	// include electron energy given in (Szabo and Ostlund)
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		E += 0.5*((H[i*nBasis+j]+FA[i*nBasis+j])*PA[i*nBasis+j] +
			      (H[i*nBasis+j]+FB[i*nBasis+j])*PB[i*nBasis+j]);
	}
	*/

	//
	// electronic energy
	//

	// core and Coulomb
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		E += (H[i*nBasis+j]  + 0.5*JT[i*nBasis+j])
		    *(PA[i*nBasis+j] + PB[i*nBasis+j]);
	}

	// exchange correlation
	switch(opt->method){
	case METHOD_HF:
	case METHOD_MP2:
	case METHOD_HF2xDFT:
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			E += (0.5*KA[i*nBasis+j])*PA[i*nBasis+j] +
			     (0.5*KB[i*nBasis+j])*PB[i*nBasis+j];
		}
	break;

	case METHOD_SVWN5:
	case METHOD_SCHACHIYO:
	case METHOD_CHACHIYO:
	case METHOD_xBECKE88:
		for(i=0; i < grid->nPoint; i++)
			E += Exc[i] * grid->w[i];
	break;

	case METHOD_HALF:
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			E += 0.5*(0.5*KA[i*nBasis+j])*PA[i*nBasis+j] +
			     0.5*(0.5*KB[i*nBasis+j])*PB[i*nBasis+j];
		}

		for(i=0; i < grid->nPoint; i++)
			E += Exc[i] * grid->w[i];
	break;

#ifdef LIBXC
	case METHOD_LIBXC:
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			E += opt->hyb_hf_coef*(0.5*KA[i*nBasis+j])*PA[i*nBasis+j] +
			     opt->hyb_hf_coef*(0.5*KB[i*nBasis+j])*PB[i*nBasis+j];
		}

		for(i=0; i < grid->nPoint; i++)
			E += Exc[i] * grid->w[i];
	break;
#endif

	}

	return E;
}

// uhf_rho : computes and return electron density at
// specified Cartesian point (x,y,z)
//
// August 15, 2010 Nanta Sophonrat
//  Fix bug when computing electron density. Forgot to multiply by 2
//
// July 12, 2010 - Teepanis Chachiyo
//  Migrate to unrestricted calculation
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//  The subroutine now takes density matrix as an argument
//  as supposed to molecular orbital coefficient. This helps
//  increase speed.
//
// March 18, 2010 - Teepanis Chachiyo
//  Store eval_chi() in memory to reduce the number of calculations
//  from nBasis*nBasis to nBasis only. Also exploit P[] symmetrical
//  nature to reduce calculation time by a factor of two.
//
double uhf_rho(int nBasis,              // number of basis function
               struct GTOBasis_t *gto,  // function structure
               double *PA, double *PB,  // density matrix
               double x, double y, double z){

	int i,j;
	double sum=0.0;
	double *chi;

/* //////// very inefficient //////////

	// evaluate density
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		sum += P[i*nBasis+j]*eval_chi(i,gto,x,y,z)
		                    *eval_chi(j,gto,x,y,z);
	}
///////////////////////////////////// */

	// allocate memory
	chi=calloc(nBasis, sizeof(double));
	if(chi==NULL){
		printf("uhf_rho - error cannot allocate memory\n");
		exit(-1);
	}

	// evaluate chi at a point
	for(i=0; i < nBasis; i++) chi[i] = eval_chi(i,gto,x,y,z);

	//
	// evaluate density
	//

	// off diagonal
	for(i=0; i < nBasis; i++)
	for(j=0; j < i; j++)
		sum += (PA[i*nBasis+j]+PB[i*nBasis+j])*chi[i]*chi[j];
	sum = 2.0*sum;

	// diagonal
	for(i=0; i < nBasis; i++)
		sum += (PA[i*nBasis+i]+PB[i*nBasis+i])*chi[i]*chi[i];

	// clean up memory
	free(chi);

	return sum;
}


//
// uhf_rho_XYPlane : compute electron density for all points on a grid
// on xy-plane. This should be faster than compute each point one-by-one.
//
// Oct 8, 2012 - Teepanis Chachiyo
//     Initial implementationa from uhf_rho(...) and testing
//
void uhf_rho_XYPlane(
	int nBasis,                         // number of basis function
    const struct GTOBasis_t *gto,       // basis function structure
    const double *PA, const double *PB, // density matrix
	double cutoff,                      // ignore prefactor
	double x0, double y0, double z0,    // origin of the grid
	double dx, double dy,               // step size
	int nx, int ny,                     // number of points
	double *rhoXY){                     // returned matrix of size ny*nx

	int p,q;                  // basis function index
	int cP, cQ;               // contracted function index
	int i,j;                  // grid point index

	double *X, *Y, Z;         // buffer in x, y, and z directions
	double x,y;               // current coordinates
	double prefactor;         // prefactor
	double rab2;              // distance between two basis function
	double x1,x2,y1,y2,z1,z2; // center of each basis functions

	// allocation
	if((X=calloc(nx,sizeof(double)))==NULL ||
	   (Y=calloc(ny,sizeof(double)))==NULL){
		printf("uhf_rho_XYPlane : error cannot allocate memory\n");
		exit(-1);
	}

	// reset values
	for(j=0; j < ny; j++)
	for(i=0; i < nx; i++)
		rhoXY[j*nx+i] = 0.0;

	// loop through all basis functions
	for(p=0; p < nBasis; p++)
	for(q=0; q <= p; q++){

		// set center variables
		x1 = gto[p].x0; y1 = gto[p].y0; z1 = gto[p].z0;
		x2 = gto[q].x0; y2 = gto[q].y0; z2 = gto[q].z0;
#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+\
                                  (y1-y2)*(y1-y2)+\
                                  (z1-z2)*(z1-z2))
		rab2 = DIST2(x1,y1,z1,x2,y2,z2);

		// loop through contractions
		for(cP=0; cP < gto[p].nContract; cP++)
		for(cQ=0; cQ < gto[q].nContract; cQ++){

			// screening
			prefactor = exp(-gto[p].exp[cP]*gto[q].exp[cQ]*rab2/
			                (gto[p].exp[cP]+gto[q].exp[cQ]));
			if(fabs(prefactor) < cutoff) continue;

			// compute prefactor which is indenpendent of points
			prefactor = (PA[p*nBasis+q] + PB[p*nBasis+q]) *
			            gto[p].norm[cP] * gto[q].norm[cQ] *
			            gto[p].coef[cP] * gto[q].coef[cQ];

			// apply symmetry
			if(p!=q) prefactor = prefactor+prefactor;

			// compute Z specific value
			Z = pow_int(z0-gto[p].z0,gto[p].n)*pow_int(z0-gto[q].z0,gto[q].n)*
			    exp( - gto[p].exp[cP] * (z0-gto[p].z0) * (z0-gto[p].z0)
			         - gto[q].exp[cQ] * (z0-gto[q].z0) * (z0-gto[q].z0));

			// absorb Z into prefactor
			prefactor = prefactor * Z;

			// compute Y specfic values
			for(j=0; j < ny; j++){
				y    = y0+j*dy;
				Y[j] = pow_int(y-gto[p].y0,gto[p].m)*pow_int(y-gto[q].y0,gto[q].m)*
				       exp( - gto[p].exp[cP] * (y-gto[p].y0) * (y-gto[p].y0)
				            - gto[q].exp[cQ] * (y-gto[q].y0) * (y-gto[q].y0));
			}

			// compute X specific values
			for(i=0; i < nx; i++){
				x    = x0+i*dx;
				X[i] = pow_int(x-gto[p].x0,gto[p].l)*pow_int(x-gto[q].x0,gto[q].l)*
				       exp( - gto[p].exp[cP] * (x-gto[p].x0) * (x-gto[p].x0)
				            - gto[q].exp[cQ] * (x-gto[q].x0) * (x-gto[q].x0));
			}

			// loop through all points on xy-plane
			for(j=0; j < ny; j++)
			for(i=0; i < nx; i++)
				rhoXY[j*nx+i] += prefactor * X[i] * Y[j];
		}
	}

	// free memory
	free(X);
	free(Y);
}


//
// uhf_potential_XYPlane : compute electric potential for all points on a grid
// on xy-plane. This should be faster than compute each point one-by-one.
//
// Oct 9, 2012 - Teepanis Chachiyo
//     Initial implementationa from uhf_rho_XYPlane(...) and testing
//
void uhf_potential_XYPlane(
	int nBasis,                         // number of basis function
	const struct GTOBasis_t *gto,       // basis function structure
	const struct Molecule_t *mol,       // molecule structure 
	const double *PA, const double *PB, // density matrix
	double cutoff,                      // ignore prefactor
	double x0, double y0, double z0,    // origin of the grid
	double dx, double dy,               // step size
	int nx, int ny,                     // number of points
	double *phiXY){                     // returned matrix of size ny*nx

	int p,q;             // basis function index
	int cP, cQ;          // contracted function index
	int i,j;             // grid point index
	int A;               // atom index
	double rr;           // distance from point to nuclei
	double r[3];         // position of nucleus

	double prefactor;    // prefactor

#define MAXL 8
	// nai integral sections
	double eta,xp,yp,zp,sum,rab2,rcp2;
	double *Ax,*Ay,Az[2*MAXL],*tAx,*tAy;
	double F[2*MAXL];
	int kX, kY, kZ;
	double x1,y1,z1,x2,y2,z2;
	int l1,m1,n1,l2,m2,n2;

	// allocation
	if((Ax=calloc(nx*2*MAXL,sizeof(double)))==NULL ||
	   (Ay=calloc(ny*2*MAXL,sizeof(double)))==NULL){
		printf("uhf_potential_XYPlane : error cannot allocate memory\n");
		exit(-1);
	}

	// reset values
	for(j=0; j < ny; j++)
	for(i=0; i < nx; i++)
		phiXY[j*nx+i] = 0.0;

	// loop through all basis functions
	for(p=0; p < nBasis; p++)
	for(q=0; q <= p; q++){

		// set angular index
		l1 = gto[p].l; m1 = gto[p].m; n1 = gto[p].n;
		l2 = gto[q].l; m2 = gto[q].m; n2 = gto[q].n;

		// set center variables
		x1 = gto[p].x0; y1 = gto[p].y0; z1 = gto[p].z0;
		x2 = gto[q].x0; y2 = gto[q].y0; z2 = gto[q].z0;

		rab2 = DIST2(x1,y1,z1,x2,y2,z2);

		// loop through contractions
		for(cP=0; cP < gto[p].nContract; cP++)
		for(cQ=0; cQ < gto[q].nContract; cQ++){

			// prepare for nai1d
			eta = 1.0/(gto[p].exp[cP]+gto[q].exp[cQ]);

			// prefactor and screening
			prefactor = exp(-gto[p].exp[cP]*gto[q].exp[cQ]*rab2*eta);

			// screening
			if(fabs(prefactor) < cutoff) continue;

			// compute prefactor which is indenpendent of points
			prefactor = prefactor                         * 
			            (PA[p*nBasis+q] + PB[p*nBasis+q]) *
			            gto[p].norm[cP] * gto[q].norm[cQ] *
			            gto[p].coef[cP] * gto[q].coef[cQ] *
	                    2.0 * M_PI * eta;

			// apply symmetry
			if(p!=q) prefactor = prefactor+prefactor;

			// compute center of two-gaussian function
			xp   = (gto[p].exp[cP]*x1+gto[q].exp[cQ]*x2)*eta;
			yp   = (gto[p].exp[cP]*y1+gto[q].exp[cQ]*y2)*eta;
			zp   = (gto[p].exp[cP]*z1+gto[q].exp[cQ]*z2)*eta;

			// compute Z specific value
			naiA_1d(Az,n1,n2,zp-z1,zp-z2,zp-z0,eta);

			// compute Y specfic values
			for(j=0; j < ny; j++)
				naiA_1d(Ay+j*2*MAXL,m1,m2,yp-y1,yp-y2,yp-y0-j*dy,eta);

			// compute X specific values
			for(i=0; i < nx; i++)
				naiA_1d(Ax+i*2*MAXL,l1,l2,xp-x1,xp-x2,xp-x0-i*dx,eta);

			// loop through all points on xy-plane
			for(j=0; j < ny; j++)
			for(i=0; i < nx; i++){

				rcp2 = DIST2(x0+i*dx,y0+j*dy,z0,xp,yp,zp);
				fgamma_set(l1+l2+m1+m2+n1+n2,rcp2/eta,F);

				sum = 0.0;
				tAx = Ax+i*2*MAXL;
				tAy = Ay+j*2*MAXL;
				for(kX=0; kX < l1+l2+1; kX++)
				for(kY=0; kY < m1+m2+1; kY++)
				for(kZ=0; kZ < n1+n2+1; kZ++)
					sum += tAx[kX]*tAy[kY]*Az[kZ]*F[kX+kY+kZ];

				phiXY[j*nx+i] -= prefactor * sum;
			}
		}
	}

	// include potential from nuclei
	for(j=0; j < ny; j++)
	for(i=0; i < nx; i++)
	for(A=0; A<mol->nAtom; A++){
		r[0]  = x0+i*dx - mol->x[A];
		r[1]  = y0+j*dy - mol->y[A];
		r[2]  = z0      - mol->z[A];

		rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

		if(rr == 0.0) continue;
		else phiXY[j*nx+i] += mol->Z[A]/rr;
	}

	// free memory
	free(Ax);
	free(Ay);
}


// uhf_mo : computes and return molecular orbital at
// specified Cartesian point (x,y,z)
//
// July 12, 2010 - Teepanis Chachiyo
//    Migrate to unrestricted calculation
//
// Dec 2009 - Theerapon Khamla
//    Initial implementation
//
// Mar 05, 2010 - Teepanis Chachiyo
//    Bug fix, invalid range of molecular orbital should be
//    "n >= nBasis" not "n > nBasis" as before.
//
double uhf_mo(int nBasis,              // number of basis function
              struct GTOBasis_t *gto,  // function structure
              double *C,               // molecular orbital
              int n,                   // orbital index (1st orbital is zero)
              double x, double y, double z){
	int i;
	double sum;

	// sanitize
	if(n >= nBasis || n < 0){
		printf("uhf_mo: error invalid molecular orbital index\n");
		exit(EXIT_FAILURE);
	}

	// evaluate molecular orbital
	sum = 0.0;
	for(i=0; i < nBasis; i++){
		sum += C[n*nBasis+i]*eval_chi(i,gto,x,y,z);
	}

	return sum;

}

// uhf_potential : computes and return electric potential at the specified 
// Cartesian point (x,y,z). It also includes the contribution due to the
// nuclei in the system, except when the point (x,y,z) is right on top of
// a nucleus, the contribution from that particular nucleus is ignored.
//
// Oct 26, 2010 - Nanta Sophonrat
//         Initial Implementation     
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//         Added to Siam Quantum source tree
//
double uhf_potential(int nBasis,                    // number of basis function
                     const struct GTOBasis_t *gto,  // function structure
                     const struct Molecule_t *mol,  // molecular structure info
                     const double *PA,              // alpha density matrix
                     const double *PB,              // beta density matrix
                     double x, double y, double z){

	int    i=0,j=0;     // loop index
	double rr;          // distance between nucleus and point(x1, x2, x3)
	double pot = 0.0;   // electric potential
	double r[3];        // distances between atoms

	// multiply by 2 due to the symmetry
	// off-diagonal elements
	for(i=1; i<nBasis; i++)
	for(j=0; j<i; j++){
		pot += 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])
		          * GTO_nai(i, j, gto, x, y, z);
	}

	// diagonal elements
	for(i=0; i<nBasis; i++)
		pot += (PA[i*nBasis+i]+PB[i*nBasis+i])*GTO_nai(i, i, gto, x, y, z);

	// include potential from nuclei
	for(i=0; i<mol->nAtom; i++){
		r[0]  = x - mol->x[i];
		r[1]  = y - mol->y[i];
		r[2]  = z - mol->z[i];

		rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

		if(rr == 0.0) continue;
		else pot += mol->Z[i]/rr;
	}
	return pot;
}

// uhf_efield: compute electric field at a specified point (x,y,z). It also 
// includes the nuclei contribution, except when the point (x,y,z) is right
// on top of the nucles, the contribution from that particular nucleus is
// ignored.
//
// Nov 22, 2010 - Nanta Sophonrat
//         Initial Implementation     
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//         Added to Siam Quantum source tree
//
void uhf_efield(int nBasis,                          // number of basis function
                const struct GTOBasis_t *gto,        // function structure
                const struct Molecule_t *mol,        // molecular structure info
                const double *PA,                    // alpha density matrix
                const double *PB,                    // beta density matrix
                double x, double y, double z,        // specified point
                double *ex, double *ey, double *ez){ // return field

	int    i=0,j=0;       // basis function index
	double rr;            // distance between the nucleus and (x,y,z)
	double r[3];          // displacement between the nucleus and (x,y,z)
	double tex, tey, tez; // temporary field

	// initial value
	*ex = 0.0; *ey = 0.0; *ez = 0.0;

	// multiply by 2 due to the symmetry
	// off-diagonal elements
	for(i=1; i<nBasis; i++)
	for(j=0; j<i; j++){
			GTO_efi(i,j,gto,x,y,z,&tex,&tey,&tez);
			*ex -= 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])*tex;
			*ey -= 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])*tey;			
			*ez -= 2.0*(PA[i*nBasis+j]+PB[i*nBasis+j])*tez;
	}

	// diagonal elements
	for(i=0; i<nBasis; i++){
		GTO_efi(i,i,gto,x,y,z,&tex,&tey,&tez);
		*ex -= (PA[i*nBasis+i]+PB[i*nBasis+i])*tex;
		*ey -= (PA[i*nBasis+i]+PB[i*nBasis+i])*tey;
		*ez -= (PA[i*nBasis+i]+PB[i*nBasis+i])*tez;
	}

	// include nuclei contribution
	for(i=0; i<mol->nAtom; i++){
		r[0]  = x - mol->x[i];
		r[1]  = y - mol->y[i];
		r[2]  = z - mol->z[i];

		rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
		rr = rr*rr*rr;

		if(rr == 0.0) continue;
		else{
			*ex += r[0]*mol->Z[i]/rr;
			*ey += r[1]*mol->Z[i]/rr;
			*ez += r[2]*mol->Z[i]/rr;
		}
	}
}


double uhf_GradRhoOverRho43(
int nBasis,              // number of basis function
struct GTOBasis_t *gto,  // function structure
double *PA, double *PB,  // density matrix
double x, double y, double z){

	int i,j,c;
	double sum;
	double *chi,*chix,*chiy,*chiz;
	double dx,dy,dz,ir2,ang,rad,radg,angx,angy,angz;

	// allocate memory
	chi =calloc(nBasis, sizeof(double));
	chix=calloc(nBasis, sizeof(double));
	chiy=calloc(nBasis, sizeof(double));
	chiz=calloc(nBasis, sizeof(double));
	if(chi==NULL || chix==NULL || chiy==NULL || chiz==NULL){
		printf("uhf_sparam - error cannot allocate memory\n");
		exit(-1);
	}

	// evaluate chi and its gradient at a point
	for(i=0; i < nBasis; i++){

		dx  = x - gto[i].x0;
		dy  = y - gto[i].y0;
		dz  = z - gto[i].z0;
		ir2 = dx*dx + dy*dy + dz*dz;

		for(rad=0.0, radg=0.0,c=0; c < gto[i].nContract; c++){
			sum  = gto[i].norm[c] * gto[i].coef[c] * exp( - gto[i].exp[c] * ir2 );
			rad  += sum;
			radg -= 2.0 * gto[i].exp [c] * sum;
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

		chi [i] = ang  * rad;
		chix[i] = angx * rad + ang * dx * radg;
		chiy[i] = angy * rad + ang * dy * radg;
		chiz[i] = angz * rad + ang * dz * radg;
	}

	//
	// evaluating densities
	//
	double suma, sumb;
	double rhoa   = 0.0, rhob   = 0.0;
	double rhoax  = 0.0, rhoay  = 0.0, rhoaz = 0.0;
	double rhobx  = 0.0, rhoby  = 0.0, rhobz = 0.0;

	// compute rho and its gradient
	for(i=0; i < nBasis; i++){

		suma = 0.0; sumb = 0.0;
		for(j=0; j < nBasis; j++){
			suma  += PA[i*nBasis+j]*chi[j];
			sumb  += PB[i*nBasis+j]*chi[j];
		}
		rhoa  +=      chi[i]*suma; rhob  +=      chi[i]*sumb;
		rhoax += 2.0*chix[i]*suma; rhobx += 2.0*chix[i]*sumb;
		rhoay += 2.0*chiy[i]*suma; rhoby += 2.0*chiy[i]*sumb;
		rhoaz += 2.0*chiz[i]*suma; rhobz += 2.0*chiz[i]*sumb;
	}

	// enforce positive value to unsure numerical stability
	rhoa = fabs(rhoa);
	rhob = fabs(rhob);

	sum = sqrt((rhoax+rhobx)*(rhoax+rhobx)+
	           (rhoay+rhoby)*(rhoay+rhoby)+
	           (rhoaz+rhobz)*(rhoaz+rhobz));
	sum = sum/pow(rhoa+rhob,4.0/3.0);

	// clean up memory
	free(chi);
	free(chix);
	free(chiy);
	free(chiz);

	return sum;
}


// uhf : carries out Unrestricted Hartree-Fock calculation until
// convergence is reached. It returns total electronic energy of the
// systems.
//
// If the scf convergence is not reached, it returns zero instead.
//
// In this version, the e-e integral is calculated using direct
// integration every time it is needed. They are not stored in
// the memory, which is not a good idea!
//
// Aug 18, 2019 - Teepanis Chachiyo
//     - Use small grid for early DFT scf cycles
//     - Call LibXC compatible XC subroutines
//
// May 6, 2016 - Teepanis Chachiyo
//    - Use JT,KA,KB instead of GA,GB
//
// March 14, 2016 - Teepanis Chachiyo
//    - Add uniform electric field
//
// Dec 2015 - Teepanis Chachiyo
//    - Report Kinetic Matrix, Nuclei Matrix, Overlap explicitly
//
// Feb 2013 - Teepanis Chachiyo
//    - No longer need to rebuild fock matrix when switching accuracy
//      This is only for the case when the initial accuracy is <= 10^-3
//
// Jan 26, 2013 - Teepanis Chachiyo
//    - Rewrite the scf convergence structure
//
// Jan 9, 2013 - Teepanis Chachiyo
//    - Forgot to check avgdP criteria when switching accuracy
//
// Nov 26, 2012 - Teepanis Chachiyo
//    - Using 3 steps integral accuracy
//    - Check point is saved even if the run is not converged
//    - Add option to save checkpoint every scf cycle
//
// Nov 19, 2012 - Teepanis Chachiyo
//    - Schwarz cutoff is calculated according to convergence requirement
//      and the density matrix
//    - do not report initial non-interacting electronic energy
// 
// Oct 21, 2010 - Teepanis Chachiyo
//    Use option SCFConv to set convergence threshold
//    The function now returns zero if the convergence failed.
//
// Sep 6, 2010 - Teepanis Chachiyo
//    Use density matrix difference (insted of absolute density matrix value)
//    to compute the change in G matrix
//    Use idensity matrix as initial density matrix
//
// Sep 6, 2010 - Teepanis Chachiyo
//    Take the convergence part outside "uhf" function
//
// July 12, 2010 - Teepanis Chachiyo
//    Migrate from RHF to UHF
//
// 2008 - Teepanis Chachiyo
// 	  Initial implementation
//
// Dec 14, 2008 - Teepanis Chachiyo
// 	  Simplify and use Direct integral for the first release
//	  version.
//
// March 11, 2010 - Teepanis Chachiyo
//    Add reading and saving density matrix
//
// May 21, 2011 - Teepanis Chachiyo
//    Add diis option
//
// Mar 6, 2013 - Teepanis Chachiyo
//     Use parallel version
//
// Mar 19, 2013 - Teepanis Chachiyo
//     Revamp convergence checking
//
// Jan 16, 2014 - Teepanis Chachiyo
//		Call someEigen instead of full spectrum of eigen value during
//		iterations. But still need to do full at the last step because
//		it is needed for post-scf
//
double uhf(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nEA,                 // total number of spin up electrons
	int nEB,                 // total number of spin down electrons
	double *CA,              // returned molecular alpha spin orbital
	double *CB,              // returned molecular beta spin orbital 
	double *eA,              // returned eigen values
	double *eB,              // returned eigen values
	struct option_t *opt){   // global option

	double Etot=0.0;      // total electronic energy
	double dE=0.0;        // energy change
	double avgdP=0.0;     // average change in density matrix
	double gamma=1.0;     // update coefficient
	double sum=0.0;       // generic summation variable
	double realCutoff;    // schwarz cutoff
	double thisCutoff;    // current cutoff
	double cutoffA=0.0;   // first stage cutoff 
	double cutoffB=0.0;   // intermediate cutoff
	double cutoffC=0.0;   // final stage cutoff
	double sumA=0.0;      // sum of alpha density matrix
	double sumB=0.0;      // sum of beta density matrix

	// matrix elements
	double *JT=NULL;     // total Coulomb matrix
	double *dJT=NULL;    // change in Coulomb matrix
	double *KA=NULL;     // exchange matrix
	double *dKA=NULL;    // change in exchange matrix
	double *KB=NULL;     // exchange matrix
	double *dKB=NULL;    // change in exchange matrix
	double *PA=NULL;     // density matrix
	double *dPA=NULL;    // change in density matrix
	double *FA=NULL;     // fock matrix
	double *PB=NULL;     // density matrix
	double *dPB=NULL;    // change in density matrix
	double *FB=NULL;     // fock matrix
	double *H=NULL;      // h core matrix
	double *T=NULL;      // kinetic matrix
	double *V=NULL;      // nuclei potential matrix
	double *S=NULL;      // overlap matrix
	double *Schwarz=NULL;// Schwarz upper bound matrix

	// DFT specific
	struct MolGrid_t *grid=NULL;      // numerical grid
	struct MolGrid_t *gridSmall=NULL; // small grid for early scf cycles
	struct MolGrid_t *thisGrid=NULL;  // pointer to current grid during scf
	struct gridBasisCenter_t *gB=NULL;// grid-basis-center screening
	//struct BGrid_t *bg=NULL;        // basis function grid
	double *Exc=NULL;            // exchange-correlation functional
	double *XCa=NULL;            // spin up exchange-correlation matrix
	double *XCb=NULL;            // spin dn exchange-correlation matrix

	int i,j,iter=0;
	int notConverged=1;

	FILE *fd;            // file descriptor for density matrix

	// report
	switch(opt->method){
	case METHOD_HF:
	case METHOD_MP2:
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----            SOLVING HARTREE-FOCK EQUATION          -----\n"
	"-------------------------------------------------------------\n"
	);
	break;

	default:
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----            SOLVING KOHN-SHAM DFT EQUATION         -----\n"
	"-------------------------------------------------------------\n"
	);
	break;
	}
	fflush(stdout);

	// report
	if(! opt->restricted){
		printf("Requested unrestricted molecular orbitals\n");
		printf("There are %d alpha spin and %d beta spin electrons\n", nEA, nEB);
	}

	if(opt->restricted){
		printf("Requested restricted molecular orbitals\n");
		printf("There are %d electrons in the density matrix\n", nEA+nEB);
		if(nEA!=nEB){
			printf("uhf - error number of electron in each spin is not the same\n");
			exit(-1);
		}
	}
	printMolecule_XYZ(mol, stdout);

#undef  ALLOCATE
#define ALLOCATE(P)                                      \
P = calloc(nBasis*nBasis, sizeof(double));               \
if(P==NULL){                                             \
	printf("uhf: Error - Cannot allocate memory\n"); \
	exit(EXIT_FAILURE);                              \
}

	// memory allocation
	ALLOCATE(JT);
	ALLOCATE(dJT);
	ALLOCATE(PA);
	ALLOCATE(KA);
	ALLOCATE(dPA);
	ALLOCATE(dKA);
	ALLOCATE(FA);
	ALLOCATE(PB);
	ALLOCATE(KB);
	ALLOCATE(dPB);
	ALLOCATE(dKB);
	ALLOCATE(FB);
	ALLOCATE(H);
	ALLOCATE(T);
	ALLOCATE(V);
	ALLOCATE(S);

#undef  ALLOCATE

	/////////////////////////////////////////////
	// Building necessary matrix elements     ///
	/////////////////////////////////////////////

	// report
	printf(
	"Computing 1-electron matrix elements .");fflush(stdout);

	// get kinetic matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		T[i*nBasis+j] = GTO_kinetic(i,j,gto);
		// symmetrize matrix
		T[j*nBasis+i] = T[i*nBasis+j];
	}

	// get nuclei matrix
	printf(".");fflush(stdout);
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		V[i*nBasis+j] = GTO_nuclei(i,j,gto,mol);
		// symmetrize matrix
		V[j*nBasis+i] = V[i*nBasis+j];
	}

	printf(".");fflush(stdout);
	// get overlap matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		S[i*nBasis+j] = GTO_overlap(i,j,gto);
		// symmetrize matrix
		S[j*nBasis+i] = S[i*nBasis+j];
	}
	printf("\n");fflush(stdout);

	// add uniform electric field term
	if(opt->Ex != 0.0 || opt->Ey != 0.0 || opt->Ez != 0.0){
		printf("Uniform electric field: (%.4f,%.4f,%.4f) AU\n",
			opt->Ex, opt->Ey, opt->Ez);
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			V[i*nBasis+j] += GTO_moment(i,j,gto,1,0,0,0.0,0.0,0.0) * opt->Ex;
			V[i*nBasis+j] += GTO_moment(i,j,gto,0,1,0,0.0,0.0,0.0) * opt->Ey;
			V[i*nBasis+j] += GTO_moment(i,j,gto,0,0,1,0.0,0.0,0.0) * opt->Ez;
		}
	}

	// build Hcore 
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		H[i*nBasis+j] = T[i*nBasis+j] + V[i*nBasis+j];
	}

	//
	// use the provided molecular coefficients
	//
	if(opt->SCFGuess == SCFGUESS_CACB){
		printf("Use provided orbitals as initial density matrix ...\n");
		fflush(stdout);
		normalizeC(nBasis, nEA, S, CA);
		normalizeC(nBasis, nEB, S, CB);
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);
	}
	
	//
	// Diagonalize core hamiltonian to guess density matrix
	//
	if(opt->SCFGuess == SCFGUESS_CORE){
		printf("Diagonalizing H for initial density matrix ...\n");
		fflush(stdout);
		gen_sym_eigen(nBasis, H, S, eA, CA);
		gen_sym_eigen(nBasis, H, S, eB, CB);
		normalizeC(nBasis, nEA, S, CA);
		normalizeC(nBasis, nEB, S, CB);
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);
	}

	//
	// set density matrix to diagonal as a guess
	//
	if(opt->SCFGuess == SCFGUESS_DIAG){
		printf("Use identity matrix as initial density matrix ...\n");
		fflush(stdout);
		sum=0.0;
		for(i=0; i<nBasis; i++) sum+=S[i*nBasis+i];
		for(i=0; i<nBasis; i++)
		for(j=0; j<nBasis; j++){
			if(i==j){
				PA[i*nBasis+j] = nEA/sum;
				PB[i*nBasis+j] = nEB/sum;
			}else{
				PA[i*nBasis+j] = 0.0;
				PB[i*nBasis+j] = 0.0;
			}
		}
	}

	//
	// read orbital from checkpoint
	//
	if(opt->SCFGuess == SCFGUESS_CHECK){
		fflush(stdout);
		printf("Read orbital from checkpoint for initial density matrix ...\n");
		guess_checkpoint_orbital(nBasis, gto, mol, opt, S, CA, CB);
		normalizeC(nBasis, nEA, S, CA);
		normalizeC(nBasis, nEB, S, CB);
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);
	}

	// load density matrix if requested
	if(opt->loadDMatrix){

		// openfile
		fd=fopen(opt->DMatrixFile,"r");
		if(fd==NULL){
			printf("uhf - error cannot open file %s for reading\n", 
			       opt->DMatrixFile);
			exit(-1);
		}

		//////////////////
		// loop and read
		//////////////////
		// handle RHF case
		if(opt->restricted){
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				if(fscanf(fd, "%lf", PA+i*nBasis+j)!=1){
					printf("uhf - error reading density matrix\n");
					exit(-1);
				}
				PB[i*nBasis+j] = PA[i*nBasis+j];
			}
		}
		// handle UHF case - read 2 times
		if(! opt->restricted){
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				if(fscanf(fd, "%lf", PA+i*nBasis+j)!=1){
					printf("uhf - error reading alpha density matrix\n");
					exit(-1);
				}
			}
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				if(fscanf(fd, "%lf", PB+i*nBasis+j)!=1){
					printf("uhf - error reading beta spin density matrix\n");
					exit(-1);
				}
			}
		}

		// close file
		fclose(fd);

		// report status
		printf("Density matrix loaded from %s\n", opt->DMatrixFile);
	}

	// compute cutoff
#define GETCUTOFF(accuracy) ((accuracy)/(sumA*sumA+sumA*sumB+sumB*sumB)/50.0)
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		sumA += PA[i*nBasis+j];
		sumB += PB[i*nBasis+j];
	}
	realCutoff = GETCUTOFF(opt->SCFConv);
	opt->SCFCutoff = realCutoff;

	// 2-electron integral information
	printf("Processing 2E integrals ...\n");
	printf("Schwarz inequality screening cut-off %.8E\n", opt->SCFCutoff);
	printf("Primitive prefactor cut-off %0.8E\n", PRIMITIVE_CUTOFF);
	Schwarz = create_Schwarz(nBasis, gto);

	//
	// generate grid for DFT
	//
	switch(opt->method){
	case METHOD_SVWN5:
	case METHOD_HALF:
	case METHOD_HF2xDFT:
	case METHOD_SCHACHIYO:
	case METHOD_CHACHIYO:
	case METHOD_xBECKE88:

#ifdef LIBXC
	case METHOD_LIBXC:
#endif

		// build full molecular grid
		printf("Basis function grid cut-off %0.8E\n", BGRID_CUTOFF);

		// always use small grid for initial stage of scf cycle
		gridSmall = genMolGrid(opt->whichGrid, opt->nradialInit, opt->nlabedevInit,
		                       nBasis, gto, mol, BGRID_CUTOFF);

		// desired grid point
		grid = genMolGrid(opt->whichGrid, opt->nradial, opt->nlabedev,
		                  nBasis, gto, mol, BGRID_CUTOFF);
		printf("Number of numerical grid points %d\n", grid->nPoint);
		gB = genGridBasisCenter(nBasis, gto, mol, grid);
		//bg = genBasisGrid(nBasis, gto, grid, BGRID_CUTOFF);
		//printf("Number of stored basis grid points %ld\n", bg->nPoint);

		// allocation
		Exc     = calloc(grid->nPoint, sizeof(double));
		XCa     = calloc(nBasis*nBasis, sizeof(double));
		XCb     = calloc(nBasis*nBasis, sizeof(double));
		if(Exc==NULL || XCa==NULL  || XCb==NULL){
			printf("uhf: Error - Cannot allocate memory\n");
			exit(EXIT_FAILURE);
		}
	break;
	}
	fflush(stdout);

	// scf loop
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   printf("Use 4-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DIIS3:   printf("Use 3-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DIIS2:   printf("Use 2-Point DIIS convergence method\n"); break;
	case CONVMETHOD_DAMPING: printf("Use simple weighting convergence method\n"); break;
	default:
		printf("uhf - error no specific request for convergence method\n");
		exit(-1);
	break;
	}

	printf("Drag coefficient %f\n", opt->SCFDrag);
	printf("SCFConv %.2E\n", opt->SCFConv);
	printf("SCFMax %d iterations\n", opt->SCFMax);
	printf("Enter SCF loop ... \n");
	printf(
	"Iteration  Total Energy [Hartrees]  RMSD Density\n"
	"------------------------------------------------\n");
	fflush(stdout);
	
	// oscillation drag coefficient
	gamma = opt->SCFDrag;

	// call convergence function for the first time to initialize it
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   conv_diis4(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DIIS3:   conv_diis3(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DIIS2:   conv_diis2(nBasis,   gamma, PA, PB); break;
	case CONVMETHOD_DAMPING: conv_damping(nBasis, gamma, PA, PB); break;
	default:
		printf("uhf - error unknown opt->convMethod\n");
		exit(-1);
	break;
	}

	// preparations
	for(i=0; i<nBasis; i++)
	for(j=0; j<nBasis; j++){

		// set Coulomb and Exchange matrix to zero
		JT[i*nBasis+j] = 0.0;
		KA[i*nBasis+j] = 0.0;
		KB[i*nBasis+j] = 0.0;

		// set delta density to the initial density
		dPA[i*nBasis+j] = PA[i*nBasis+j];
		dPB[i*nBasis+j] = PB[i*nBasis+j];
	}

	// manage integral accuracy
	switch(opt->SCFAccuracy){
	case SCFACCURACY_1STEP:
		thisCutoff = realCutoff;
		thisGrid   = grid;
	break;
	case SCFACCURACY_3STEP:
		cutoffA    = GETCUTOFF(SCFACCURACY_3STEP_A);
		cutoffB    = GETCUTOFF(SCFACCURACY_3STEP_B);
		cutoffC    = GETCUTOFF(SCFACCURACY_3STEP_C);
		thisCutoff = cutoffA;
		thisGrid   = gridSmall;
	break;
	default:
		printf("uhf - error unknown SCFAccuracy\n");
		exit(-1);
	break;
	}

	do{
		iter++;

		// compute delta J,KA,KB matrix
		if(opt->restricted)
			uhf_getGMatrix(nBasis, gto, Schwarz, thisCutoff, dPA, dPA, dJT, dKA, dKB, opt);
		else
			uhf_getGMatrix(nBasis, gto, Schwarz, thisCutoff, dPA, dPB, dJT, dKA, dKB, opt);

		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){

			// update J,K matrix
			JT[i*nBasis+j] = JT[i*nBasis+j] + dJT[i*nBasis+j];
			KA[i*nBasis+j] = KA[i*nBasis+j] + dKA[i*nBasis+j];
			KB[i*nBasis+j] = KB[i*nBasis+j] + dKB[i*nBasis+j];

			// saving current density matrix to dPA and dPB
			dPA[i*nBasis+j] = PA[i*nBasis+j];
			dPB[i*nBasis+j] = PB[i*nBasis+j];
		}

		// build fock matrix
		switch(opt->method){
		case METHOD_HF:
		case METHOD_MP2:
		case METHOD_HF2xDFT:
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				FA[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + KA[i*nBasis+j];
				FB[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + KB[i*nBasis+j];
			}
		break;

		case METHOD_SVWN5:
		case METHOD_SCHACHIYO:
		case METHOD_CHACHIYO:
		case METHOD_xBECKE88:
			// build exchange-correlation matrix
			memset(Exc, 0, grid->nPoint *sizeof(double));
			memset(XCa, 0, nBasis*nBasis*sizeof(double));
			memset(XCb, 0, nBasis*nBasis*sizeof(double));
			//getXC(nBasis, PA, PB, grid, bg, opt, gto, mol,
			//      thisCutoff,1,0,
			//      Exc, XCa, XCb,
			//      NULL, NULL, NULL);
			//getXCDirect(nBasis, nEA, nEB, CA, CB, PA, PB, thisGrid, gB, opt, gto, mol,
			//            thisCutoff,1,0,
			//            Exc, XCa, XCb,
			//            NULL, NULL, NULL);
			getXCDirect_LibXCcompat(nBasis, nEA, nEB, CA, CB, PA, PB, thisGrid, gB, opt, gto, mol,
			            thisCutoff,1,0,
			            Exc, XCa, XCb,
			            NULL, NULL, NULL);

			// compute fock matrix
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				FA[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + XCa[i*nBasis+j];
				FB[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + XCb[i*nBasis+j];
			}
		break;

		case METHOD_HALF:
			// build exchange-correlation matrix
			memset(Exc, 0, grid->nPoint *sizeof(double));
			memset(XCa, 0, nBasis*nBasis*sizeof(double));
			memset(XCb, 0, nBasis*nBasis*sizeof(double));
			//getXC(nBasis, PA, PB, grid, bg, opt, gto, mol,
			//      thisCutoff,1,0,
			//      Exc, XCa, XCb,
			//      NULL, NULL, NULL);
			//getXCDirect(nBasis, nEA, nEB, CA, CB, PA, PB, thisGrid, gB, opt, gto, mol,
			//            thisCutoff,1,0,
			//            Exc, XCa, XCb,
			//            NULL, NULL, NULL);
			getXCDirect_LibXCcompat(nBasis, nEA, nEB, CA, CB, PA, PB, thisGrid, gB, opt, gto, mol,
			            thisCutoff,1,0,
			            Exc, XCa, XCb,
			            NULL, NULL, NULL);

			// compute fock matrix
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				FA[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + 0.5*KA[i*nBasis+j] + XCa[i*nBasis+j];
				FB[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + 0.5*KB[i*nBasis+j] + XCb[i*nBasis+j];
			}
		break;

#ifdef LIBXC
		case METHOD_LIBXC:
			// build exchange-correlation matrix
			memset(Exc, 0, grid->nPoint *sizeof(double));
			memset(XCa, 0, nBasis*nBasis*sizeof(double));
			memset(XCb, 0, nBasis*nBasis*sizeof(double));
			getXCDirect_LibXCcompat(nBasis, nEA, nEB, CA, CB, PA, PB, thisGrid, gB, opt, gto, mol,
			            thisCutoff,1,0,
			            Exc, XCa, XCb,
			            NULL, NULL, NULL);

			// compute fock matrix
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				FA[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + opt->hyb_hf_coef*KA[i*nBasis+j] + XCa[i*nBasis+j];
				FB[i*nBasis+j] = H[i*nBasis+j] + JT[i*nBasis+j] + opt->hyb_hf_coef*KB[i*nBasis+j] + XCb[i*nBasis+j];
			}
		break;
#endif

		}

		// solve generalized eigen value problem and normalize orbital
		gen_sym_SomeEigen(nBasis, nEA, FA, S, eA, CA);
		gen_sym_SomeEigen(nBasis, nEB, FB, S, eB, CB);
		normalizeC(nBasis, nEA, S, CA);
		normalizeC(nBasis, nEB, S, CB);

		// get new P matrix
		uhf_getDMatrix(nBasis, nEA, CA, PA);
		uhf_getDMatrix(nBasis, nEB, CB, PB);

		// compute energy and energ difference
		dE     = Etot;
		Etot   = uhf_getEtotal(nBasis, mol, opt, PA, PB, H, JT, KA, KB, thisGrid, Exc);
		dE     = Etot - dE;

		// update P matrix using convergence method
		switch(opt->convMethod){
		case CONVMETHOD_DIIS4:   avgdP = conv_diis4(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DIIS3:   avgdP = conv_diis3(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DIIS2:   avgdP = conv_diis2(nBasis,   gamma, PA, PB); break;
		case CONVMETHOD_DAMPING: avgdP = conv_damping(nBasis, gamma, PA, PB); break;
		}

		// check convergence
		notConverged = fabs(dE) > opt->SCFConv || avgdP    > opt->SCFConv;

		// compute delta density matrix for the next step
		for(i=0; i<nBasis; i++)
		for(j=0; j<nBasis; j++){
			dPA[i*nBasis+j] = PA[i*nBasis+j] - dPA[i*nBasis+j];
			dPB[i*nBasis+j] = PB[i*nBasis+j] - dPB[i*nBasis+j];
		}

		printf(" %5d %20.8f %20.4E\n", iter, Etot, avgdP);
		
		// check if we have reached scfmax limit
		if(iter >= opt->SCFMax) break;

		// flush output
		fflush(stdout);

		// manage integral accuracy
		if(opt->SCFAccuracy == SCFACCURACY_3STEP){
			
			// check convergence
			if(!notConverged) break;

			// check if we need to switch accuracy
			if((thisCutoff == cutoffA && fabs(dE) <= SCFACCURACY_3STEP_A && avgdP <= SCFACCURACY_3STEP_A) ||
			   (thisCutoff == cutoffB && fabs(dE) <= SCFACCURACY_3STEP_B && avgdP <= SCFACCURACY_3STEP_B) ||
			   (thisCutoff == cutoffC && fabs(dE) <= SCFACCURACY_3STEP_C && avgdP <= SCFACCURACY_3STEP_C)){

				// switch accuracy
				     if(thisCutoff == cutoffA) {thisCutoff = cutoffB; thisGrid = grid; }
				else if(thisCutoff == cutoffB) {thisCutoff = cutoffC;                  }
				else if(thisCutoff == cutoffC) {thisCutoff = realCutoff;               }

				// rebuilding fock matrix
				//for(i=0; i < nBasis; i++)
				//for(j=0; j < nBasis; j++){
				//	 GA[i*nBasis+j] = 0.0;
				//	 GB[i*nBasis+j] = 0.0;
				//	dPA[i*nBasis+j] = PA[i*nBasis+j];
				//	dPB[i*nBasis+j] = PB[i*nBasis+j];  
				//}
				printf("................. switch accuracy ..............\n");
				//printf("........ switch accuracy and rebuild fock matrix\n");
				fflush(stdout);

			}

		}

		// save checkpoint if requested
		if(opt->saveCheckAll){
			//printf("\n[Begin] saving checkpoint file to %s\n", opt->CheckFile);
			save_checkpoint(nBasis, gto, mol, Etot, dE, avgdP, CA, CB, eA, eB, opt);
			//printf("[Done]  saving checkpoint file to %s\n\n", opt->CheckFile);
			//fflush(stdout);
		}

	}while(notConverged);

	// report
	if(notConverged){
		printf("SCF have not converged because iter >= SCFMax\n");
	}else{
		printf("Done SCF Loop Total Energy is %20.8f Hartrees\n", Etot);

		// report DFT exchange energy at Hartree-Fock densities
		if(opt->method==METHOD_HF2xDFT){

			// print header
			printf(
			"                                                             \n"
			"-------------------------------------------------------------\n"
			"-----    DFT Exchange Energy @Hartree-Fock Density      -----\n"
			"-------------------------------------------------------------\n"
			"                                                             \n");

			// build exchange-correlation matrix
			memset(Exc, 0, grid->nPoint *sizeof(double));
			memset(XCa, 0, nBasis*nBasis*sizeof(double));
			memset(XCb, 0, nBasis*nBasis*sizeof(double));
			//getXCDirect(nBasis, nEA, nEB, CA, CB, PA, PB, thisGrid, gB, opt, gto, mol,
			//            thisCutoff,1,0,
			//            Exc, XCa, XCb,
			//            NULL, NULL, NULL);
			getXCDirect_LibXCcompat(nBasis, nEA, nEB, CA, CB, PA, PB, thisGrid, gB, opt, gto, mol,
			            thisCutoff,1,0,
			            Exc, XCa, XCb,
			            NULL, NULL, NULL);
	
			double DFTex=0.0;
			for(i=0; i < thisGrid->nPoint; i++)
				DFTex += Exc[i] * thisGrid->w[i];
			switch(opt->whichHF2xDFT){
				case HF2xDFT_xSLATER:  printf("  Slater Exchange = %10.6f Hartree\n",DFTex); break;
				case HF2xDFT_xPY86:    printf("    PY86 Exchange = %10.6f Hartree\n",DFTex); break;
				case HF2xDFT_xPBE:     printf("     PBE Exchange = %10.6f Hartree\n",DFTex); break;
				case HF2xDFT_xBECKE88: printf(" Becke88 Exchange = %10.6f Hartree\n",DFTex); break;
				case HF2xDFT_xMVS:     printf("     MVS Exchange = %10.6f Hartree\n",DFTex); break;
				case HF2xDFT_xCHACHIYO:printf("Chachiyo Exchange = %10.6f Hartree\n",DFTex); break;
			}

			double HFex=0.0;
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				HFex += (0.5*KA[i*nBasis+j])*PA[i*nBasis+j] +
				        (0.5*KB[i*nBasis+j])*PB[i*nBasis+j];
			}
			printf("      HF Exchange = %10.6f Hartree\n",HFex);
			printf("            Error = %10.6f\n",(DFTex-HFex));
		}
	}
	fflush(stdout);

	// re-calculate full eigen vector spectrum
	gen_sym_eigen(nBasis, FA, S, eA, CA);
	gen_sym_eigen(nBasis, FB, S, eB, CB);
	normalizeC(nBasis, nBasis, S, CA);
	normalizeC(nBasis, nBasis, S, CB);

	// save density matrix if requested
	if(opt->saveDMatrix){

		// openfile
		fd=fopen(opt->DMatrixFile,"w");
		if(fd==NULL){
			printf("uhf - error cannot open file %s for writing\n", 
			       opt->DMatrixFile);
		}

		///////////////////////////////////////////////////
		// loop and save each element, always write twice
		///////////////////////////////////////////////////
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			if((i*nBasis+j)%5==0) fprintf(fd,"\n");
			fprintf(fd," %15.8E ", PA[i*nBasis+j]);
		}
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			if((i*nBasis+j)%5==0) fprintf(fd,"\n");
			fprintf(fd," %15.8E ", PB[i*nBasis+j]);
		}
		
		// close file
		fclose(fd);

		// report status
		printf("Density matrix saved to %s\n", opt->DMatrixFile);
	}

	// save checkpoint if requested
	if(opt->saveCheck){
		printf("\n[Begin] saving checkpoint file to %s\n", opt->CheckFile);
		save_checkpoint(nBasis, gto, mol, Etot, dE, avgdP, CA, CB, eA, eB, opt);
		printf("[Done]  saving checkpoint file to %s\n\n", opt->CheckFile);
		fflush(stdout);
	}

	// clear convergence routine
	switch(opt->convMethod){
	case CONVMETHOD_DIIS4:   conv_diis4(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DIIS3:   conv_diis3(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DIIS2:   conv_diis2(0,   0.0, NULL, NULL); break;
	case CONVMETHOD_DAMPING: conv_damping(0, 0.0, NULL, NULL); break;
	}

	// clean memory
	free(Schwarz);
	free(JT);
	free(dJT);
	free(PA);
	free(KA);
	free(dPA);
	free(dKA);
	free(FA);
	free(PB);
	free(KB);
	free(dPB);
	free(dKB);
	free(FB);
	free(H);
	free(T);
	free(V);
	free(S);

	if(grid!=NULL) grid = cleanMolGrid(grid);
	if(gridSmall!=NULL) gridSmall = cleanMolGrid(gridSmall);
	//if(bg!=NULL) bg = cleanBGrid(bg);
	if(gB !=NULL) cleanGridBasisCenter(nBasis, gB);
	if(Exc!=NULL) free(Exc);
	if(XCa!=NULL) free(XCa);
	if(XCb!=NULL) free(XCb);

	// clean electron storage in JK integral subroutine
	//GTO_JK_Matrix_ShellSet(0, NULL, NULL, NULL, NULL, 0.0, NULL, NULL, NULL);
	//GTO_JK_Matrix_Quartet(0, NULL, NULL, NULL, NULL, 0.0, NULL, NULL, NULL);

	//
	// parallel version of cleaning 
	//
	int *status;             // status for each cpu
	int alldone;             // all idle flag

	// allocate memory
	status=calloc(opt->nCPU,sizeof(int));
	if(status==NULL){
		printf("uhf - error cannot allocate memory\n");
		exit(-1);
	}

	// reset status to idle
	for(i=(opt->nCPU-1);i>=0;i--) status[i] = RPC_IDLE;

	// clean all child processes
	do{
		// remote call to all cpu except childID=0
		for(i=(opt->nCPU-1);i>0;i--)
			status[i] = rpc_GTO_JK_Matrix_Quartet_Parallel(status[i],i,0,
			                                               NULL, NULL, NULL, NULL,
			                                               0.0, NULL, NULL, NULL, opt);
		// local call for childID=0
		if(status[0]==RPC_IDLE){
			GTO_JK_Matrix_Quartet_Parallel(0, 0, NULL, NULL, NULL, NULL,
			                               0.0, NULL, NULL, NULL, opt);
			status[0] = RPC_DONE;
		}

		// check if all done
		alldone=1;
		for(i=(opt->nCPU-1);i>=0;i--) if(status[i] != RPC_DONE) alldone=0;

	}while(!alldone);

	// clean memory
	free(status);

	if(notConverged) return 0.0;
	else             return Etot;
}

