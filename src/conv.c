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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "conv.h"
#include "lin.h"
//
// conv_damping : update density matrix using a simple damping algorithm
// Pn = alpha*P_n + (1-alpha)P_n-1
// where P is an element of the density matrix
// 
// it returns average error between previous density and the current one
//
// call with nBasis=0, alpha=0.0, PA==NULL, PB==NULL to clear all memory
// and iteration information
//
// Sep 6, 2010 - Teepanis Chachiyo
//    Initial implementation and testing
//
// Oct 22, 2010 - Teepanis Chachiyo
//    Bugfix, forgot to set prevPA and prevPB to NULL upon reseting
//
double conv_damping(int nBasis,     // number of basis
                    double alpha,   // drag coefficient between 0 and 1
                    double *PA,     // alpha density matrix
                    double *PB){    // beta density matrix

	static double *prevPA=NULL;  // previous iteration alpha density matrix
	static double *prevPB=NULL;  // previous iteration beta density matrix
	static int nIter=0;          // number of iteration
	int i,j;                     // generic loop indexes
	double avgdP;                // average deviation

	// reset if requested
	if(nBasis==0 && alpha==0.0 && PA==NULL && PB==NULL){
		nIter=0;
		if(prevPA!=NULL){ free(prevPA); prevPA=NULL; }
		if(prevPB!=NULL){ free(prevPB); prevPB=NULL; }
		return 0.0;
	}

	// check range of alpha
	if(alpha < 0.0 || alpha > 1.0){
		printf("conv_damping : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix 
	////////////////////////////////////////////////////////
	if(nIter==0){

		// check if array is NULL
		if(prevPA!=NULL || prevPB!=NULL){
			printf("conv_damping : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		prevPA=calloc(nBasis*nBasis,sizeof(double));
		prevPB=calloc(nBasis*nBasis,sizeof(double));
		if(prevPA==NULL || prevPB==NULL){
			printf("conv_damping : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			prevPA[i*nBasis+j] = PA[i*nBasis+j];
			prevPB[i*nBasis+j] = PB[i*nBasis+j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	//////////////////////////////////////////////////////////////////
	/// for subsequent iteration we first compute average deviation
	//////////////////////////////////////////////////////////////////
	avgdP     = 0.0;
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		avgdP +=  (PA[i*nBasis+j]-prevPA[i*nBasis+j])
		         *(PA[i*nBasis+j]-prevPA[i*nBasis+j]);
		avgdP +=  (PB[i*nBasis+j]-prevPB[i*nBasis+j])
		         *(PB[i*nBasis+j]-prevPB[i*nBasis+j]);
	}
	avgdP = sqrt( avgdP/nBasis/nBasis/2.0);

	// update
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		PA[i*nBasis+j] = prevPA[i*nBasis+j] + 
		                 alpha *(PA[i*nBasis+j]-prevPA[i*nBasis+j]);

		PB[i*nBasis+j] = prevPB[i*nBasis+j] + 
		                 alpha *(PB[i*nBasis+j]-prevPB[i*nBasis+j]);
	}

	// store previous density matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		prevPA[i*nBasis+j] = PA[i*nBasis+j];
		prevPB[i*nBasis+j] = PB[i*nBasis+j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}


//
// conv_diis2 : is a simple 2 point diis scheme as described in "The 
// Maththematic of DIIS" by C. David Sherrill. In this algorithm we view
// the density matrix as a set of variables, call it P, which is iteratively
// produced in sequence. For example,
//
// P0(input) --> Solving Hartree-Fock --> P0(output)
// P1(input) --> Solving Hartree-Fock --> P1(output)
//   ...                                     ...
// Pi(input) --> Solving Hartree-Fock --> Pi(output)
//
// When we say the calculation has converged during the SCF cycle, we mean
// that the (input) and (output) are the same. In other words, 
//
//                 Pi(output) - Pi(input) = dPi = 0
//
// The whole question behind DIIS algorithm is this:
//
//           "How to find Pi(input) such that its dPi = 0 ?"
//
// To better quantify the residual dPi, which we want it to vanish, we define
// its norm as 
//                          Bii = dPi * dPi,      ---------------- (1)
//
// which is very similar to dot product of two vectors because Pi is simply a 
// collection of variables.
//
// DIIS assumes that we can write Pi(input) as a linear combination of the
// previous value, For example
//
//           P4(input) = c3*P3(output) + c2*P2(output) + c1*P1(output)
//
// Here, {c3,c2,c1) is a set of coefficients which we also assume that they
// sum to exaxtly 1. In mathetical form, (c3+c2+c1) = 1.
//
// In this implementation, we choose to go as far as 2 points back from the
// current iteration, so that 
//
//           Pi(input) = c1*Pi-1(input) + c0*Pi-2(input) ----------- (2)
//
// The question is how to find a suitable {c2,c1} such that the norm
// in Eq{1} is minimized subject to the contraint 
//
//             c1 + c0 = 1       ----------------------------------- (3)
//
// To answer such question we use a fairly standard method of minimization
// with constraint called "Lagrange Multiplier Method" which will not be
// explained here (See Wikipedia on the subject).
//
// As it turns out, we can find the {c1,c0} that satisfies the above conditions 
// by solving a linear matrix equation like this
//
//              -           -  -        -      -     - 
//              | B00 B01 -1|  |   c0   |      |  0  |
//              | B10 B11 -1|  |   c1   |   =  |  0  |
//              | -1  -1   0|  | lambda |      | -1  |
//              -           -  -        -      -     -
//
// Here Bij is a dot product dPi*dPj. We can see that c1 and c2 is very
// simple to calculate for a 2-point DIIS. They are
//
//   c0 = (B11-B01)/(B00+B11-2.0*B01) and  c1 = (B00-B01)/(B00+B11-2.0*B01)
//                                     ------------------------------ (4)
//
// In general, the DIIS algorithm is very good at accelerating the convergence.
// If the system is already converging, DIIS makes it faster. But the DIIS
// does not help if it is already not converging at all. But, how to deal
// with bad convergence behavior is another matter for another time.
//
// Sep 30, 2012 - Teepanis Chachiyo
//      1) Check if DIIS fails using criteria B11>B00 and use simple damping
//      in this case
//      2) Rewrite the formula for computing c0 and c1 to avoid dividing by
//      zero.
// 
// May 21, 2011 - Teepanis Chachiyo
//      Initial implementation and testing
//
double conv_diis2(int nBasis,      // number of basis
                  double alpha,    // drag coefficient between 0 and 1
                  double *PA,      // output alpha density matrix
                  double *PB){     // output beta  density matrix

	static double *PA0=NULL;     // input alpha density matrix
	static double *PB0=NULL;     // input beta  density matrix
	static double *PA1=NULL;     // input alpha density matrix
	static double *PB1=NULL;     // input beta  density matrix
	static double *dPA0=NULL;    // error alpha density matrix
	static double *dPB0=NULL;    // error beta  density matrix
	static double *dPA1=NULL;    // error alpha density matrix
	static double *dPB1=NULL;    // error beta  density matrix
	static int nIter=0;          // number of iteration
	int i,j;                     // generic loop indexes
	double avgdP;                // average deviation
	double c1,c0,B00,B01,B11;    // diis variables

	// reset if requested
	if(nBasis==0 && alpha==0.0 && PA==NULL && PB==NULL){
		nIter=0;
		if(PA0!=NULL){ free(PA0); PA0=NULL; }
		if(PB0!=NULL){ free(PB0); PB0=NULL; }
		if(PA1!=NULL){ free(PA1); PA1=NULL; }
		if(PB1!=NULL){ free(PB1); PB1=NULL; }

		if(dPA0!=NULL){ free(dPA0); dPA0=NULL; }
		if(dPB0!=NULL){ free(dPB0); dPB0=NULL; }
		if(dPA1!=NULL){ free(dPA1); dPA1=NULL; }
		if(dPB1!=NULL){ free(dPB1); dPB1=NULL; }
		return 0.0;
	}

	// check range of alpha
	if(alpha < 0.0 || alpha > 1.0){
		printf("conv_diis2 : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix 
	////////////////////////////////////////////////////////
	if(nIter==0){

		// check if array is NULL
		if( PA0!=NULL ||  PB0!=NULL ||
		    PA1!=NULL ||  PB1!=NULL ||
		   dPA0!=NULL || dPB0!=NULL || 
		   dPA1!=NULL || dPB1!=NULL){
			printf("conv_diis2 : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		PA0 =calloc(nBasis*nBasis,sizeof(double));
		PB0 =calloc(nBasis*nBasis,sizeof(double));
		PA1 =calloc(nBasis*nBasis,sizeof(double));
		PB1 =calloc(nBasis*nBasis,sizeof(double));
		dPA0=calloc(nBasis*nBasis,sizeof(double));
		dPB0=calloc(nBasis*nBasis,sizeof(double));
		dPA1=calloc(nBasis*nBasis,sizeof(double));
		dPB1=calloc(nBasis*nBasis,sizeof(double));
		if( PA0==NULL ||  PB0==NULL ||
		    PA1==NULL ||  PB1==NULL ||
		   dPA0==NULL || dPB0==NULL ||
		   dPA1==NULL || dPB1==NULL){
			printf("conv_diis2 : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			PA1[i*nBasis+j] = PA[i*nBasis+j];
			PB1[i*nBasis+j] = PB[i*nBasis+j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	////////////////////////////////////////////////////////
	/// subsequent iteration
	////////////////////////////////////////////////////////

	// compute rms deviation
	avgdP     = 0.0;
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		avgdP +=  (PA[i*nBasis+j]-PA1[i*nBasis+j])
		         *(PA[i*nBasis+j]-PA1[i*nBasis+j]);
		avgdP +=  (PB[i*nBasis+j]-PB1[i*nBasis+j])
		         *(PB[i*nBasis+j]-PB1[i*nBasis+j]);
	}
	avgdP = sqrt( avgdP/nBasis/nBasis/2.0);

	// compute dPA1 and dPB1
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
			dPA1[i*nBasis+j] = PA[i*nBasis+j] - PA1[i*nBasis+j];
			dPB1[i*nBasis+j] = PB[i*nBasis+j] - PB1[i*nBasis+j];
	}

	switch(nIter){

	case 1: // update using simple damping for the first iteration
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			PA[i*nBasis+j]  = alpha * PA[i*nBasis+j] +
			                 (1.0-alpha) * PA1[i*nBasis+j];

			PB[i*nBasis+j]  = alpha * PB[i*nBasis+j] +
			                 (1.0-alpha) * PB1[i*nBasis+j];
		}
	break;

	default: // using 2-point DIIS for the rest

		// compute diis variables
		B00 = B11 = B01 = 0.0;

		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
#define DIIS_B(X,Y) B##X##Y += dPA##X[i*nBasis+j]*dPA##Y[i*nBasis+j]; \
                    B##X##Y += dPB##X[i*nBasis+j]*dPB##Y[i*nBasis+j];

			DIIS_B(0,0); DIIS_B(1,1); DIIS_B(0,1);
		}
		c1 = (B00-B01)/(B11-B01);
		c0 = 1.0/(1.0+c1);
		c1 = c1*c0;

		// use simple damping if DIIS fails
		if(B11>B00){
			c1 = alpha;
			c0 = 1.0-alpha;
		}

		// apply DIIS updates
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			PA[i*nBasis+j]  = c1*(PA[i*nBasis+j]) + 
			                  c0*(PA0[i*nBasis+j]+dPA0[i*nBasis+j]);
			PB[i*nBasis+j]  = c1*(PB[i*nBasis+j]) + 
			                  c0*(PB0[i*nBasis+j]+dPB0[i*nBasis+j]);
		}
	break;

	}

	// store previous density matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		dPA0[i*nBasis+j] = dPA1[i*nBasis+j];
		dPB0[i*nBasis+j] = dPB1[i*nBasis+j];

		PA0[i*nBasis+j] = PA1[i*nBasis+j];
		PB0[i*nBasis+j] = PB1[i*nBasis+j];

		PA1[i*nBasis+j] = PA[i*nBasis+j];
		PB1[i*nBasis+j] = PB[i*nBasis+j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}


//
// conv_diis3 : is a the 3-point version of the conv_diis2
//
// Oct, 2012 - Teepanis Chachiyo
//      Initial implementation and testing
//
double conv_diis3(int nBasis,      // number of basis
                  double alpha,    // drag coefficient between 0 and 1
                  double *PA,      // output alpha density matrix
                  double *PB){     // output beta  density matrix

	static double *PA0=NULL;        // input alpha density matrix
	static double *PB0=NULL;        // input beta  density matrix
	static double *PA1=NULL;        // input alpha density matrix
	static double *PB1=NULL;        // input beta  density matrix
	static double *PA2=NULL;        // input alpha density matrix
	static double *PB2=NULL;        // input beta  density matrix
	static double *dPA0=NULL;       // error alpha density matrix
	static double *dPB0=NULL;       // error beta  density matrix
	static double *dPA1=NULL;       // error alpha density matrix
	static double *dPB1=NULL;       // error beta  density matrix
	static double *dPA2=NULL;       // error alpha density matrix
	static double *dPB2=NULL;       // error beta  density matrix
	static int nIter=0;             // number of iteration
	int i,j;                        // generic loop indexes
	double avgdP;                   // average deviation
	double c2,c1,c0;                // diis coefficients
	double B00,B11,B22,B01,B02,B12; // diis variables

	// reset if requested
	if(nBasis==0 && alpha==0.0 && PA==NULL && PB==NULL){
		nIter=0;
		if(PA0!=NULL){ free(PA0); PA0=NULL; }
		if(PB0!=NULL){ free(PB0); PB0=NULL; }
		if(PA1!=NULL){ free(PA1); PA1=NULL; }
		if(PB1!=NULL){ free(PB1); PB1=NULL; }
		if(PA2!=NULL){ free(PA2); PA2=NULL; }
		if(PB2!=NULL){ free(PB2); PB2=NULL; }

		if(dPA0!=NULL){ free(dPA0); dPA0=NULL; }
		if(dPB0!=NULL){ free(dPB0); dPB0=NULL; }
		if(dPA1!=NULL){ free(dPA1); dPA1=NULL; }
		if(dPB1!=NULL){ free(dPB1); dPB1=NULL; }
		if(dPA2!=NULL){ free(dPA2); dPA2=NULL; }
		if(dPB2!=NULL){ free(dPB2); dPB2=NULL; }
		return 0.0;
	}

	// check range of alpha
	if(alpha < 0.0 || alpha > 1.0){
		printf("conv_diis3 : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix 
	////////////////////////////////////////////////////////
	if(nIter==0){

		// check if array is NULL
		if( PA0!=NULL ||  PB0!=NULL ||
		    PA1!=NULL ||  PB1!=NULL ||
		    PA2!=NULL ||  PB2!=NULL ||
		   dPA0!=NULL || dPB0!=NULL ||
		   dPA1!=NULL || dPB1!=NULL ||  
		   dPA2!=NULL || dPB2!=NULL){
			printf("conv_diis3 : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		PA0 =calloc(nBasis*nBasis,sizeof(double));
		PB0 =calloc(nBasis*nBasis,sizeof(double));
		PA1 =calloc(nBasis*nBasis,sizeof(double));
		PB1 =calloc(nBasis*nBasis,sizeof(double));
		PA2 =calloc(nBasis*nBasis,sizeof(double));
		PB2 =calloc(nBasis*nBasis,sizeof(double));
		dPA0=calloc(nBasis*nBasis,sizeof(double));
		dPB0=calloc(nBasis*nBasis,sizeof(double));
		dPA1=calloc(nBasis*nBasis,sizeof(double));
		dPB1=calloc(nBasis*nBasis,sizeof(double));
		dPA2=calloc(nBasis*nBasis,sizeof(double));
		dPB2=calloc(nBasis*nBasis,sizeof(double));
		if( PA0==NULL ||  PB0==NULL ||
		    PA1==NULL ||  PB1==NULL ||
		    PA2==NULL ||  PB2==NULL ||
		   dPA0==NULL || dPB0==NULL ||
		   dPA1==NULL || dPB1==NULL ||
		   dPA2==NULL || dPB2==NULL){
			printf("conv_diis3 : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for(i=0; i < nBasis; i++) 
		for(j=0; j < nBasis; j++){
			PA2[i*nBasis+j] = PA[i*nBasis+j];
			PB2[i*nBasis+j] = PB[i*nBasis+j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	////////////////////////////////////////////////////////
	/// subsequent iteration
	////////////////////////////////////////////////////////

	// compute rms deviation
	avgdP     = 0.0;
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		avgdP +=  (PA[i*nBasis+j]-PA2[i*nBasis+j])
		         *(PA[i*nBasis+j]-PA2[i*nBasis+j]);
		avgdP +=  (PB[i*nBasis+j]-PB2[i*nBasis+j])
		         *(PB[i*nBasis+j]-PB2[i*nBasis+j]);
	}
	avgdP = sqrt( avgdP/nBasis/nBasis/2.0);

	// compute dPA2 and dPB2
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
			dPA2[i*nBasis+j] = PA[i*nBasis+j] - PA2[i*nBasis+j];
			dPB2[i*nBasis+j] = PB[i*nBasis+j] - PB2[i*nBasis+j];
	}

	switch(nIter){

	case 1: // update using simple damping for the first iteration
	case 2: // update using simple damping for the second iteration
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			PA[i*nBasis+j]  = alpha * PA[i*nBasis+j] +
			                 (1.0-alpha) * PA2[i*nBasis+j];

			PB[i*nBasis+j]  = alpha * PB[i*nBasis+j] +
			                 (1.0-alpha) * PB2[i*nBasis+j];
		}
	break;

	default: // using 3-point DIIS for the rest

		// compute diis variables
		B00 = B11 = B22 = B01 = B02 = B12 = 0.0;

		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			DIIS_B(0,0); DIIS_B(1,1); DIIS_B(2,2);
			DIIS_B(0,1); DIIS_B(0,2); DIIS_B(1,2);
		}
		c1 = (B02*B02+B00*B12-B01*B02-B00*B22+B01*B22-B02*B12)/
		     (B12*B12-B01*B12+B02*B11+B01*B22-B02*B12-B11*B22); 
		c2 = (B01*B01-B00*B11-B01*B02+B02*B11+B12*B00-B12*B01)/
		     (B12*B12-B01*B12+B02*B11+B01*B22-B02*B12-B11*B22); 
		c0 = 1.0/(1.0+c1+c2);
		c1 = c1*c0;
		c2 = c2*c0;

		// use simple damping if DIIS fails
		if(B22>B11){
			c2 = alpha;
			c1 = 1.0-alpha;
			c0 = 0.0;
		}
		
		// apply DIIS updates
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			PA[i*nBasis+j]  = c2*(PA[i*nBasis+j]) + 
			                  c1*(PA1[i*nBasis+j]+dPA1[i*nBasis+j]) +
			                  c0*(PA0[i*nBasis+j]+dPA0[i*nBasis+j]);
			PB[i*nBasis+j]  = c2*(PB[i*nBasis+j]) + 
			                  c1*(PB1[i*nBasis+j]+dPB1[i*nBasis+j]) +
			                  c0*(PB0[i*nBasis+j]+dPB0[i*nBasis+j]);
		}

	break;

	}

	// store previous density matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		dPA0[i*nBasis+j] = dPA1[i*nBasis+j];
		dPB0[i*nBasis+j] = dPB1[i*nBasis+j];

		dPA1[i*nBasis+j] = dPA2[i*nBasis+j];
		dPB1[i*nBasis+j] = dPB2[i*nBasis+j];

		PA0[i*nBasis+j] = PA1[i*nBasis+j];
		PB0[i*nBasis+j] = PB1[i*nBasis+j];

		PA1[i*nBasis+j] = PA2[i*nBasis+j];
		PB1[i*nBasis+j] = PB2[i*nBasis+j];

		PA2[i*nBasis+j] = PA[i*nBasis+j];
		PB2[i*nBasis+j] = PB[i*nBasis+j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}

//
// conv_diis4 : is a the 4-point version of the conv_diis2
//
// Jan, 2013 - Teepanis Chachiyo
//      Initial implementation and testing
//
double conv_diis4(int nBasis,      // number of basis
                  double alpha,    // drag coefficient between 0 and 1
                  double *PA,      // output alpha density matrix
                  double *PB){     // output beta  density matrix

	static double *PA0=NULL;        // input alpha density matrix
	static double *PB0=NULL;        // input beta  density matrix
	static double *PA1=NULL;        // input alpha density matrix
	static double *PB1=NULL;        // input beta  density matrix
	static double *PA2=NULL;        // input alpha density matrix
	static double *PB2=NULL;        // input beta  density matrix
	static double *PA3=NULL;        // input alpha density matrix
	static double *PB3=NULL;        // input beta  density matrix
	static double *dPA0=NULL;       // error alpha density matrix
	static double *dPB0=NULL;       // error beta  density matrix
	static double *dPA1=NULL;       // error alpha density matrix
	static double *dPB1=NULL;       // error beta  density matrix
	static double *dPA2=NULL;       // error alpha density matrix
	static double *dPB2=NULL;       // error beta  density matrix
	static double *dPA3=NULL;       // error alpha density matrix
	static double *dPB3=NULL;       // error beta  density matrix
	static int nIter=0;             // number of iteration
	int i,j;                        // generic loop indexes
	double avgdP;                   // average deviation
	double c3,c2,c1,c0;             // diis coefficients
	double B00,B01,B02,B03,B11,
	       B12,B13,B22,B23,B33;     // diis variables

	// reset if requested
	if(nBasis==0 && alpha==0.0 && PA==NULL && PB==NULL){
		nIter=0;
		if(PA0!=NULL){ free(PA0); PA0=NULL; }
		if(PB0!=NULL){ free(PB0); PB0=NULL; }
		if(PA1!=NULL){ free(PA1); PA1=NULL; }
		if(PB1!=NULL){ free(PB1); PB1=NULL; }
		if(PA2!=NULL){ free(PA2); PA2=NULL; }
		if(PB2!=NULL){ free(PB2); PB2=NULL; }
		if(PA3!=NULL){ free(PA3); PA3=NULL; }
		if(PB3!=NULL){ free(PB3); PB3=NULL; }

		if(dPA0!=NULL){ free(dPA0); dPA0=NULL; }
		if(dPB0!=NULL){ free(dPB0); dPB0=NULL; }
		if(dPA1!=NULL){ free(dPA1); dPA1=NULL; }
		if(dPB1!=NULL){ free(dPB1); dPB1=NULL; }
		if(dPA2!=NULL){ free(dPA2); dPA2=NULL; }
		if(dPB2!=NULL){ free(dPB2); dPB2=NULL; }
		if(dPA3!=NULL){ free(dPA3); dPA3=NULL; }
		if(dPB3!=NULL){ free(dPB3); dPB3=NULL; }
		return 0.0;
	}

	// check range of alpha
	if(alpha < 0.0 || alpha > 1.0){
		printf("conv_diis3 : invalid range of alpha\n");
		exit(-1);
	}

	////////////////////////////////////////////////////////
	/// at the zeroth iteration we simply allocate memory
	/// and store previous density matrix 
	////////////////////////////////////////////////////////
	if(nIter==0){

		// check if array is NULL
		if( PA0!=NULL ||  PB0!=NULL ||
		    PA1!=NULL ||  PB1!=NULL ||
		    PA2!=NULL ||  PB2!=NULL ||
		    PA3!=NULL ||  PB3!=NULL ||
		   dPA0!=NULL || dPB0!=NULL ||
		   dPA1!=NULL || dPB1!=NULL ||  
		   dPA2!=NULL || dPB2!=NULL ||
		   dPA3!=NULL || dPB3!=NULL ){
			printf("conv_diis4 : detect prior initialized pointers\n");
			exit(-1);
		}

		// allocate memory
		PA0 =calloc(nBasis*nBasis,sizeof(double));
		PB0 =calloc(nBasis*nBasis,sizeof(double));
		PA1 =calloc(nBasis*nBasis,sizeof(double));
		PB1 =calloc(nBasis*nBasis,sizeof(double));
		PA2 =calloc(nBasis*nBasis,sizeof(double));
		PB2 =calloc(nBasis*nBasis,sizeof(double));
		PA3 =calloc(nBasis*nBasis,sizeof(double));
		PB3 =calloc(nBasis*nBasis,sizeof(double));
		dPA0=calloc(nBasis*nBasis,sizeof(double));
		dPB0=calloc(nBasis*nBasis,sizeof(double));
		dPA1=calloc(nBasis*nBasis,sizeof(double));
		dPB1=calloc(nBasis*nBasis,sizeof(double));
		dPA2=calloc(nBasis*nBasis,sizeof(double));
		dPB2=calloc(nBasis*nBasis,sizeof(double));
		dPA3=calloc(nBasis*nBasis,sizeof(double));
		dPB3=calloc(nBasis*nBasis,sizeof(double));
		if( PA0==NULL ||  PB0==NULL ||
		    PA1==NULL ||  PB1==NULL ||
		    PA2==NULL ||  PB2==NULL ||
		    PA3==NULL ||  PB3==NULL ||
		   dPA0==NULL || dPB0==NULL ||
		   dPA1==NULL || dPB1==NULL ||
		   dPA2==NULL || dPB2==NULL ||
		   dPA3==NULL || dPB3==NULL){
			printf("conv_diis4 : cannot allocate memory\n");
			exit(-1);
		}

		// store previous density matrix
		for(i=0; i < nBasis; i++) 
		for(j=0; j < nBasis; j++){
			PA3[i*nBasis+j] = PA[i*nBasis+j];
			PB3[i*nBasis+j] = PB[i*nBasis+j];
		}

		// increase interation index
		nIter++;

		// simply return 0
		return 0.0;
	}

	////////////////////////////////////////////////////////
	/// subsequent iteration
	////////////////////////////////////////////////////////

	// compute rms deviation
	avgdP     = 0.0;
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		avgdP +=  (PA[i*nBasis+j]-PA3[i*nBasis+j])
		         *(PA[i*nBasis+j]-PA3[i*nBasis+j]);
		avgdP +=  (PB[i*nBasis+j]-PB3[i*nBasis+j])
		         *(PB[i*nBasis+j]-PB3[i*nBasis+j]);
	}
	avgdP = sqrt( avgdP/nBasis/nBasis/2.0);

	// compute dPA3 and dPB3
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
			dPA3[i*nBasis+j] = PA[i*nBasis+j] - PA3[i*nBasis+j];
			dPB3[i*nBasis+j] = PB[i*nBasis+j] - PB3[i*nBasis+j];
	}

	switch(nIter){

	case 1: // update using simple damping for the first iteration
	case 2: // update using simple damping for the second iteration
	case 3: // update using simple damping for the second iteration
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			PA[i*nBasis+j]  = alpha * PA[i*nBasis+j] +
			                 (1.0-alpha) * PA3[i*nBasis+j];

			PB[i*nBasis+j]  = alpha * PB[i*nBasis+j] +
			                 (1.0-alpha) * PB3[i*nBasis+j];
		}
	break;

	default: // using 4-point DIIS for the rest

		// compute diis variables
		B00 = B01 = B02 = B03 = B11 = B12 = B13 = B22 = B23 = B33 = 0.0;

		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			DIIS_B(0,0); DIIS_B(0,1); DIIS_B(0,2); DIIS_B(0,3);
			             DIIS_B(1,1); DIIS_B(1,2); DIIS_B(1,3);
			                          DIIS_B(2,2); DIIS_B(2,3);
			                                       DIIS_B(3,3);
		}

		double A[25] = {  0.0,  0.0,  0.0,  0.0, -1.0,
		                  0.0,  0.0,  0.0,  0.0, -1.0,
		                  0.0,  0.0,  0.0,  0.0, -1.0,
		                  0.0,  0.0,  0.0,  0.0, -1.0,
		                 -1.0, -1.0, -1.0, -1.0,  0.0};
		double B[5]  = {  0.0,  0.0,  0.0,  0.0, -1.0};
		A[ 0] = B00; A[ 1] = B01; A[ 2] = B02; A[ 3] = B03;
		A[ 5] = B01; A[ 6] = B11; A[ 7] = B12; A[ 8] = B13;
		A[10] = B02; A[11] = B12; A[12] = B22; A[13] = B23;
		A[15] = B03; A[16] = B13; A[17] = B23; A[18] = B33;
		linear_solver(5,1,A,B);
		c0 = B[0]; c1 = B[1]; c2 = B[2]; c3 = B[3];


		// use simple damping if DIIS fails
		if(B33>B22){
			c3 = alpha;
			c2 = 1.0-alpha;
			c1 = 0.0;
			c0 = 0.0;
		}
		
		// apply DIIS updates
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			PA[i*nBasis+j]  = c3*(PA[i*nBasis+j]) + 
			                  c2*(PA2[i*nBasis+j]+dPA2[i*nBasis+j]) +
			                  c1*(PA1[i*nBasis+j]+dPA1[i*nBasis+j]) +
			                  c0*(PA0[i*nBasis+j]+dPA0[i*nBasis+j]);
			PB[i*nBasis+j]  = c3*(PB[i*nBasis+j]) + 
			                  c2*(PB2[i*nBasis+j]+dPB2[i*nBasis+j]) +
			                  c1*(PB1[i*nBasis+j]+dPB1[i*nBasis+j]) +
			                  c0*(PB0[i*nBasis+j]+dPB0[i*nBasis+j]);
		}

	break;

	}

	// store previous density matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		dPA0[i*nBasis+j] = dPA1[i*nBasis+j];
		dPB0[i*nBasis+j] = dPB1[i*nBasis+j];

		dPA1[i*nBasis+j] = dPA2[i*nBasis+j];
		dPB1[i*nBasis+j] = dPB2[i*nBasis+j];

		dPA2[i*nBasis+j] = dPA3[i*nBasis+j];
		dPB2[i*nBasis+j] = dPB3[i*nBasis+j];

		PA0[i*nBasis+j] = PA1[i*nBasis+j];
		PB0[i*nBasis+j] = PB1[i*nBasis+j];

		PA1[i*nBasis+j] = PA2[i*nBasis+j];
		PB1[i*nBasis+j] = PB2[i*nBasis+j];

		PA2[i*nBasis+j] = PA3[i*nBasis+j];
		PB2[i*nBasis+j] = PB3[i*nBasis+j];

		PA3[i*nBasis+j] = PA[i*nBasis+j];
		PB3[i*nBasis+j] = PB[i*nBasis+j];
	}

	// increment iteration number
	nIter++;

	// return average deviation
	return avgdP;
}
