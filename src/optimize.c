/*
	Siam Quantum [SQ] is a program pacakage that performs quantum
	mechanical calculations on atomic systems.

	Copyright (C) 2020
        Teepanis Chachiyo <teepanisc@nu.ac.th>
        Hathaithip Chachiyo <hathaithip.chachiyo@gmail.com>

      - Department of Physics, Faculty of Science,
        Naresuan University, Phitsanulok, Thailand

      - Thailand Center of Excellence in Physics,
        Ministry of Higher Education, Science, Research
        and Innovation, Thailand

	This file is a part of SQ.
	                                                       
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2, or (at your option) 
	any later version.                                                  
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "uhf.h"
#include "grad.h"
#include "optimize.h"
#include "util.h"

// invHessian_BFGS: compute inverse of Hessian matrix using BFBS schemme
// as explianed in the equation C.25b Szabo & Ostland "Modern Quantum Chemistry,
// Introduction to Advanced Electronic Structure Theory"
//
//                T       T
// G_n = B * G * B   + q*q 
//                     ---
//                     alpha
//
//                         T                      T
// where matrix B = I - q*d      and     alpha = q*d
//                      ---        
//                     alpha
//
// here, q is the steping vector dR and d is the change of gradient vector dGrad
//
// July 5,2010 - Theerapon Khamla
//     Initial implementation
//
// July 14, 2010 - Theerapon Khamla
//     Debuging
//
// Oct 21, 2010 - Teepanis Chachiyo & Theerapon Khamla
//     Change from updating Hessian to updating its inverse instead
//
void invHessian_BFGS(
	int nDim,            // matrix dimension
	const double *dR,    // displacement vector
	const double *dGrad, // change of gradient vector,
	double *invHessian){ // original inverse of the Hessian

	int i,j,p,q;       // loop index
	double alpha;      // dR DOT dGrad
	double *B;         // auxilary matrix
	double *G;         // updated value

	// memory allocation
	G=calloc(nDim*nDim, sizeof(double));
	B=calloc(nDim*nDim, sizeof(double));
	if(G==NULL || B==NULL){
		printf("invHessian_BFGS: error cannnot allocate memory\n");
		exit(-1);
	}

	// compute alpha
	alpha = 0.0;
	for(i=0; i < nDim; i++)
		alpha += dR[i] * dGrad[i];

	// check if alpha is zero
	if(alpha==0.0){
		printf("invHessian_BFGS: error dR perpendicular to dGrad\n");
		exit(-1);
	}

	// compute B matrix
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++)
		if(i==j) B[i*nDim+j] = 1.0 - dR[i] * dGrad[j] / alpha;
		else     B[i*nDim+j] =     - dR[i] * dGrad[j] / alpha;

	// update inverse of the Hessian
	for(p=0; p < nDim; p++)
	for(q=0; q < nDim; q++){
		G[p*nDim+q] = 0.0;

		for(i=0; i < nDim; i++)
		for(j=0; j < nDim; j++)
			G[p*nDim+q] +=   B[p*nDim+i]*B[q*nDim+j]*invHessian[i*nDim+j];

		G[p*nDim+q] += dR[p]*dR[q]/alpha;
	}

	// update value
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++)
		invHessian[i*nDim+j] = G[i*nDim+j];

	// free memory
	free(G);
	free(B);
}

// stepVector: calculate stepping vector using Newton method 
//         -1
// dR = - H  * Grad
// 
// Oct 21, 2010 - Teepanis Chachiyo & Theerapon Khamla
//    Take inverse of Hesssian instead of the actual Hessian itself as
//    an argument
//
// July 7,2010 - Theerapon Khamla
//    Initial implementation
//
void stepVector_Newton(
	int nDim,            // dimension of the matrix
	const double *invH,  // pointer to inverse of the Hessian matrix
	const double *grad,  // pointer to gradient vector
	double maxStepSize,  // maximum step size in each direction
	double *dR){         // returned step vector

	int i,j;    // loop index
	double r;   // step size

	// multiply inverse of the hessian to a negative of the gradient
	for(i=0; i < nDim; i++){
		dR[i] = 0.0;
		for(j=0; j < nDim; j++)
			dR[i] -= invH[i*nDim+j] * grad[j];
	}

	// validate step size
	if(maxStepSize <= 0.0){
		printf("stepVector_Newton: invalid range of maxStepSize");
		exit(-1);
	}

	// truncate the step size to maxStepSize
	r = 0;
	for(i=0; i < nDim; i++) r+= dR[i]*dR[i];
	r = sqrt(r);

	if(r > maxStepSize)
	for(i=0; i < nDim; i++)
		dR[i] = dR[i]*maxStepSize/r;
}


// deleteTranslation: delete simple translation from the geometry step
//
// Jan 1, 2016 - Teepanis Chachiyo
//               Initial implementation
//
void deleteTranslation(
	int nDim,     // vector dimension
	double *dR){ // displacement 

	double dx=0.0,dy=0.0,dz=0.0;
	int N;             // number of atom
	int i;             // loop index
	// validation
	if(nDim%3 != 0){
		printf("deleteTranslation: error nDim not multiple of 3\n");
		exit(-1);
	}
	
	N = nDim/3;
	// computer average displacement
	for(i=0; i < N; i++){
		dx += dR[i*3+0]/N;
		dy += dR[i*3+1]/N;
		dz += dR[i*3+2]/N;
	}

	// shift by this amount
	for(i=0; i < N; i++){
		dR[i*3+0] -= dx;
		dR[i*3+1] -= dy;
		dR[i*3+2] -= dz;
	}
}


//
// optimize: optmize molecular structure of the molecule using uhf module
//
// July 5, 2010 - Theerapon Khamla
//   Initial implementation
//
// Oct 22, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Added to Siam Quantum source code
//
// Dec 17, 2015 - Teepanis Chachiyo
//    Use previous molecular orbital to build next step guess density
//
// Jan 1, 2015 - Teepanis Chachiyo
//    delete translation from dR
//
void optimize(
	int dbSize,                         // number of record in basis set db
	const struct GTOBasisSet_t *basisDB,// pointer to basis set database
	struct Molecule_t *mol,             // returned molecular coordinate
	struct option_t *opt){              // options

	int nBasis;               // number of basis function
	struct GTOBasis_t *gto;   // pointer to basis function storage
	int nDim;                 // degree of freedoms
	double *invH;             // inverse of Hessian
	double *G;                // current gradient vector
	double *dR, *dG;          // steping vector and change of gradient
	double *fx,*fy,*fz;       // forces acting on nuclei
	double *CA, *eA, *CB, *eB;// molecular orbitals and eigen values
	double Etot;              // total energy
	int savedSCFGuess;        // global guess option

	int nIter;                // iteration index
	double force_max;         // maximum force on nuclei
	double force_rms;         // root mean square of force
	double step_max;          // maximum step size
	double step_rms;          // root mean square of step size

	int i,j;                  // loop index

	// save global guess option
	savedSCFGuess = opt->SCFGuess;

	// set degree of freedom
	nDim = mol->nAtom * 3;

	// memory allocation
	invH = calloc(nDim*nDim, sizeof(double));
	dR   = calloc(nDim,      sizeof(double));
	dG   = calloc(nDim,      sizeof(double));
	G    = calloc(nDim,      sizeof(double));
	fx   = calloc(mol->nAtom,sizeof(double));
	fy   = calloc(mol->nAtom,sizeof(double));
	fz   = calloc(mol->nAtom,sizeof(double));

	// initiale inverse of Hessian to identity matrix
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++)
		if(i==j) invH[i*nDim+j] = 1.0; else invH[i*nDim+j] = 0.0;

	// allocate memory for molecular orbitals and their eigen values
	gto   = genBasis(mol, &nBasis, dbSize, basisDB);
	cleanGTOBasis(gto,nBasis);
	CA = calloc(nBasis*nBasis,sizeof(double));
	eA = calloc(nBasis,sizeof(double));
	CB = calloc(nBasis*nBasis,sizeof(double));
	eB = calloc(nBasis,sizeof(double));
	if(CA==NULL || eA==NULL || CB==NULL || eB==NULL){
		printf("optimize: error - cannot allocate memory\n");
		exit(-1);
	}

	nIter=0;
	do{
		///////////////////////////////////////////////
		// perform scf calculation and compute forces
		///////////////////////////////////////////////
		printf(
		"                                                             \n"
		"                                                             \n"
		"-------------------------------------------------------------\n"
		"-----       GEOMETRY OPTIMIZATION Step %5d            -----\n"
		"-------------------------------------------------------------\n",
		nIter+1);
		fflush(stdout);

		// generate basis function
		gto   = genBasis(mol, &nBasis, dbSize, basisDB);

		//
		// scf calculation
		//

		// load orbital from previos cycle to form initial guess
		if(nIter>0) opt->SCFGuess = SCFGUESS_CACB;

		Etot = uhf(nBasis, gto, mol, 
		    get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		    CA, CB, eA, eB, opt);
		if(Etot==0.0){
			printf("optimize: error SCF calculation did not converge\n");
			exit(-1);
		}

		// compute force
		uhf_force(nBasis, gto, mol, 
		      get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		      CA, CB, eA, eB, opt, fx, fy, fz);

		////////////////////////////////////
		// update new molecular coordinate 
		////////////////////////////////////

		// construct gradient vector
		for(i=0; i < mol->nAtom; i++){
			// change of gradient vector
			dG[i*3+0] = -fx[i] - G[i*3+0];
			dG[i*3+1] = -fy[i] - G[i*3+1];
			dG[i*3+2] = -fz[i] - G[i*3+2];

			// current gradient vector
			G [i*3+0] = -fx[i];
			G [i*3+1] = -fy[i];
			G [i*3+2] = -fz[i];
		}

		// update inverse of Hessian if all info are available
		if(nIter > 0)
			invHessian_BFGS(nDim, dR, dG, invH);

		// compute steping vector
		stepVector_Newton(nDim, invH, G, MAXSTEPSIZE, dR);

		// delete translation
		//deleteTranslation(nDim, dR);

		// update molecular gemetry
		for(i=0; i < mol->nAtom; i++){
			mol->x[i] = mol->x[i] + dR[i*3+0];
			mol->y[i] = mol->y[i] + dR[i*3+1];
			mol->z[i] = mol->z[i] + dR[i*3+2];
		}

		// free memory
		cleanGTOBasis(gto,nBasis);

		// compute convergence criteria
		force_max = 0.0; force_rms = 0.0; step_max = 0.0; step_rms = 0.0;
		for(i=0; i < nDim; i++){
			if(fabs(G[i])  > force_max) force_max = fabs(G[i]);
			if(fabs(dR[i]) > step_max)  step_max  = fabs(dR[i]);
			force_rms += G[i]*G[i];
			step_rms  += dR[i]*dR[i];
		}
		force_rms = sqrt(force_rms/nDim);
		step_rms  = sqrt(step_rms/nDim);

		// report convergence status
		printf("Convergence Criterion    Value        Threshold\n");
		printf("  Maximum Force        %10.6f   %10.6f  ", force_max, OPT_CONV_FORCEMAX);
		if(force_max > OPT_CONV_FORCEMAX) printf("NO\n"); else printf("YES\n");
		printf("  RMS     Force        %10.6f   %10.6f  ", force_rms, OPT_CONV_FORCERMS);
		if(force_rms > OPT_CONV_FORCERMS) printf("NO\n"); else printf("YES\n");
		printf("  Maximum Displacement %10.6f   %10.6f  ", step_max, OPT_CONV_DISPMAX);
		if(step_max > OPT_CONV_DISPMAX) printf("NO\n"); else printf("YES\n");
		printf("  RMS     Displacement %10.6f   %10.6f  ", step_rms, OPT_CONV_DISPRMS);
		if(step_rms > OPT_CONV_DISPRMS) printf("NO\n"); else printf("YES\n");
		fflush(stdout);

		nIter++;

		if(nIter > opt->optMax){
			printf("Geometry optimization has not converged, nIter >= %d\n",opt->optMax);
			exit(-1);
		}

	}while(force_max > OPT_CONV_FORCEMAX || force_rms > OPT_CONV_FORCERMS ||
	       step_max  > OPT_CONV_DISPMAX  || step_rms  > OPT_CONV_DISPRMS);


	// print optimized geometry
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----                  OPTIMIZED GEOMETRY               -----\n"
	"-------------------------------------------------------------\n"
	"                                                             \n");
	printMolecule_XYZ(mol, stdout);

	// print distance matrix
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----                  DISTANCE MATRIX                  -----\n"
	"-------------------------------------------------------------\n"
	"                                                             \n"
	"                       Output Sequence                       \n"
/*
	"          +----------------------------------------          \n"
	"          |   Atom1   Atom2   Atom3   Atom4   ...            \n"
	"          +----------------------------------------          \n"
	"    Atom1 |    1                                             \n"
	"    Atom2 |    2       3                                     \n"
	"    Atom3 |    4       5       6                             \n"
	"    Atom4 |    7       8       9      10                     \n"
	"      :   |    .       .       .       .       .             \n"
	"          +----------------------------------------          \n"
	"                                                             \n"
	"                                                             \n"
*/

	"             +---------------------------------              \n"
	"             |   Atom1   Atom2   Atom3   ...                 \n"
	"             +---------------------------------              \n"
	"       Atom1 |    1                                          \n"
	"       Atom2 |    2       3                                  \n"
	"       Atom3 |    4       5       6                          \n"
	"         :   |    7       8       9      10                  \n"
	"                                                             \n"
	"                                                             \n"

/*
	"             +---------------------------------              \n"
	"             |   Atom1   Atom2   Atom3   Atom4               \n"
	"             +---------------------------------              \n"
	"       Atom1 |    1                                          \n"
	"       Atom2 |    2       3                                  \n"
	"       Atom3 |    4       5       6                          \n"
	"       Atom4 |    .       .       .       .                  \n"
	"                                                             \n"
*/
	"                      Output (Angstroms)                    \n"
	"                      ------------------                    \n");
	int nItem=0;
	for(i=0; i < mol->nAtom; i++)
	for(j=0; j <=i; j++){
		double dist=(mol->x[i]-mol->x[j])*(mol->x[i]-mol->x[j]) +
		            (mol->y[i]-mol->y[j])*(mol->y[i]-mol->y[j]) +
		            (mol->z[i]-mol->z[j])*(mol->z[i]-mol->z[j]);
		dist = sqrt(dist);
		printf("%11.5f ",dist/ANGSTROM2BOHR);
		nItem++;
		if(nItem%5==0) printf("\n");
	}
	// next line if not already done it
	if(nItem%5!=0) printf("\n");




	// free memory
	free(CA);
	free(eA);
	free(CB);
	free(eB);

	free(invH);
	free(dR);
	free(dG);
	free(G);
	free(fx);
	free(fy);
	free(fz);

	// set the guess option to original value
	opt->SCFGuess = savedSCFGuess;
}
