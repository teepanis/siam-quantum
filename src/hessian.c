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
#include "hessian.h"
#include "util.h"
#include "lin.h"
#include "qmd.h"

struct atomic_mass_t{
	int Z;
	int isotope;
	double amu_mass;
};

double Z2_amu_mass(int Z, int *isotope){
int maxZ = 18;
struct atomic_mass_t data[] ={
{  0,  0, 0.0},               // none
{  1,  1, 1.0078250321},      // h
{  2,  4, 4.0026032497},      // he
{  3,  7, 7.0160040},         // li
{  4,  9, 9.0121821},         // be
{  5, 11, 11.00930555},       // b
{  6, 12, 12.00000000},       // c
{  7, 14, 14.00307400529},    // n
{  8, 16, 15.994914622115},   // o
{  9, 19, 18.998403207},      // f
{ 10, 20, 19.992440175920},   // ne

{ 11,  23, 22.9897696723},    // na
{ 12,  24, 23.9850419020},    // mg
{ 13,  27, 26.9815384414},    // al
{ 14,  28, 27.976926532720},  // si
{ 15,  31, 30.9737615120},    // p
{ 16,  32, 31.9720706912},    // s
{ 17,  35, 34.968852714},     // cl
{ 18,  40, 39.9623831233}     // ar

};

	if(Z > maxZ){
		printf("Z2MASS : error Z too large\n");
		exit(-1);
	}
	if(Z != data[Z].Z){
		printf("Z2MASS : error inconsistent data\n");
		exit(-1);
	}

	double m = data[Z].amu_mass;

	if(isotope != NULL) *isotope = data[Z].isotope;
	return m;
}


//
// hessian_output: print hessian matrix and compute subsequent
//                 quantities like vibration frequencies, 
//                 or zero point energy
//
// March 26, 2020 - Teepanis Chachiyo
//   Initial implementation and testing
//
void hessian_output(
	int nDim,
	const double *H,
	const struct Molecule_t *mol){

	double *mwH;    // mass weighted hessian
	double *S;      // identity matrix
	double *C;      // normal mode eigen vector
	double *lambda; // eigen values

	int i,j,nItem;        // index

	// allocate arrays
	mwH    = calloc(nDim*nDim,sizeof(double));
	S      = calloc(nDim*nDim,sizeof(double));
	C      = calloc(nDim*nDim,sizeof(double));
	lambda = calloc(nDim,sizeof(double));
	if(mwH==NULL || S==NULL || C==NULL || lambda==NULL){
		printf("hessian_output: error - cannot allocate mwH,S,C,lambda\n");
		exit(-1);
	}
	for(i=0; i < nDim; i++) S[i*nDim+i] = 1.0;

	// print hessian matrix
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------- HESSIAN NUMERICAL STEPS COMPLETED ------------\n"
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----                  HESSIAN MATRIX                   -----\n"
	"-------------------------------------------------------------\n"
	"                                                             \n"
	"                       Output Sequence                       \n"
	"             +--------------------------------               \n"
	"             |   X1    Y1    Z1    X2    ...                 \n"
	"             +--------------------------------               \n"
	"         X1  |    1                                          \n"
	"         Y1  |    2     3                                    \n"
	"         Z1  |    4     5     6                              \n"
	"         X2  |    7     8     9    10                        \n"
	"         :   |    .     .     .     .     .                  \n"
	"             +--------------------------------               \n"
//	"         :   |    .     .     .     .     .                  \n"
//	"         :   |   11    12    13    14    15                  \n"
	"                                                             \n"
	"                                                             \n"
	"                     Output (Hartree/Bohr^2)                 \n"
	"                     -----------------------                 \n");
	nItem=0;
	for(i=0; i < nDim; i++)
	for(j=0; j <=i; j++){
		printf("%11.5f ",H[i*nDim+j]);
		nItem++;
		if(nItem%5==0) printf("\n");
	}
	// next line if not already done it
	if(nItem%5!=0) printf("\n");

	////////////////////////////
	// Molecular Vibration
	////////////////////////////
	printf(
	"                                                             \n"
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----                 MOLECULAR VIBRATION               -----\n"
	"-------------------------------------------------------------\n");

	// transform to mass-weighted hessian in SI
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++){
		mwH[i*nDim+j] = H[i*nDim+j]*HARTREE2J
		                           /BOHR2ANGSTROM/1.0E-10
		                           /BOHR2ANGSTROM/1.0E-10
		                /sqrt( Z2_amu_mass(mol->Z[i/3], NULL)*AMU2KG )
		                /sqrt( Z2_amu_mass(mol->Z[j/3], NULL)*AMU2KG );
	}

	// compute normal modes and frequencies
	gen_sym_eigen(nDim, mwH, S, lambda, C);

	// compute mass-weighted eigen vector
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++)
		C[i*nDim+j] = C[i*nDim+j]/sqrt( Z2_amu_mass(mol->Z[j/3], NULL));

	// convert lambda to frequencies in cm^-1
	for(i=0; i < nDim; i++){
		if(lambda[i] >= 0.0) lambda[i] =  sqrt( lambda[i]);
		else                 lambda[i] = -sqrt(-lambda[i]);

		lambda[i] = lambda[i]*1.0E-6*MHZ2CM/2.0/M_PI;

	}

	// starting mode, exculding translation and rotation
	int iStart=0;
	     if(mol->nAtom == 1) iStart=3;
	else if(mol->nAtom == 2) iStart=5;
	else                     iStart=6;

	printf(
	"                                                             \n"
	"                                                             \n"
	"                   Normal Mode Output Format                 \n"
	"                   -------------------------                 \n"
	"                   ID   Freq   X1  Y1  Z1                    \n"
	"                               X2  Y2  Z2                    \n"
	"                               X3  Y3  Z3                    \n"
	"                               :   :   :                     \n"
	"                                                             \n"
	"                            Output                          \n"
	"                            ------                          \n");

	// print normal modes
	for(i=iStart; i < nDim; i++){
		printf("%5d %11.1lf %12.4lf %12.4lf %12.4lf\n",
		       i-iStart+1, lambda[i], C[i*nDim+0], C[i*nDim+1], C[i*nDim+2]);
		for(j=1; j < mol->nAtom; j++)
		printf("                  %12.4lf %12.4lf %12.4lf\n",
		       C[i*nDim+j*3+0],C[i*nDim+j*3+1],C[i*nDim+j*3+2]);
	}


	// print
	printf(
	"                                                             \n"
	"                                                             \n"
	"                      Frequencies (cm^-1)                    \n"
	"                      -------------------                    \n");

	nItem=0;
	// start from non-degenerate eigen values
	double ZPE=0.0;
	for(i=iStart; i < nDim; i++){

		// print frequecy
		printf("%11.1f ",lambda[i]);

		// compute zero point energy
		ZPE += lambda[i]/2.0;

		// 5 items each line
		nItem++; if(nItem%5==0) printf("\n");
	}
	// next line if not already done it
	if(nItem%5!=0) printf("\n");

	printf(
	"                                                             \n"
	"                                                             \n"
	"                    Zero Point Energy (ZPE)                  \n"
	"                    -----------------------                  \n"
	"               %18.5lf  Hartrees                  \n"
	"               %18.5lf  kcal/mol                  \n",

	ZPE*CM2HARTREE,ZPE*CM2KCALMOL);

	// free memory
	free(mwH);
	free(S);
	free(C);
	free(lambda);
}


//
// hessian_numerical: compute hesssian matrix using finite difference
//                    of the forces
//
// March 22, 2020 - Teepanis Chachiyo
//   Initial implementation
//
void hessian_numerical(
	int dbSize,                         // number of record in basis set db
	const struct GTOBasisSet_t *basisDB,// pointer to basis set database
	int nBasis,              // number of basis functions
	struct Molecule_t * mol, // pointer to molecule structure
	int nEA,                 // total number of spin up electrons
	int nEB,                 // total number of spin down electrons
	const double *CA,        // molecular alpha spin orbital
	const double *CB,        // molecular beta spin orbital 
	const double *eA,        // eigen values
	const double *eB,        // eigen values
	struct option_t *opt){   // options

	struct GTOBasis_t *gto;       // pointer to basis function storage
	int nDim;                     // degree of freedoms
	double *H;                    // Hessian
	double *dG;                   // change of gradient
	double *fx,*fy,*fz;           // forces acting on nuclei
	double *tCA, *teA, *tCB, *teB;// current molecular orbitals and eigen values
	double Etot;                  // total energy
	int savedSCFGuess;            // global guess option

	int nIter;                    // iteration index
	int i;                        // loop index

	// save global guess option
	savedSCFGuess = opt->SCFGuess;

	// load orbital from previos cycle to form initial guess
	opt->SCFGuess = SCFGUESS_CACB;

	// set degree of freedom
	nDim = mol->nAtom * 3;

	// memory allocation
	H      = calloc(nDim*nDim, sizeof(double));
	dG     = calloc(nDim,      sizeof(double));
	fx     = calloc(mol->nAtom,sizeof(double));
	fy     = calloc(mol->nAtom,sizeof(double));
	fz     = calloc(mol->nAtom,sizeof(double));

	// allocate memory for molecular orbitals and their eigen values
	gto   = genBasis(mol, &i, dbSize, basisDB);
	if(i!=nBasis){
		printf("hessian_numerical: error - inconsistent number of basis\n");
		exit(-1);
	}
	cleanGTOBasis(gto,nBasis);

	tCA = calloc(nBasis*nBasis,sizeof(double));
	teA = calloc(nBasis,sizeof(double));
	tCB = calloc(nBasis*nBasis,sizeof(double));
	teB = calloc(nBasis,sizeof(double));
	if(tCA==NULL || teA==NULL || tCB==NULL || teB==NULL){
		printf("hessian_numerical: error - cannot allocate memory\n");
		exit(-1);
	}

	///////////////////////////////////////////////
    // perform scf calculation and compute forces
	///////////////////////////////////////////////
	for(nIter=0; nIter < nDim; nIter++){
		printf(
		"                                                             \n"
		"                                                             \n"
		"-------------------------------------------------------------\n"
		"-----       HESSIAN NUMERICAL Step %5d  of %5d      -----\n"
		"-------------------------------------------------------------\n",
		nIter+1, nDim);
		fflush(stdout);

		///////////////////
		// forward step
		///////////////////
		switch(nIter%3){
		/* x */ case 0: mol->x[nIter/3] += HESSIAN_STEPSIZE; break;
		/* y */ case 1: mol->y[nIter/3] += HESSIAN_STEPSIZE; break;
		/* z */ case 2: mol->z[nIter/3] += HESSIAN_STEPSIZE; break;
		}

		// generate basis function
		gto   = genBasis(mol, &nBasis, dbSize, basisDB);

		// local initial guess
		for(i=0; i < nBasis*nBasis; i++){
			tCA[i] = CA[i];
			tCB[i] = CB[i];
		}

		// scf calculation
		Etot = uhf(nBasis, gto, mol, 
		    get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		    tCA, tCB, teA, teB, opt);
		if(Etot==0.0){
			printf("hessian_numerical: error SCF calculation did not converge\n");
			exit(-1);
		}

		// compute force
		uhf_force(nBasis, gto, mol, 
		      get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		      tCA, tCB, teA, teB, opt, fx, fy, fz);

		// free memory
		cleanGTOBasis(gto,nBasis);

		// construct change of gradient vector
		for(i=0; i < mol->nAtom; i++){
			dG[i*3+0] = -fx[i];
			dG[i*3+1] = -fy[i];
			dG[i*3+2] = -fz[i];
		}

		// reset molecular geometry
		switch(nIter%3){
		/* x */ case 0: mol->x[nIter/3] -= HESSIAN_STEPSIZE; break;
		/* y */ case 1: mol->y[nIter/3] -= HESSIAN_STEPSIZE; break;
		/* z */ case 2: mol->z[nIter/3] -= HESSIAN_STEPSIZE; break;
		}

		///////////////////
		// backward step
		///////////////////
		switch(nIter%3){
		/* x */ case 0: mol->x[nIter/3] -= HESSIAN_STEPSIZE; break;
		/* y */ case 1: mol->y[nIter/3] -= HESSIAN_STEPSIZE; break;
		/* z */ case 2: mol->z[nIter/3] -= HESSIAN_STEPSIZE; break;
		}

		// generate basis function
		gto   = genBasis(mol, &nBasis, dbSize, basisDB);

		// local initial guess
		for(i=0; i < nBasis*nBasis; i++){
			tCA[i] = CA[i];
			tCB[i] = CB[i];
		}

		// scf calculation
		Etot = uhf(nBasis, gto, mol, 
		    get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		    tCA, tCB, teA, teB, opt);
		if(Etot==0.0){
			printf("hessian_numerical: error SCF calculation did not converge\n");
			exit(-1);
		}

		// compute force
		uhf_force(nBasis, gto, mol, 
		      get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		      tCA, tCB, teA, teB, opt, fx, fy, fz);

		// free memory
		cleanGTOBasis(gto,nBasis);

		// construct change of gradient vector
		for(i=0; i < mol->nAtom; i++){
			dG[i*3+0] -= -fx[i];
			dG[i*3+1] -= -fy[i];
			dG[i*3+2] -= -fz[i];
		}

		// finite difference to compute hessian
		for(i=0; i < nDim; i++)
		H[nIter*nDim + i] = dG[i]/2.0/HESSIAN_STEPSIZE;

		// reset molecular geometry
		switch(nIter%3){
		/* x */ case 0: mol->x[nIter/3] += HESSIAN_STEPSIZE; break;
		/* y */ case 1: mol->y[nIter/3] += HESSIAN_STEPSIZE; break;
		/* z */ case 2: mol->z[nIter/3] += HESSIAN_STEPSIZE; break;
		}
	}

	// output
	hessian_output(nDim, H, mol);

	// free memory
	free(tCA);
	free(teA);
	free(tCB);
	free(teB);

	free(H);
	free(dG);
	free(fx);
	free(fy);
	free(fz);

	// set the guess option to original value
	opt->SCFGuess = savedSCFGuess;
}
