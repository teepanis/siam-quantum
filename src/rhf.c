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
#include "basis.h"
#include "matrix.h"
#include "lin.h"
#include "int.h"
#include "util.h"
#include "option.h"
#include "rhf.h"

// normalizeC : normalized eigen vector 
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static void normalizeC(int nBasis, double *S, double *C){
	double sum;
	int i,j,k;

	for(k=0; k < nBasis; k++){
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

// getDMatrix : compute density matrix
//
// 2008 - Teepanis Chachiyo
//    Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//    Not using static anymore
//
void rhf_getDMatrix(int nBasis, int nOcc, double *C, double *P){
	int i,j,k;
	double sum;

	for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			sum = 0.0;
			for(k=0; k < nOcc; k++){
				sum += C[k*nBasis+i] * C[k*nBasis+j];
			}
			P[i*nBasis+j] = 2.0*sum;
		}
	return;
}

// getGMatrix : compute G = 2J-K matrix element
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// March 13, 2019 - Teepanis Chachiyo
//    Add parallel support
//
static void getGMatrix(
	int nBasis,
	struct GTOBasis_t *gto,
	double *Schwarz,
	double *P,
	double *G){

	int i,j;
	
	// reset to zero
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++)
		G[i*nBasis+j] = 0.0;

	///////////////////////////////////
	// perform G matrix computation  //
	///////////////////////////////////
	// There are several choices available //
	//
	// 1) Simple straightforward by slow
	//GTO_2JK_Matrix_NoSymm(nBasis, P, gto, Schwarz, SCHWARZ_CUTOFF, G);
	//
	// 2) Expoit permutation symmetry of 2 electron integrals
	//GTO_2JK_Matrix(nBasis, P, gto, Schwarz, SCHWARZ_CUTOFF, G);
	//
	// 3) Group 2e integral into shells and compute the integral
	// in batches at the same time, and use permutation symmetry too.
	// This method is the most efficient we have so far.
	GTO_2JK_Matrix_Symm_Shell(nBasis, P, gto, Schwarz, SCHWARZ_CUTOFF, G);
			
	// divided by 2 according to canonical definition
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++)
		G[i*nBasis+j] = 0.5*G[i*nBasis+j];
}

// getEtotal : compute total energy this is equal to
// electronic part + nuclei part
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
static double getEtotal(
	int nBasis, 
	struct Molecule_t *mol,
	double *P, 
	double *F, 
	double *H){

	double E=0.0;
	int i,j;

	// compute nuclei repulsion energy
	E += nuclei_coulomb(mol);

	// include electron energy
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		E += (H[i*nBasis+j]+F[i*nBasis+j])*P[i*nBasis+j]*0.5;
	}

	return E;
}


// eval_chi : evaluate ith basis function at
// a Cartesian point (x,y,z)
//
// Dec 26, 2009 - Theerapon Kamla
// Fix potential bug with pow(0,0) tobe pow_int(0,0) instead
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
double eval_chi(int i, struct GTOBasis_t *gto,
                double x, double y, double z){

	int n;
	double sum=0;
	double r2;

	r2 = (x-gto[i].x0)*(x-gto[i].x0) + 
	     (y-gto[i].y0)*(y-gto[i].y0) +
	     (z-gto[i].z0)*(z-gto[i].z0);

	for(n=0; n < gto[i].nContract; n++)
		sum += gto[i].norm[n]*gto[i].coef[n]*exp(-gto[i].exp[n]*r2);
	sum = sum * pow_int(x-gto[i].x0, gto[i].l)
	          * pow_int(y-gto[i].y0, gto[i].m)
	          * pow_int(z-gto[i].z0, gto[i].n);
	return sum;
}

// rhf_rho : computes and return electron density at
// specified Cartesian point (x,y,z)
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
double rhf_rho(int nBasis,              // number of basis function
               struct GTOBasis_t *gto,  // function structure
               double *P,               // density matrix
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
		printf("rhf_rho - error cannot allocate memory\n");
		exit(-1);
	}

	// evaluate chi at a point
	for(i=0; i < nBasis; i++) chi[i] = eval_chi(i,gto,x,y,z);

	// evaluate density
	for(i=0; i < nBasis; i++)
	for(j=0; j < i; j++)
		sum += 2.0*P[i*nBasis+j]*chi[i]*chi[j];

	for(i=0; i < nBasis; i++)
		sum += P[i*nBasis+i]*chi[i]*chi[i];

	// clean up memory
	free(chi);

	return sum;
}

// rhf_mo : computes and return molecular orbital at
// specified Cartesian point (x,y,z)
//
// Dec 2009 - Theerapon Khamla
//    Initial implementation
//
// Mar 05, 2010 - Teepanis Chachiyo
//    Bug fix, invalid range of molecular orbital should be
//    "n >= nBasis" not "n > nBasis" as before.
//
double rhf_mo(int nBasis,              // number of basis function
              struct GTOBasis_t *gto,  // function structure
              double *C,               // molecular orbital
              int n,                   // orbital index (1st orbital is zero)
              double x, double y, double z){
	int i;
	double sum;

	// sanitize
	if(n >= nBasis || n < 0){
		printf("rho_mo: error invalid molecular orbital index\n");
		exit(EXIT_FAILURE);
	}

	// evaluate molecular orbital
	sum = 0.0;
	for(i=0; i < nBasis; i++){
		sum += C[n*nBasis+i]*eval_chi(i,gto,x,y,z);
	}

	return sum;

}

// rhf : carries out Restricted Hartree-Fock calculation until
// convergence is reached. It returns total electronic of the
// systems.
//
// In this version, the e-e integral is calculated using direct
// integration every time it is needed. They are not stored in
// the memory, which is not a good idea!
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
double rhf(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nE,                  // total number of electrons
	double *C,               // returned molecular orbitals
	double *e,               // returned eigen values
	struct option_t *opt){   // global option

	int nOcc;             // number of occupy orbital
	double Etot=0.0;      // total electronic energy
	double dE=0.0;        // energy change
	double avgdP=0.0;     // average change in density matrix
	double gamma=1.0;     // update coefficient

	// matrix elements
	double *P=NULL;      // density matrix
	double *prevP=NULL;  // previous iteration P matrix
	double *G=NULL;      // EE matrix
	double *F=NULL;      // fock matrix
	double *H=NULL;      // h core matrix
	double *T=NULL;      // kinetic matrix
	double *V=NULL;      // nuclei potential matrix
	double *S=NULL;      // overlap matrix
	double *Schwarz=NULL;// Schwarz upper bound matrix

	int i,j,iter=0;

	FILE *fd;            // file descriptor for density matrix

	// report
	printf("Requested Restricted Hartree-Fock calculations\n");
	printf("There are %d electrons in the density matrix\n", nE);

	// check number of electrons
	if(nE % 2 != 0){
		printf("RHF requires even number of electrons\n");
		exit(EXIT_FAILURE);
	}

	// number of occupy orbital
	nOcc = nE/2;

#define ALLOCATE(P)                                      \
P = calloc(nBasis*nBasis, sizeof(double));               \
if(P==NULL){                                             \
	printf("rhf: Error - Cannot allocate memory\n"); \
	exit(EXIT_FAILURE);                              \
}

	// memory allocation
	ALLOCATE(P);
	ALLOCATE(prevP);
	ALLOCATE(G);
	ALLOCATE(F);
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
	"Computing 1-electron matrix elements ... \n"
	"    H - Core Hamiltonian                 \n"
	"    S - Overlap Matrix                   \n");

	// get kinetic matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		T[i*nBasis+j] = GTO_kinetic(i,j,gto);
		// symmetrize matrix
		T[j*nBasis+i] = T[i*nBasis+j];
	}

	// get nuclei matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		V[i*nBasis+j] = GTO_nuclei(i,j,gto,mol);
		// symmetrize matrix
		V[j*nBasis+i] = V[i*nBasis+j];
	}

	// get overlap matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j <=i; j++){
		// compute explicitly
		S[i*nBasis+j] = GTO_overlap(i,j,gto);
		// symmetrize matrix
		S[j*nBasis+i] = S[i*nBasis+j];
	}

	// build Hcore 
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
		H[i*nBasis+j] = T[i*nBasis+j] + V[i*nBasis+j];
	}

	// It might be possible, as initial guess, to reduce the strength of
	// V matrix so that the electrons are farther apart even in the
	// non-interacting calculation.
	// Teepanis Chachiyo - July 10, 2010
	//

	//
	// Diagonalize core hamiltonian to guess density matrix
	//
	fflush(stdout);
	printf("Diagonalizing H for initial density matrix ...\n");
	gen_sym_eigen(nBasis, H, S, e, C);
	normalizeC(nBasis, S, C);
	rhf_getDMatrix(nBasis, nOcc, C, P);

	// compute initial E total
	Etot = getEtotal(nBasis, mol, P, H, H);
	printf("Non-interacting electronic energy  %.8f Hartree\n", Etot);

	printf("Processing 2E integrals ...\n");
	printf("Schwarz inequality screening cut-off %.8E\n", SCHWARZ_CUTOFF);
	Schwarz = create_Schwarz(nBasis, gto);

	// load density matrix if requested
	if(opt->loadDMatrix){

		// openfile
		fd=fopen(opt->DMatrixFile,"r");
		if(fd==NULL){
			printf("rhf - error cannot open file %s for reading\n", 
			       opt->DMatrixFile);
		}

		// loop and read
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			if(fscanf(fd, "%lf", P+i*nBasis+j)!=1){
				printf("rhf - error reading density matrix\n");
				exit(-1);
			}
		}

		// close file
		fclose(fd);

		// report status
		printf("Density matrix loaded from %s\n", opt->DMatrixFile);
	}

	// scf loop
	printf("Drag coefficient %f\n", opt->SCFDrag);
	printf("SCFMax %d iterations\n", opt->SCFMax);
	printf("Enter SCF loop ... \n");
	printf(
	"Iteration  Total Energy [Hartrees]  RMSD Density\n"
	"------------------------------------------------\n");
	fflush(stdout);
	
	// oscillation drag coefficient
	gamma = opt->SCFDrag;

	do{
		iter++;
		
		getGMatrix(nBasis, gto, Schwarz, P, G);

		// update Fock-matrix
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++)
			F[i*nBasis+j] = H[i*nBasis+j] + G[i*nBasis+j];
		gen_sym_eigen(nBasis, F, S, e, C);
		normalizeC(nBasis, S, C);

		// save previous P matrix
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++)
			prevP[i*nBasis+j] = P[i*nBasis+j];

		// get new P matrix
		rhf_getDMatrix(nBasis, nOcc, C, P);

		// compute average change
		avgdP     = 0.0;
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++)
			avgdP +=  (P[i*nBasis+j]-prevP[i*nBasis+j])
			         *(P[i*nBasis+j]-prevP[i*nBasis+j]);
		avgdP = sqrt( avgdP/nBasis/nBasis);

		// compute energy and energ difference
		dE     = Etot;
		Etot   = getEtotal(nBasis, mol, P, F, H);
		dE     = Etot - dE;

		printf(" %5d %20.8f %20.4E\n", iter, Etot, avgdP);

		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++)
			P[i*nBasis+j] = prevP[i*nBasis+j] + 
			                gamma *(     P[i*nBasis+j]
			                        -prevP[i*nBasis+j]
			                       );
		
		// check if we have reached scfmax limit
		if(iter >= opt->SCFMax) break;

		// flush output
		fflush(stdout);

	}while(fabs(dE) > 1.0E-6 || avgdP > 1.0E-6);

	// report
	if(iter >= opt->SCFMax){
		printf("SCF have not converged because iter >= SCFMax\n");
	}else{
		printf("Done SCF Loop Total Energy is %20.8f Hartrees\n", Etot);
		printf("Koopmann Ionization Energy is %20.8f Hartrees  [ %.4f eV ]\n",
		       e[nOcc-1], e[nOcc-1] * HARTREE2EV);
	}
	fflush(stdout);

	// save density matrix if requested
	if(opt->saveDMatrix){

		// openfile
		fd=fopen(opt->DMatrixFile,"w");
		if(fd==NULL){
			printf("rhf - error cannot open file %s for writing\n", 
			       opt->DMatrixFile);
		}

		// loop and save each element
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++){
			if((i*nBasis+j)%5==0) fprintf(fd,"\n");
			fprintf(fd," %15.8E ", P[i*nBasis+j]);
		}
		
		// close file
		fclose(fd);

		// report status
		printf("Density matrix saved to %s\n", opt->DMatrixFile);
	}

	// clean memory
	free(Schwarz);
	free(P);
	free(prevP);
	free(G);
	free(F);
	free(H);
	free(T);
	free(V);
	free(S);

	return Etot;
}
