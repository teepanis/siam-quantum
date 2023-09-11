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
#include <math.h>
#include <string.h>

#include "check.h"

#include "lin.h"
#include "util.h"
#include "int.h"

#define CHECKPOINT_VERSION "1.2.5"

//
// load_checkpoint_molecule : Load molecular structure from checkpoint file
//     
//
// Oct 5, 2012 - Teepanis Chachiyo
//      Initial implementation and testing
//
struct Molecule_t *load_checkpoint_molecule(const struct option_t *opt){

	FILE *fd;                // checkpoint file pointer
	struct Molecule_t *mol;  // pointer to created molecule
	int Q;                   // total charge
	int nAtom;               // number of atoms
	int i;                   // generic loop

	// open file
	fd=fopen(opt->CheckFile,"r");
	if(fd==NULL){
		printf("load_checkpoint_molecule - error "
		       "cannot open file %s for reading\n",
		       opt->CheckFile);
		exit(-1);
	}

	// search for total charge keyword
	if(findf(fd,1,"TOTAL_CHARGE")==EOF){
		printf("load_checkpoint_molecule - "
		       "error Cannot find TOTAL_CHARGE keyword\n");
		       exit(-1);
	}

	// read total charge
	if(fscanf(fd,"%d", &Q)!=1){
		printf("load_checkpoint_molecule - error cannot read total charge\n");
		exit(-1);
	}

	// search for molecular structure keyword
	if(findf(fd,1,"XYZ_FORMAT")==EOF){
		printf("load_checkpoint_molecule - "
		       "error cannot find XYZ_FORMAT keyword\n");
		       exit(-1);
	}

	// allocate memory
	mol = (struct Molecule_t *) calloc(1, sizeof(struct Molecule_t));
	if(mol == NULL){
		printf("load_checkpoint_molecule - error annot allocate memory\n");
		exit(-1);
	}

	// read number of atoms
	if(fscanf(fd,"%d",&nAtom)!=1){
		printf("load_checkpont_molecule - error cannot read number of atom\n");
		exit(-1);
	}

	// validate the number of atoms
	if(nAtom <= 0){
		printf("load_checkpont_molecule - error invalid range of nAtom\n");
		exit(-1);
	}

	// assign the number of atom
	mol->nAtom = nAtom;

	// allocate atomic number and coordinates
	mol->Z = calloc(mol->nAtom, sizeof(int));
	if(mol->Z == NULL){
		printf("load_checkpont_molecule - error annot allocate atomic number\n");
		exit(EXIT_FAILURE);
	}
	mol->x = calloc(mol->nAtom, sizeof(double));
	mol->y = calloc(mol->nAtom, sizeof(double));
	mol->z = calloc(mol->nAtom, sizeof(double));
	if(mol->x == NULL || mol->y == NULL || mol->z == NULL){
		printf("load_checkpoint_molecule - error annot allocate coordinates\n");
		exit(EXIT_FAILURE);
	}

	// read atom type and position 
	for(i=0; i < mol->nAtom; i++){
		if(fscanf(fd, "%d %lf %lf %lf", mol->Z+i,
		                                mol->x+i,
		                                mol->y+i,
		                                mol->z+i)!=4){
			printf("load_checkpoint_molecule - error reading coordinate\n");
			exit(EXIT_FAILURE);
		}
	}

	// set total charge
	mol->Q = Q;	

	// close file
	fclose(fd);

	return mol;
}

//
// load_checkpoint_multiplicity : Load and return multiplicity from checkpoint
//
// Oct 5, 2012 - Teepanis Chachiyo
//      Initial implementation and testing
//
int load_checkpoint_multiplicity(const struct option_t *opt){
	FILE *fd;                // checkpoint file pointer
	int M;                   // multiplicity

	// open file
	fd=fopen(opt->CheckFile,"r");
	if(fd==NULL){
		printf("load_checkpoint_multiplicity - error "
		       "cannot open file %s for reading\n",
		       opt->CheckFile);
		exit(-1);
	}

	// search for multiplicity keyword
	if(findf(fd,1,"SPIN_MULTIPLICITY")==EOF){
		printf("load_checkpoint_multiplicity - "
		       "error Cannot find SPIN_MULTIPLICITY keyword\n");
		       exit(-1);
	}

	// read spin multiplicity
	if(fscanf(fd, "%d", &M)!=1){
		printf("load_checkpoint_multiplicity - error reading multiplicity\n");
		exit(-1);
	}

	// close file
	fclose(fd);

	return M;
}

//
// load_checkpoint_basis : Load and return pointer to basis function data
// structure. The number of basis function is also set through *nBasis
//
// Oct 5, 2012 - Teepanis Chachiyo
//     Initial implementation and testing
//
struct GTOBasis_t *load_checkpoint_basis(const struct option_t *opt, int *nBasis){
	FILE *fd;                // checkpoint file pointer
	struct GTOBasis_t *gto;  // pointer to basis data structure

	// open file
	fd=fopen(opt->CheckFile,"r");
	if(fd==NULL){
		printf("load_checkpoint_basis - error cannot open file %s for reading\n",
		       opt->CheckFile);
		exit(-1);
	}

	// search for basis functions keyword
	if(findf(fd,1,"SECTION_BASIS_FUNCTION_BEGIN")==EOF){
		printf("load_checkpoint_basis - "
		       "error Cannot find SECTION_BASIS_FUNCTION_BEGIN keyword\n");
		       exit(-1);
	}

	// read basis function
	gto = readGTOInput(fd, nBasis);

	// close file
	fclose(fd);

	return gto;
}

//
// load_checkpoint_orbital : Load molecular orbital from checkpoint.
// It also sets the UHF or RHF flag inside struct optiont_t
//
// Oct 5, 2012 - Teepanis Chachiyo
//      Initial implementation and testing
//
double load_checkpoint_orbital(int nBasis,                  // number of basis function
                               double *CA,                  // alpha orbitals
                               double *CB,                  // beta orbitals 
                               double *eA,                  // alpha eigen values
                               double *eB,                  // beta eigen values
                               const struct option_t *opt){ // global options

	FILE *fd;                // checkpoint file pointer
	int i,j;                 // generic loop
	double Etot;             // total energy;

	// open file
	fd=fopen(opt->CheckFile,"r");
	if(fd==NULL){
		printf("load_checkpoint_orbital - error "
		       "cannot open file %s for reading\n",
		       opt->CheckFile);
		exit(-1);
	}

	// search for total energy keyword
	if(findf(fd,1,"TOTAL_ENERGY")==EOF){
		printf("load_checkpoint_orbital - "
		       "error Cannot find TOTAL_ENERGY keyword\n");
		       exit(-1);
	}
	if(fscanf(fd,"%lf",&Etot)!=1){
		printf("load_checkpoint_orbital - error reading total energy\n");
		exit(-1);
	}

	// search for alpha eigen value keyword
	if(findf(fd,1,"EIGEN_VALUE_ALPHA")==EOF){
		printf("load_checkpoint_orbital - "
		       "error Cannot find EIGEN_VALUE_ALPHA keyword\n");
		       exit(-1);
	}
	for(i=0; i < nBasis; i++)
		if(fscanf(fd,"%lf",eA+i)!=1){
			printf("load_checkpoint_orbital - error reading alpha eigen\n");
			exit(-1);
		}

	// search for beta eigen value keyword
	if(findf(fd,1,"EIGEN_VALUE_BETA")==EOF){
		printf("load_checkpoint_orbital - "
		       "error Cannot find EIGEN_VALUE_BETA keyword\n");
		       exit(-1);
	}
	for(i=0; i < nBasis; i++)
		if(fscanf(fd,"%lf",eB+i)!=1){
			printf("load_checkpoint_orbital - error reading beta eigen\n");
			exit(-1);
		}

	// search for alpha orbital keyword
	if(findf(fd,1,"ORBITAL_ALPHA")==EOF){
		printf("load_checkpoint_orbital - "
		       "error Cannot find ORBITAL_ALPHA keyword\n");
		       exit(-1);
	}
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++)
		if(fscanf(fd,"%lf",CA+(i*nBasis+j))!=1){
			printf("load_checkpoint_orbital - error reading alpha orbital\n");
			exit(-1);
		}

	// search for beta orbital keyword
	if(findf(fd,1,"ORBITAL_BETA")==EOF){
		printf("load_checkpoint_orbital - "
		       "error Cannot find ORBITAL_BETA keyword\n");
		       exit(-1);
	}
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++)
		if(fscanf(fd,"%lf",CB+(i*nBasis+j))!=1){
			printf("load_checkpoint_orbital - error reading beta orbital\n");
			exit(-1);
		}

	// close file
	fclose(fd);

	return Etot;
}

// guess_checkpoint_orbital : read molecular orbital from checkpoint file
// and fit it to current basis set
//
// Oct 5, 2012 - Teepanis Chachiyo
//     Initial implementation and testing
//
// Orbital Projection Technique [See. Documentation]
// ----------------------------
//
void guess_checkpoint_orbital(
	int nBasis,                   // number of basis function
	const struct GTOBasis_t *gto, // basis function
    const struct Molecule_t *mol, // molecular structure structure
	const struct option_t *opt,   // global options
	const double *S,              // overlap matrix 
	double *CA,                   // returned guess alpha orbital
	double *CB){                  // returned guess beta  orbital

	int nChkBasis;             // checkpoint number of basis
	struct Molecule_t *chkMol; // checkpoint molecular structure
	struct GTOBasis_t *chkGTO; // checkpoint basis functions
	double *Y;                 // overlap matrix between current and checkpoint
	double *B;                 // right-hand-side linear vector
	double *chkCA, *chkCB;     // checkpoint orbitals
	double *chkEA, *chkEB;     // checkoint eigenvalues
	int i,j,k;                 // generic loop
	int iCnt, jCnt;            // contracted function 
	double sum=0.0;            // integral sum

	// load molecule from checkpoint and check if it is the same
	chkMol = load_checkpoint_molecule(opt);
	if(chkMol->nAtom != mol->nAtom){
		printf("guess_checkpoint_orbital - error "
		       "checkpoint contains different number of atoms\n");
		exit(-1);
	}
	for(i=0; i < chkMol->nAtom; i++)
		if(chkMol->Z[i] != mol->Z[i]){
			printf("guess_checkpoint_orbital - error "
			       "molecule types in checkpoint do not match\n");
			exit(-1);
		}

	// load basis function from checkpoint
	chkGTO = load_checkpoint_basis(opt, &nChkBasis);

	// validate the number of basis
	if(nChkBasis > nBasis){
		printf("guess_checkpoint_orbital - error "
		       "cannot use larger basis as guess\n");
		exit(-1);
	}

	// allocate overlap matrix
	Y     = calloc(nBasis*nChkBasis,    sizeof(double));
	B     = calloc(nBasis*nChkBasis,    sizeof(double));
	chkCA = calloc(nChkBasis*nChkBasis, sizeof(double));
	chkCB = calloc(nChkBasis*nChkBasis, sizeof(double));
	chkEA = calloc(nChkBasis,           sizeof(double));
	chkEB = calloc(nChkBasis,           sizeof(double));
	if(  Y==NULL ||     B==NULL || 
	 chkCA==NULL || chkCB==NULL || chkEA==NULL || chkEB==NULL){
		printf("guess_checkpoint_orbital - error cannot allocate memory\n");
		exit(-1);
	}

	// load checkpoint orbitals
	load_checkpoint_orbital(nChkBasis, chkCA, chkCB, chkEA, chkEB, opt);

	// compute overlap matrix between current and checkpoint basis
	for(i=0; i < nBasis; i++)
	for(j=0; j < nChkBasis; j++){

		// looping over contracted functions
		sum = 0.0;
		for(iCnt=0; iCnt <    gto[i].nContract; iCnt++)
		for(jCnt=0; jCnt < chkGTO[j].nContract; jCnt++){
			sum += gto[i].coef[iCnt] * chkGTO[j].coef[jCnt] *
			       gto[i].norm[iCnt] * chkGTO[j].norm[jCnt] *
		   overlap(gto[i].exp[iCnt],    gto[i].l,     gto[i].m,     gto[i].n,
	                                    gto[i].x0,    gto[i].y0,    gto[i].z0,
		        chkGTO[j].exp[jCnt], chkGTO[j].l,  chkGTO[j].m,  chkGTO[j].n,
			                         chkGTO[j].x0, chkGTO[j].y0, chkGTO[j].z0);
		}
		Y[i*nChkBasis+j] = sum;
	}

	// compute B matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j < nChkBasis; j++){
		sum = 0.0;
		for(k=0; k < nChkBasis; k++)
			sum += Y[i*nChkBasis+k]*chkCA[j*nChkBasis+k];
		B[i*nChkBasis+j] = sum;
	}

	// compute CA
	linear_solver(nBasis, nChkBasis, S, B);

	for(i=0; i < nChkBasis; i++)
	for(j=0; j < nBasis; j++)
		CA[i*nBasis+j] = B[j*nChkBasis+i];

	// compute B matrix
	for(i=0; i < nBasis; i++)
	for(j=0; j < nChkBasis; j++){
		sum = 0.0;
		for(k=0; k < nChkBasis; k++)
			sum += Y[i*nChkBasis+k]*chkCB[j*nChkBasis+k];
		B[i*nChkBasis+j] = sum;
	}

	// compute CB
	linear_solver(nBasis, nChkBasis, S, B);

	for(i=0; i < nChkBasis; i++)
	for(j=0; j < nBasis; j++)
		CB[i*nBasis+j] = B[j*nChkBasis+i];

	// clean memory usage
	free(Y);
	free(B);
	free(chkCA);
	free(chkCB);
	free(chkEA);
	free(chkEB);
	cleanGTOBasis(chkGTO, nChkBasis);
	cleanMolecule(chkMol);
}

// save_checkpoint : Save information to checkpoint file
//
// Oct 5, 2012 - Teepanis Chachiyo
//     Initial implementaion and testing
//
void save_checkpoint(int nBasis,             // number of basis
                     struct GTOBasis_t *gto, // pointer to basis function
                     struct Molecule_t *mol, // pointer to molecule structure
                     double Etot,            // total energy
                     double dE,              // energy convergence
                     double rmsDP,           // density convergence
                   	 double *CA,             // molecular alpha spin orbital
	                 double *CB,             // molecular beta spin orbital 
	                 double *eA,             // eigen values for alpha spin
	                 double *eB,             // eigen values for beta  spin
                     struct option_t *opt ){ // global options

	FILE *fd;       // checkpoint file pointer
	int i;          // generic loop
	char time[256]; // time stamp

	// open file
	fd=fopen(opt->CheckFile,"w");
	if(fd==NULL){
		printf("save_checkpoint - error cannot open file %s for writing\n", 
		       opt->CheckFile);
	}

	////////////
	// header
	////////////
	fprintf(fd,"SIAM_QUANTUM_CHECKPOINT\n");
	fprintf(fd,"FORMAT_VERSION    %s\n", CHECKPOINT_VERSION);
	time_str(256, time);
	fprintf(fd, "TIME %s\n", time);
	
	/////////////////////////////
	// molecular specification	
	/////////////////////////////
	
	// keyword
	fprintf(fd, "SECTION_MOLECULAR_SPECIFICATION_BEGIN\n");

	// charge
	fprintf(fd, "TOTAL_CHARGE   %d\n", mol->Q);

	// structure in Bohr
	fprintf(fd, "XYZ_FORMAT\n");
	fprintf(fd,"%d\n\n",mol->nAtom);
	for(i=0; i < mol->nAtom; i++){
		fprintf(fd, "%3d %16.8E %16.8E %16.8E\n",
		        mol->Z[i],
		        mol->x[i],
		        mol->y[i],
		        mol->z[i]);
	}

	fprintf(fd, "SECTION_MOLECULAR_SPECIFICATION_END\n");

	/////////////////////////////////////
	// save basis function information
	/////////////////////////////////////

	// keyword
	fprintf(fd, "SECTION_BASIS_FUNCTION_BEGIN\n");

	// info
	printGTOBasis(fd, nBasis, gto);

	fprintf(fd, "SECTION_BASIS_FUNCTION_END\n");

	///////////////////////////
	// electronic properties
	///////////////////////////
	
	// keyword
	fprintf(fd, "SECTION_ELECTRONIC_STRUCTURE_BEGIN\n");

	// spin multiplicity
	fprintf(fd, "SPIN_MULTIPLICITY   %d\n", opt->multiplicity);

	// scf method
	     if(opt->method==METHOD_HF)        fprintf(fd, "SCF_METHOD   HF\n");
	else if(opt->method==METHOD_MP2)       fprintf(fd, "SCF_METHOD   HF\n");
	else if(opt->method==METHOD_HF2xDFT)   fprintf(fd, "SCF_METHOD   HF\n");
	else if(opt->method==METHOD_SVWN5)     fprintf(fd, "SCF_METHOD   LDA - exchange: Slater   correlation: VWN5\n");
	else if(opt->method==METHOD_SCHACHIYO) fprintf(fd, "SCF_METHOD   LDA - exchange: Slater   correlation: Chachiyo\n");
	else if(opt->method==METHOD_CHACHIYO)  fprintf(fd, "SCF_METHOD   GGA - exchange: Chachiyo correlation: Chachiyo\n");
#ifdef LIBXC
	else if(opt->method==METHOD_LIBXC)     fprintf(fd, "SCF_METHOD   GGA - LIBXC\n");
#endif
	else if(opt->method==METHOD_HALF)      fprintf(fd, "SCF_METHOD   HALF\n");
	else{
		printf("save_checkpoint: Error - not recognize scf method\n");
		exit(-1);
	}

	// total energy
	fprintf(fd,"TOTAL_ENERGY   %20.8f\n", Etot);

	// energy convergence
	fprintf(fd,"ENERGY_CONVERGENCE   %15.3E\n", dE);

	// density convergence
	fprintf(fd,"DENSITY_CONVERGENCE   %15.3E\n", rmsDP);

	// alpha eigen values
	fprintf(fd, "EIGEN_VALUE_ALPHA");
	for(i=0; i < nBasis; i++){
		if(i%5==0) fprintf(fd,"\n");
		fprintf(fd, " %16.8E ", eA[i]);
	}
	fprintf(fd,"\n");

	// beta eigen value
	fprintf(fd, "EIGEN_VALUE_BETA");
	for(i=0; i < nBasis; i++){
		if(i%5==0) fprintf(fd,"\n");
		fprintf(fd, " %16.8E ", eB[i]);
	}
	fprintf(fd,"\n");

	// alpha orbitals
	fprintf(fd, "ORBITAL_ALPHA");
	for(i=0; i < nBasis*nBasis; i++){
		if(i%5==0) fprintf(fd,"\n");
		fprintf(fd, " %16.8E ", CA[i]);
	}
	fprintf(fd,"\n");

	// beta orbitals
	fprintf(fd, "ORBITAL_BETA");
	for(i=0; i < nBasis*nBasis; i++){
		if(i%5==0) fprintf(fd,"\n");
		fprintf(fd, " %16.8E ", CB[i]);
	}
	fprintf(fd,"\n");

	fprintf(fd, "SECTION_ELECTRONIC_STRUCTURE_END\n");

	// close file
	fclose(fd);
}

