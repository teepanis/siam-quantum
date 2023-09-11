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
#include <stdio.h>
#include <string.h>

#include "mp2.h"
#include "int.h"
#include "matrix.h"
#include "lin.h"
#include "rpc.h"
#include "util.h"

// the number of core electrons for each atom type
#define MP2_MaxNobleGasCoreAtom 54
int MP2_NobleGasCore[] = {0,
 0,                                                 0,
 2, 2,                               2, 2, 2, 2, 2, 2,
10,10,                              10,10,10,10,10,10,
18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36};


// getFrozenCore_RHF : computes the number of frozen core orbitals
// for the closed-shell system
//
// Nove 5, 2012 - Teepanis Chachiyo
//      Initial implementation and testing
//
int getFrozenCore_RHF(const struct Molecule_t *mol){

	int nCore = 0;               // number of frozen core orbitals
	int a;                       // atom index

	for(a=0; a < mol->nAtom; a++){

		// do we have the number of core available ?
		if(mol->Z[a] > MP2_MaxNobleGasCoreAtom){
			printf("getFrozenCore_RHF - error no available noble gas atom info.\n");
			exit(-1);
		}

		// yes we do, then accumulate
		nCore += MP2_NobleGasCore[mol->Z[a]];
	}

	// divided by 2 to get the number of orbital for closed-shell case
	nCore = nCore/2;

	return nCore;
}


// mp2_gen_aqbj_NoSymm : generate coulomb integral between (aq|bj)
// here a,b are molecular orbitals and q,j are basis functions.
//
// Oct 20, 2012 - Teepanis Chachiyo
//    Initial implementation and testing
//
double * mp2_gen_aqbj_NoSymm(
	int nBasis,                     // number of basis function
	int nOcc,                       // number of occpupied orbitals
	int nCore,                      // number of core orbitals to ignore
	const double *C,                // molecular orbitals
	const struct GTOBasis_t *gto){  // basis function structure

	int nCorr;        // number of correlated orbitals to include
	double *schwarz;  // schwarz inequality for basis function
	double *aqbjBank; // pointer to integral storage
	int p,q,i,j;      // basis function index
	int a,b;          // orbital index
	double *ptr;      // pointer to current value
	double EE;        // calculated 2e integral

	// compute number of correlated orbitals
	nCorr = nOcc - nCore;

	// generate schwarz inequality between basis function
	schwarz = create_Schwarz(nBasis, gto);

	// memory allocation
	if((aqbjBank=calloc((nCorr*(nCorr+1))/2*nBasis*nBasis,sizeof(double)))==NULL){
		printf("mp2_gen_aqbj_NoSymm - error cannot allocate memory\n");
		exit(-1); 
	}

	for(p=0; p < nBasis; p++)
	for(i=0; i < nBasis; i++)
	for(q=0; q < nBasis; q++)
	for(j=0; j < nBasis; j++)

		// schwarz screening at basis function level
		if(schwarz[p*nBasis+q]*schwarz[i*nBasis+j] > MP2_SCHWARZ_CUTOFF){

			// compute two-electron integral
			EE = contr_eri(
			     gto[p].nContract,
			     gto[p].exp, gto[p].coef, gto[p].norm,
			     gto[p].x0,        gto[p].y0,  gto[p].z0,
			     gto[p].l,         gto[p].m,   gto[p].n,

			     gto[q].nContract,
			     gto[q].exp, gto[q].coef, gto[q].norm,
			     gto[q].x0,        gto[q].y0,  gto[q].z0,
			     gto[q].l,         gto[q].m,   gto[q].n,

			     gto[i].nContract,
			     gto[i].exp, gto[i].coef, gto[i].norm,
			     gto[i].x0,        gto[i].y0,  gto[i].z0,
			     gto[i].l,         gto[i].m,   gto[i].n,

			     gto[j].nContract,
			     gto[j].exp, gto[j].coef, gto[j].norm,
			     gto[j].x0,        gto[j].y0,  gto[j].z0,
			     gto[j].l,         gto[j].m,   gto[j].n);

			// loop through correlated orbitals
			ptr = aqbjBank + (q*nBasis+j)*nCorr*(nCorr+1)/2;
			for(a=0; a < nCorr; a++)
			for(b=0; b <= a; b++){
				*ptr += C[(a+nCore)*nBasis+p]*C[(b+nCore)*nBasis+i]*EE;
				 ptr++;
			}
		}
	

	// free memory
	free(schwarz);

	return aqbjBank;
}


// mp2_report_memory : reports the how much memory is required
//
// Nov 5, 2012 - Teepanis Chachiyo
//     Initial implementation and testing
//
void mp2_report_memory(double nByte){
	if(nByte > 1.0E12)
	printf("The number of memory required (in terabytes):  %4.1f\n", nByte/1.0E12);
	else if(nByte > 1.0E9)
	printf("The number of memory required (in gigabytes):  %4.1f\n", nByte/1.0E9);
	else if(nByte > 1.0E6)
	printf("The number of memory required (in megabytes):  %4.1f\n", nByte/1.0E6);
	else
	printf("The number of memory required (in kilobytes):  %4.1f\n", nByte/1.0E3);
}


// mp2_rhf_aqbj : computes MP2 correction energy for the case of closed-shell
// orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for the closed-shell system is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[arbs] { (ar|bs) * ( 2*(ar|bs) - (as|br) ) / (ea+eb-er-es)   }
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. See Document for details.
//
// Oct 20, 2012 - Teepanis Chachiyo
//    Initial implementation and testing
//
double mp2_rhf_aqbj(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorr;                 // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int q,j;                   // basis function index
	double *ptr;               // pointer to current aqbj value
	double *arbs,*asbr;        // computed integral
	double *arbsPtr, *asbrPtr; // pointer to above matrix
	int nVir;                  // number of virtual orbitals
	const double *eVir;        // pointer to virtual orbital eigen value
	const double *eCor;        // correlated orbital eigen value
	double sum;                // accumulated mp2 correction
	double contrib;            // current contribution to sum
	double *aqbjBank;          // 2e storage between orbitals and basis functions
	double aqbj;               // current aqbj value
	double rsCoef, srCoef;     // coefficient for compute arbs and asbr

	// check if we have restricted hartree-fock
	if(opt->restricted==0){
		printf("mp2_rhf - error only restricted MP2 is supported\n");
		exit(-1);
	}

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVir  = nBasis - nOcc;
	nCorr = nOcc   - nCore; 

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----               MP2 ENERGY CORRECTION               -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:          %4d\n", nCore);
	printf("The number of included occupied orbitals:      %4d\n", nCorr);
	printf("The number of included virtual  orbitals:      %4d\n", nVir);

	// memory requirement
	sum = 0.5*(nCorr*(nCorr+1))*nBasis*nBasis*sizeof(double);
	mp2_report_memory(sum);
	fflush(stdout);

	// allocate computed integral
	if((arbs=calloc((nCorr*(nCorr+1))/2,sizeof(double)))==NULL ||
	   (asbr=calloc((nCorr*(nCorr+1))/2,sizeof(double)))==NULL){
		printf("mp2_rhf - error cannot allocate memory\n");
		exit(-1);
	}

	// locate virtual orbital and correlated orbitals eigen value
	eVir = e + nOcc;
	eCor = e + nCore;

	// generate aqbj integral
	//aqbjBank = mp2_gen_aqbj_NoSymm(nBasis, nOcc, nCore, C, gto);
	aqbjBank = mp2_gen_aqbj_ShellSet(nBasis, nOcc, nCore, C, gto);

	// compute mp2 energy
	sum = 0.0;
	for(r=0; r < nVir; r++)
	for(s=0; s <= r; s++){

		// reset arbs and asbr matrix to zero
		memset(arbs,0,sizeof(double)*(nCorr*(nCorr+1))/2);
		memset(asbr,0,sizeof(double)*(nCorr*(nCorr+1))/2);

		// compute arbs and asbr
		ptr = aqbjBank;
		for(q=0; q < nBasis; q++)
		for(j=0; j < nBasis; j++){

			rsCoef = C[(r+nOcc)*nBasis+q]*C[(s+nOcc)*nBasis+j];
			srCoef = C[(s+nOcc)*nBasis+q]*C[(r+nOcc)*nBasis+j];

			// screening
			if(fabs(rsCoef) < MP2_SCHWARZ_CUTOFF && 
			   fabs(srCoef) < MP2_SCHWARZ_CUTOFF){
				ptr += nCorr*(nCorr+1)/2;
				continue;
			}

			arbsPtr = arbs;
			asbrPtr = asbr;
			for(a=0; a < nCorr; a++)
			for(b=0; b <= a; b++){
				aqbj = *ptr;
				ptr++;
				*arbsPtr += rsCoef*aqbj; arbsPtr++;
				*asbrPtr += srCoef*aqbj; asbrPtr++;
			}
		}

		// add contribution to mp2 energy
		arbsPtr = arbs;
		asbrPtr = asbr;
		for(a=0; a < nCorr; a++)
		for(b=0; b <= a; b++){
			contrib =   (*arbsPtr)*(*arbsPtr)
			          + (*asbrPtr)*(*asbrPtr)
			          - (*arbsPtr)*(*asbrPtr);
			contrib = contrib*2.0/(eCor[a]+eCor[b]-eVir[r]-eVir[s]);
			if(r==s) contrib = contrib/2.0;
			if(a!=b) contrib = contrib*2.0;
			sum += contrib;
			arbsPtr++; asbrPtr++;
		}
	}

	printf("\nMP2 Energy Correction is %.6f Hartrees\n", sum);
	printf(
       "-------------------------------------------------------------\n");

	// free memory
	free(aqbjBank);
	free(arbs);
	free(asbr);

	return sum;
}


// mp2_rhf_aqbj : computes MP2 correction energy for the case of closed-shell
// orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for the closed-shell system is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[arbs] { (ar|bs) * ( 2*(ar|bs) - (as|br) ) / (ea+eb-er-es)   }
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. See Document for details.
//
// Nov 5, 2012 - Teepanis Chachiyo
//    Initial implementation and testing. As of today, the version aqij
//    is faster and more memory friendly than the aqbj
//
// Nov 5, 2012 - Teepanis Chachiyo
//    The qaij array is now using hierarchy i->j->a->q which will make 
//    the loop locally accessing memory and is much faster.
//
// Feb 2013 - Teepanis Chachiyo
//    Use the quartet version instead
//
// Sep 17 2014 - Teepanis Chachiyo
//    Check is maxMem is to low to do even one round
//
double mp2_rhf_aqij(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorr;                 // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int q,i,j;                 // basis function index
	double *aqbj, *arbj;       // intermediate integrals
	double Cbi, arbjSum;       // intermediate values
	double arbs, asbr;         // integral needed to compute MP2 energy
	double *ptr;               // pointer to current aqbj value
	int nVir;                  // number of virtual orbitals
	const double *eVir;        // pointer to virtual orbital eigen value
	const double *eCor;        // correlated orbital eigen value
	double sum;                // accumulated mp2 correction
	double contrib;            // current contribution to sum
	double *aqijBank;          // 2e storage between orbitals and basis functions
	double **IJm;              // IJm mapping matrix
	double *Schwarz;           // schwarz inequality
	int nEE;                   // number of integral index
	int maxCorr;               // maximum correlated orbitals fits in memory
	int aCycle;                // a index for this round
	int aCount;                // counter for a index round

	// check if we have restricted hartree-fock
	if(opt->restricted==0){
		printf("mp2_rhf_aqij - error only restricted MP2 is supported\n");
		exit(-1);
	}

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVir  = nBasis - nOcc;
	nCorr = nOcc   - nCore; 

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----               MP2 ENERGY CORRECTION               -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:          %4d\n", nCore);
	printf("The number of included occupied orbitals:      %4d\n", nCorr);
	printf("The number of included virtual  orbitals:      %4d\n", nVir);

	// memory requirement
	Schwarz = create_Schwarz(nBasis, gto);
	for(nEE=0,i=0; i < nBasis; i++)
	for(j=0; j <= i; j++)
		if(Schwarz[i*nBasis+j] >= MP2_SCHWARZ_CUTOFF) nEE++;

	// Check the case where opt->maxMem is too low to do even one round
	// Teepanis - 17 September 2014
	if(opt->maxMem*1.0E6 < nEE*nBasis*sizeof(double)){
		printf("mp2_rhf_aqij - error maxMem must be at least %.1f MB\n",
		      nEE*nBasis*sizeof(double)*1.0E-6);
		exit(-1);
	}

	maxCorr = (int)floor(opt->maxMem*1.0E6/nEE/nBasis/sizeof(double));
	if(maxCorr > nCorr) maxCorr = nCorr;

	sum = nEE*nBasis*maxCorr*sizeof(double);
	mp2_report_memory(sum);

	// report the number of cycles
	printf("The number of cycles to compute (pq|ij):       %4d\n",
	       nCorr==0?0: (int)ceil((float)nCorr/maxCorr));

	fflush(stdout);

	// locate virtual orbital and correlated orbitals eigen value
	eVir = e + nOcc;
	eCor = e + nCore;

	// allocate computed integral
	if((aqbj=calloc(nBasis*nBasis,sizeof(double )))==NULL ||
	   (arbj=calloc(nVir*nBasis,sizeof(double )))==NULL ||
	   (IJm=calloc(nBasis*nBasis,sizeof(double*)))==NULL){
		printf("mp2_rhf_aqij - error cannot allocate memory\n");
		exit(-1);
	}

	/*
	// Assuming we have unlimited memory, below is the simple code
	// that computes MP2 energy in a single cycles

	// generate aqbj integral
	aqijBank = mp2_gen_aqij_ShellSet(nBasis, nCorr, nCore, C, gto);

	// build IJm mapping matrix
	for(nEE=0,i=0; i < nBasis; i++)
	for(j=0; j <= i; j++){
		if(Schwarz[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) IJm[i*nBasis+j] = NULL;
		else{
			IJm[i*nBasis+j] = aqijBank + nEE*nBasis*maxCorr;
			nEE++;
		}
		IJm[j*nBasis+i] = IJm[i*nBasis+j];
	}

	sum = 0.0;
	for(a=0; a < nCorr; a++)
	for(b=0; b <= a; b++){

		memset(aqbj,0,sizeof(double)*nBasis*nBasis);
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++)
			if((ptr=IJm[i*nBasis+j]))
				for(q=0; q < nBasis; q++)
					aqbj[j*nBasis+q] += C[(b+nCore)*nBasis+i]*ptr[a*nBasis+q];

		memset(arbj,0,sizeof(double)*nVir*nBasis);
		for(r=0; r < nVir; r++)
		for(j=0; j < nBasis; j++)
		for(q=0; q < nBasis; q++)
			arbj[r*nBasis+j] += C[(r+nOcc)*nBasis+q]*aqbj[j*nBasis+q];

		for(r=0; r < nVir; r++)
		for(s=0; s <= r; s++){
			arbs = asbr = 0.0;
			for(j=0; j < nBasis; j++){
				arbs += C[(s+nOcc)*nBasis+j]*arbj[r*nBasis+j];
				asbr += C[(r+nOcc)*nBasis+j]*arbj[s*nBasis+j];
			}

			contrib = arbs*arbs + asbr*asbr - arbs*asbr;
			contrib = contrib*2.0/(eCor[a]+eCor[b]-eVir[r]-eVir[s]);
			if(r==s) contrib = contrib/2.0;
			if(a!=b) contrib = contrib*2.0;
			sum += contrib;
		}

	}
	*/

	// compute MP2 energy in multiple cycles
	sum = 0.0;
	for(a=0,aCount=0; a < nCorr; a+=maxCorr,aCount++){

		// compute number of orbitals to compute for this cycles
		if((nCorr-a) < maxCorr) maxCorr = (nCorr-a);

		// generate aqbj integral
		//aqijBank = mp2_gen_aqij_ShellSet(nBasis, maxCorr, nCore+a, C, gto);
		aqijBank = mp2_gen_aqij_Quartet(nBasis, maxCorr, nCore+a, C, gto, Schwarz);

		// build IJm mapping matrix
		for(nEE=0,i=0; i < nBasis; i++)
		for(j=0; j <= i; j++){
			if(Schwarz[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) IJm[i*nBasis+j] = NULL;
			else{
				IJm[i*nBasis+j] = aqijBank + nEE*nBasis*maxCorr;
				nEE++;
			}
			IJm[j*nBasis+i] = IJm[i*nBasis+j];
		}

		// compute MP2 energy this this cycle
		for(aCycle=0; aCycle < maxCorr; aCycle++)
		for(b=0; b <= (aCycle+a); b++){
	
			memset(aqbj,0,sizeof(double)*nBasis*nBasis);
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++)
			if((ptr=IJm[i*nBasis+j])){
				Cbi = C[(b+nCore)*nBasis+i];
				if(fabs(Schwarz[i*nBasis+j]*Cbi) > MP2_SCHWARZ_CUTOFF)
				for(q=0; q < nBasis; q++)
					aqbj[j*nBasis+q] += Cbi*ptr[aCycle*nBasis+q];
			}

			memset(arbj,0,sizeof(double)*nVir*nBasis);
			for(r=0; r < nVir; r++)
			for(j=0; j < nBasis; j++){
				for(arbjSum=0.0,q=0; q < nBasis; q++)
					arbjSum += C[(r+nOcc)*nBasis+q]*aqbj[j*nBasis+q];
				arbj[r*nBasis+j] = arbjSum;
			}

			for(r=0; r < nVir; r++)
			for(s=0; s <= r; s++){
				arbs = asbr = 0.0;
				for(j=0; j < nBasis; j++){
					arbs += C[(s+nOcc)*nBasis+j]*arbj[r*nBasis+j];
					asbr += C[(r+nOcc)*nBasis+j]*arbj[s*nBasis+j];
				}
	
				contrib = arbs*arbs + asbr*asbr - arbs*asbr;
				contrib = contrib*2.0/(eCor[a+aCycle]+eCor[b]-eVir[r]-eVir[s]);
				if(r==s) contrib = contrib/2.0;
				if((a+aCycle)!=b) contrib = contrib*2.0;
				sum += contrib;
			}
	
		}

		// free memory
		free(aqijBank);

		// progress report
		printf("The number of cycles already finished:         %4d\n",aCount+1);
		fflush(stdout);

	}

	printf("\nMP2 Energy Correction is %.6f Hartrees\n", sum);
	printf(
       "-------------------------------------------------------------\n");

	// free memory
	mp2_gen_aqij_Quartet(0, 0, 0, NULL, NULL, NULL);
	free(aqbj);
	free(arbj);
	free(IJm);
	free(Schwarz);

	return sum;
}


// mp2_rhf_aqij_Parallel_Contribute: this is the parallel version of the 
// mp2_rhf_aqij subroutine. It only compute partial mp2 contributions.
//
// Mar 6, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
double mp2_rhf_aqij_Parallel_Contribute(
	int a,                         // starting orbital index
	int maxCorr,                   // number of correlated orbitals
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const double *Schwarz,         // schwarz inequality
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorr;                 // number of correlated orbitals
	int b,r,s;                 // molecular orbitals
	int q,i,j;                 // basis function index
	double *aqbj, *arbj;       // intermediate integrals
	double Cbi, arbjSum;       // intermediate values
	double arbs, asbr;         // integral needed to compute MP2 energy
	double *ptr;               // pointer to current aqbj value
	int nVir;                  // number of virtual orbitals
	const double *eVir;        // pointer to virtual orbital eigen value
	const double *eCor;        // correlated orbital eigen value
	double sum;                // accumulated mp2 correction
	double contrib;            // current contribution to sum
	double *aqijBank;          // 2e storage between orbitals and basis functions
	double **IJm;              // IJm mapping matrix
	int nEE;                   // number of integral index
	int aCycle;                // a index for this round

	// reset call
	if(nBasis==0){
		mp2_gen_aqij_Quartet(0, 0, 0, NULL, NULL, NULL);
		return 0.0;
	}

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVir  = nBasis - nOcc;
	nCorr = nOcc   - nCore; 

	// locate virtual orbital and correlated orbitals eigen value
	eVir = e + nOcc;
	eCor = e + nCore;

	// allocate computed integral
	if((aqbj=calloc(nBasis*nBasis,sizeof(double )))==NULL ||
	   (arbj=calloc(nVir*nBasis,sizeof(double )))==NULL ||
	   (IJm=calloc(nBasis*nBasis,sizeof(double*)))==NULL){
		printf("mp2_rhf_aqij_Parallel_Contribute - error cannot allocate memory\n");
		exit(-1);
	}

	// compute MP2 energy in a single cycle
	sum = 0.0;

	// generate aqbj integral
	aqijBank = mp2_gen_aqij_Quartet(nBasis, maxCorr, nCore+a, C, gto, Schwarz);

	// build IJm mapping matrix
	for(nEE=0,i=0; i < nBasis; i++)
	for(j=0; j <= i; j++){
		if(Schwarz[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) IJm[i*nBasis+j] = NULL;
		else{
			IJm[i*nBasis+j] = aqijBank + nEE*nBasis*maxCorr;
			nEE++;
		}
		IJm[j*nBasis+i] = IJm[i*nBasis+j];
	}

	// compute MP2 energy this this cycle
	for(aCycle=0; aCycle < maxCorr; aCycle++)
	for(b=0; b <= (aCycle+a); b++){

		memset(aqbj,0,sizeof(double)*nBasis*nBasis);
		for(i=0; i < nBasis; i++)
		for(j=0; j < nBasis; j++)
		if((ptr=IJm[i*nBasis+j])){
			Cbi = C[(b+nCore)*nBasis+i];
			if(fabs(Schwarz[i*nBasis+j]*Cbi) > MP2_SCHWARZ_CUTOFF)
			for(q=0; q < nBasis; q++)
				aqbj[j*nBasis+q] += Cbi*ptr[aCycle*nBasis+q];
		}

		memset(arbj,0,sizeof(double)*nVir*nBasis);
		for(r=0; r < nVir; r++)
		for(j=0; j < nBasis; j++){
			for(arbjSum=0.0,q=0; q < nBasis; q++)
				arbjSum += C[(r+nOcc)*nBasis+q]*aqbj[j*nBasis+q];
			arbj[r*nBasis+j] = arbjSum;
		}

		for(r=0; r < nVir; r++)
		for(s=0; s <= r; s++){
			arbs = asbr = 0.0;
			for(j=0; j < nBasis; j++){
				arbs += C[(s+nOcc)*nBasis+j]*arbj[r*nBasis+j];
				asbr += C[(r+nOcc)*nBasis+j]*arbj[s*nBasis+j];
			}

			contrib = arbs*arbs + asbr*asbr - arbs*asbr;
			contrib = contrib*2.0/(eCor[a+aCycle]+eCor[b]-eVir[r]-eVir[s]);
			if(r==s) contrib = contrib/2.0;
			if((a+aCycle)!=b) contrib = contrib*2.0;
			sum += contrib;
		}

	}

	// free memory
	free(aqijBank);

	// free memory
	free(aqbj);
	free(arbj);
	free(IJm);

	return sum;
}


// mp2_rhf_aqij_Parallel: this is the parallel version of the mp2_rhf_aqij 
// subroutine. 
//
// Mar 6, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
// Sep 17, 2014 - Teepanis Chachiyo
//     Check maxMem
//
double mp2_rhf_aqij_Parallel(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorr;                 // number of correlated orbitals
	int a;                     // molecular orbitals
	int i,j;                   // basis function index
	int nVir;                  // number of virtual orbitals
	double sum;                // accumulated mp2 correction
	double mp2;                // mp2 contribution
	double *Schwarz;           // schwarz inequality
	int nEE;                   // number of integral index
	int maxCorr;               // maximum correlated orbitals fits in memory
	int aCount;                // counter for a index round

	// check if we have restricted hartree-fock
	if(opt->restricted==0){
		printf("mp2_rhf_aqij - error only restricted MP2 is supported\n");
		exit(-1);
	}

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVir  = nBasis - nOcc;
	nCorr = nOcc   - nCore; 

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----               MP2 ENERGY CORRECTION               -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:          %4d\n", nCore);
	printf("The number of included occupied orbitals:      %4d\n", nCorr);
	printf("The number of included virtual  orbitals:      %4d\n", nVir);

	// memory requirement
	Schwarz = create_Schwarz(nBasis, gto);
	for(nEE=0,i=0; i < nBasis; i++)
	for(j=0; j <= i; j++)
		if(Schwarz[i*nBasis+j] >= MP2_SCHWARZ_CUTOFF) nEE++;

	// Check the case where opt->maxMem is too low to do even one round
	// Teepanis - 17 September 2014
	if(opt->maxMem*1.0E6 < nEE*nBasis*sizeof(double)){
		printf("mp2_rhf_aqij_Parallel - error maxMem must be at least %.1f MB\n",
		      nEE*nBasis*sizeof(double)*1.0E-6);
		exit(-1);
	}

	// the number of rounds per one cpu
	for(i=1; ceil((double)nCorr/(i*opt->nCPU))*nEE*nBasis*sizeof(double)
	         > opt->maxMem*1.0E6; i++);
	
	// each correlation orbital per round per cpu
	maxCorr = (int)ceil((float)nCorr/(i*opt->nCPU));

	sum = nEE*nBasis*maxCorr*sizeof(double);
	mp2_report_memory(sum);

	// report the number of cycles
	printf("The number of cycles to compute (pq|ij):       %4d\n",
	       nCorr==0?0: (int)ceil((float)nCorr/maxCorr));

	fflush(stdout);

	//
	// parallel version
	//
	int *status;             // status for each cpu
	int alldone;             // all idle flag

	// allocate memory
	status=calloc(opt->nCPU,sizeof(int));
	if(status==NULL){
		printf("uhf - error cannot allocate memory\n");
		exit(-1);
	}

	sum = 0.0;
	for(a=0,aCount=0; a < nCorr; ){

		// reset status to idle
		for(i=(opt->nCPU-1);i>=0;i--) status[i] = RPC_IDLE;

		// go thru all cpu
		do{
			//
			// remote call to all cpu except childID=0
			//
			for(i=(opt->nCPU-1);i>0;i--){

				// skip if done
				if(status[i]==RPC_DONE) continue;

				if(status[i]==RPC_IDLE){
					// compute number of orbitals to compute for this cycles
					if((nCorr-a) < maxCorr) maxCorr = (nCorr-a);

					// call remote proxy
					status[i] = rpc_mp2_rhf_aqij_Parallel_Contribute(status[i], i, &mp2,
					            a, maxCorr, nBasis, nOcc, e, C, Schwarz, gto, mol, opt);
	
					// increment starting orbital
					a += maxCorr;

					// increment counter
					if(maxCorr != 0) aCount++;
				}else
					status[i] = rpc_mp2_rhf_aqij_Parallel_Contribute(status[i], i, &mp2,
					            a, maxCorr, nBasis, nOcc, e, C, Schwarz, gto, mol, opt);

				// add contribution
				if(status[i]==RPC_DONE) sum += mp2;
			}
	
			//
			// local call for chilID=0
			//
			if(status[0]==RPC_IDLE){
	
				// compute number of orbitals to compute for this cycles
				if((nCorr-a) < maxCorr) maxCorr = (nCorr-a);
	
				sum += mp2_rhf_aqij_Parallel_Contribute(a, maxCorr, nBasis, nOcc,
				                                        e, C, Schwarz, gto, mol, opt);
				status[0] = RPC_DONE;
	
				// increment starting orbital
				a += maxCorr;

				// increment counter
				if(maxCorr != 0) aCount++;
			}
	
			// check if all done
			alldone=1;
			for(i=(opt->nCPU-1);i>=0;i--) if(status[i] != RPC_DONE) alldone=0;

		}while(!alldone);

		// progress report
		printf("The number of cycles already finished:         %4d\n",aCount);
		fflush(stdout);

	}

	printf("\n"
	       "MP2 Energy Correction is %12.6f Hartree\n"
	       "MP2 Energy Correction is %12.6f eV\n", sum, HARTREE2EV*sum);
	printf("-------------------------------------------------------------\n");

	//
	// clean memory
	//

	// reset status to idle
	for(i=(opt->nCPU-1);i>=0;i--) status[i] = RPC_IDLE;

	// go thru all cpu
	do{
		// remote call to all cpu except childID=0
		for(i=(opt->nCPU-1);i>0;i--){
			status[i] = rpc_mp2_rhf_aqij_Parallel_Contribute(status[i], i, &mp2,
			            0, 0, 0, 0, NULL, NULL, NULL, NULL, mol, opt);
		}

		// local call for chilID=0
		if(status[0]==RPC_IDLE){
			mp2_rhf_aqij_Parallel_Contribute(0, 0, 0, 0, NULL, NULL, NULL, NULL, mol, opt);
			status[0] = RPC_DONE;
		}

		// check if all done
		alldone=1;
		for(i=(opt->nCPU-1);i>=0;i--) if(status[i] != RPC_DONE) alldone=0;

	}while(!alldone);

	free(status);
	free(Schwarz);

	return sum;
}


// mp2_rhf_direct : computes MP2 correction energy for the case of closed-shell
// orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for the closed-shell system is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[arbs] { (ar|bs) * ( 2*(ar|bs) - (as|br) ) / (ea+eb-er-es)   }
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. See Document for details.
//
// Nov 5, 2012 - Teepanis Chachiyo
//    Initial implementation and testing. The "direct" version means that
//    no integral are stored in memory at all making it the slowest possible
//    method of computing MP2. But it is good for testing and comparing with
//    other methods.
//
double mp2_rhf_direct(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorr;                 // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int p,q,i,j;               // basis function index
	double EE;                 // compute electron integral
	int nVir;                  // number of virtual orbitals
	const double *eVir;        // pointer to virtual orbital eigen value
	const double *eCor;        // correlated orbital eigen value
	double arbs, asbr;         // current integral value
	double *schwarz;           // schwarz inequality at basis level
	double Cap,Crq,Cbi,Csj,Csq,Crj; // transform coefficients
	double Spq,Sij;                 // schwarz screening
	double contrib, sum=0.0;        // acuumulated MP2 energy
	unsigned long long nEE=0;

	// check if we have restricted hartree-fock
	if(opt->restricted==0){
		printf("mp2_rhf_direct - error only restricted MP2 is supported\n");
		exit(-1);
	}

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVir  = nBasis - nOcc;
	nCorr = nOcc   - nCore; 

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----               MP2 ENERGY CORRECTION               -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:         %4d\n", nCore);
	printf("The number of included occupied orbitals:     %4d\n", nCorr);
	printf("The number of included virtual  orbitals:     %4d\n", nVir);
	fflush(stdout);

	// locate virtual orbital and correlated orbitals eigen value
	eVir = e + nOcc;
	eCor = e + nCore;

	// generate schwarz inequality between basis function
	schwarz = create_Schwarz(nBasis, gto);

	// compute mp2 energy
	for(a=0; a < nCorr; a++)
	for(b=0; b <= a; b++)
	for(r=0; r < nVir; r++)
	for(s=0; s <= r; s++){

		// reset values
		arbs = asbr = 0.0;

		for(p=0; p < nBasis; p++)
		for(q=0; q < nBasis; q++){
			Cap = C[(a+nCore)*nBasis+p];
			Crq = C[(r+ nOcc)*nBasis+q];
			Csq = C[(s+ nOcc)*nBasis+q];
			Spq = schwarz[p*nBasis+q];
			if(fabs(Cap*Crq*Spq) > MP2_SCHWARZ_CUTOFF ||
			   fabs(Cap*Csq*Spq) > MP2_SCHWARZ_CUTOFF){

				for(i=0; i < nBasis; i++)
				for(j=0; j < nBasis; j++){
					Cbi = C[(b+nCore)*nBasis+i];
					Csj = C[(s+ nOcc)*nBasis+j];
					Crj = C[(r+ nOcc)*nBasis+j];
					Sij = schwarz[i*nBasis+j];
					if( fabs(Cap*Crq*Cbi*Csj*Spq*Sij) > MP2_SCHWARZ_CUTOFF  ||
					    fabs(Cap*Csq*Cbi*Crj*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

						// compute two-electron integral
						EE = contr_eri(
						     gto[p].nContract,
						     gto[p].exp, gto[p].coef, gto[p].norm,
						     gto[p].x0,        gto[p].y0,  gto[p].z0,
						     gto[p].l,         gto[p].m,   gto[p].n,
	
						     gto[q].nContract,
						     gto[q].exp, gto[q].coef, gto[q].norm,
						     gto[q].x0,        gto[q].y0,  gto[q].z0,
						     gto[q].l,         gto[q].m,   gto[q].n,
	
						     gto[i].nContract,
						     gto[i].exp, gto[i].coef, gto[i].norm,
						     gto[i].x0,        gto[i].y0,  gto[i].z0,
						     gto[i].l,         gto[i].m,   gto[i].n,
	
						     gto[j].nContract,
						     gto[j].exp, gto[j].coef, gto[j].norm,
						     gto[j].x0,        gto[j].y0,  gto[j].z0,
						     gto[j].l,         gto[j].m,   gto[j].n);
	
						// counter
						nEE++;
	
						arbs += EE*Cap*Crq*Cbi*Csj;
						asbr += EE*Cap*Csq*Cbi*Crj;
					}
				}
			}
		}

		contrib = arbs*arbs + asbr*asbr - arbs*asbr;
		contrib = contrib*2.0/(eCor[a]+eCor[b]-eVir[r]-eVir[s]);
		if(r==s) contrib = contrib/2.0;
		if(a!=b) contrib = contrib*2.0;
		sum += contrib;		
	}

	printf("The number of (pq|ij) integral computed:  %lld\n", nEE);

	printf("\nMP2 Energy Correction is %.6f Hartrees\n", sum);
	printf(
       "-------------------------------------------------------------\n");

	// free memory
	free(schwarz);

	return sum;
}


// mp2_rhf_semi_direct : computes MP2 correction energy for the case of closed-shell
// orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for the closed-shell system is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[arbs] { (ar|bs) * ( 2*(ar|bs) - (as|br) ) / (ea+eb-er-es)   }
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. See Document for details.
//
// Nov 5, 2012 - Teepanis Chachiyo
//    Initial implementation and testing. The "semi_direct_aqbj" version means
//    the integral [aqbj]qj are stored. This method has nBasis*nBasis memory
//    requirement which is very small. It is a little bit faster than the fully
//    direct but still very slow compared to mp2_rhf_aqbj
//
double mp2_rhf_semi_direct_aqbj(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorr;                 // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int p,q,i,j;               // basis function index
	double EE;                 // compute electron integral
	int nVir;                  // number of virtual orbitals
	const double *eVir;        // pointer to virtual orbital eigen value
	const double *eCor;        // correlated orbital eigen value
	double arbs, asbr;         // current integral value
	double *aqbj;              // intermediate storage
	double *aqbjPtr;           // pointer to current aqbj 
	double *schwarz;           // schwarz inequality at basis level
	double Cap,Crq,Cbi,Csj,Csq,Crj; // transform coefficients
	double Spq,Sij;                 // schwarz screening
	double contrib, sum=0.0;        // acuumulated MP2 energy
	unsigned long long nEE=0;

	// check if we have restricted hartree-fock
	if(opt->restricted==0){
		printf("mp2_rhf_semi_direct_aqbj - error only restricted MP2 is supported\n");
		exit(-1);
	}

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVir  = nBasis - nOcc;
	nCorr = nOcc   - nCore; 

	// allocate intermediate storage
	if((aqbj=calloc(nBasis*nBasis,sizeof(double)))==NULL){
		printf("mp2_rhf_semi_direct_aqbj - error cannot allocate memory\n");
		exit(-1);
	}

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----               MP2 ENERGY CORRECTION               -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:     %3d\n", nCore);
	printf("The number of included occupied orbitals: %3d\n", nCorr);
	printf("The number of included virtual  orbitals: %3d\n", nVir);
	fflush(stdout);

	// locate virtual orbital and correlated orbitals eigen value
	eVir = e + nOcc;
	eCor = e + nCore;

	// generate schwarz inequality between basis function
	schwarz = create_Schwarz(nBasis, gto);

	for(a=0; a < nCorr; a++)
	for(b=0; b <= a; b++){

		// reset (aqbj)
		memset(aqbj,0,nBasis*nBasis*sizeof(double));

		// generate (aqbj)
		for(p=0; p < nBasis; p++)
		for(i=0; i < nBasis; i++){

			// load coefficients
			Cap = C[(a+nCore)*nBasis+p];
			Cbi = C[(b+nCore)*nBasis+i];

			// reset pointer
			aqbjPtr = aqbj;

			// screening
			if(fabs(Cap*Cbi) > MP2_SCHWARZ_CUTOFF)
			for(q=0; q < nBasis; q++)
			for(j=0; j < nBasis; j++){

				// load schwarz inequality values
				Spq = schwarz[p*nBasis+q];
				Sij = schwarz[i*nBasis+j];

				// screening
				if(fabs(Cap*Cbi*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

					// compute two-electron integral
					EE = contr_eri(
					     gto[p].nContract,
					     gto[p].exp, gto[p].coef, gto[p].norm,
					     gto[p].x0,        gto[p].y0,  gto[p].z0,
					     gto[p].l,         gto[p].m,   gto[p].n,

					     gto[q].nContract,
					     gto[q].exp, gto[q].coef, gto[q].norm,
					     gto[q].x0,        gto[q].y0,  gto[q].z0,
					     gto[q].l,         gto[q].m,   gto[q].n,

					     gto[i].nContract,
					     gto[i].exp, gto[i].coef, gto[i].norm,
					     gto[i].x0,        gto[i].y0,  gto[i].z0,
					     gto[i].l,         gto[i].m,   gto[i].n,

					     gto[j].nContract,
					     gto[j].exp, gto[j].coef, gto[j].norm,
					     gto[j].x0,        gto[j].y0,  gto[j].z0,
					     gto[j].l,         gto[j].m,   gto[j].n);

					// counter
					nEE++;

					*aqbjPtr += EE*Cap*Cbi;
				}

				// increment pointer index
				aqbjPtr++;
			}
		}

		// use (aq|bj) to generate (ar|bs) and (as|br)
		for(r=0; r < nVir; r++)
		for(s=0; s <= r; s++){

			//reset value
			arbs = asbr = 0.0;

			// reset pointer 
			aqbjPtr = aqbj;

			for(q=0; q < nBasis; q++)
			for(j=0; j < nBasis; j++){
				// load coefficient
				Crq = C[(r+ nOcc)*nBasis+q];
				Csq = C[(s+ nOcc)*nBasis+q];
				Csj = C[(s+ nOcc)*nBasis+j];
				Crj = C[(r+ nOcc)*nBasis+j];

				// compute integral
				arbs += (*aqbjPtr)*Crq*Csj;
				asbr += (*aqbjPtr)*Csq*Crj;

				// step up pointer
				aqbjPtr++;
			}

			// compute MP2 contributions
			contrib = arbs*arbs + asbr*asbr - arbs*asbr;
			contrib = contrib*2.0/(eCor[a]+eCor[b]-eVir[r]-eVir[s]);
			if(r==s) contrib = contrib/2.0;
			if(a!=b) contrib = contrib*2.0;
			sum += contrib;
		}

	}

	printf("The number of (pq|ij) integral computed:  %lld\n", nEE);

	printf("\nMP2 Energy Correction is %.6f Hartrees\n", sum);
	printf(
       "-------------------------------------------------------------\n");

	// free memory
	free(schwarz);
	free(aqbj);

	return sum;
}


// mp2_rhf_semi_direct : computes MP2 correction energy for the case of closed-shell
// orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for the closed-shell system is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[arbs] { (ar|bs) * ( 2*(ar|bs) - (as|br) ) / (ea+eb-er-es)   }
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. See Document for details.
//
// Nov 5, 2012 - Teepanis Chachiyo
//    Initial implementation and testing. The "semi_direct_aqij" version means
//    the integral [aqij]qij are stored. This method has nBasis*nBasis*nBasis
//    requirement which is sizable. It is faster than the semi_direct aqbj 
//    owing to the large memory requirement.
//
double mp2_rhf_semi_direct_aqij(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorr;                 // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int p,q,i,j;               // basis function index
	double EE;                 // compute electron integral
	int nVir;                  // number of virtual orbitals
	const double *eVir;        // pointer to virtual orbital eigen value
	const double *eCor;        // correlated orbital eigen value
	double arbs, asbr;         // current integral value
	double *aqij;              // nBasis^3 intermediate storage
	double *aqijPtr;           // pointer to current aqbj 
	double *schwarz;           // schwarz inequality at basis level
	double Cap,Crq,Cbi,Csj,Csq,Crj; // transform coefficients
	double Spq,Sij;                 // schwarz screening
	double contrib, sum=0.0;        // acuumulated MP2 energy
	unsigned long long nEE=0;

	// check if we have restricted hartree-fock
	if(opt->restricted==0){
		printf("mp2_rhf_semi_direct_aqij - error only restricted MP2 is supported\n");
		exit(-1);
	}

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVir  = nBasis - nOcc;
	nCorr = nOcc   - nCore; 

	// allocate intermediate storage
	if((aqij=calloc(nBasis*nBasis*nBasis,sizeof(double)))==NULL){
		printf("mp2_rhf_semi_direct_aqij - error cannot allocate memory\n");
		exit(-1);
	}

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----               MP2 ENERGY CORRECTION               -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:     %3d\n", nCore);
	printf("The number of included occupied orbitals: %3d\n", nCorr);
	printf("The number of included virtual  orbitals: %3d\n", nVir);
	fflush(stdout);

	// locate virtual orbital and correlated orbitals eigen value
	eVir = e + nOcc;
	eCor = e + nCore;

	// generate schwarz inequality between basis function
	schwarz = create_Schwarz(nBasis, gto);

	for(a=0; a < nCorr; a++){

		// reset (aqij)
		memset(aqij,0,nBasis*nBasis*nBasis*sizeof(double));

		// generate (aqij)
		for(p=0; p < nBasis; p++){

			// load coefficients
			Cap = C[(a+nCore)*nBasis+p];

			// reset pointer
			aqijPtr = aqij;

			// screening
			if(fabs(Cap) > MP2_SCHWARZ_CUTOFF)
			for(q=0; q < nBasis; q++)
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){

				// load schwarz inequality values
				Spq = schwarz[p*nBasis+q];
				Sij = schwarz[i*nBasis+j];

				// screening
				if(fabs(Cap*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

					// compute two-electron integral
					EE = contr_eri(
					     gto[p].nContract,
					     gto[p].exp, gto[p].coef, gto[p].norm,
					     gto[p].x0,        gto[p].y0,  gto[p].z0,
					     gto[p].l,         gto[p].m,   gto[p].n,

					     gto[q].nContract,
					     gto[q].exp, gto[q].coef, gto[q].norm,
					     gto[q].x0,        gto[q].y0,  gto[q].z0,
					     gto[q].l,         gto[q].m,   gto[q].n,

					     gto[i].nContract,
					     gto[i].exp, gto[i].coef, gto[i].norm,
					     gto[i].x0,        gto[i].y0,  gto[i].z0,
					     gto[i].l,         gto[i].m,   gto[i].n,

					     gto[j].nContract,
					     gto[j].exp, gto[j].coef, gto[j].norm,
					     gto[j].x0,        gto[j].y0,  gto[j].z0,
					     gto[j].l,         gto[j].m,   gto[j].n);

					// counter
					nEE++;

					*aqijPtr += EE*Cap;
				}

				// increment pointer index
				aqijPtr++;
			}
		}

		// use (aq|ij) to generate (ar|bs) and (as|br)
		for(b=0; b <= a; b++)
		for(r=0; r < nVir; r++)
		for(s=0; s <= r; s++){

			//reset value
			arbs = asbr = 0.0;

			// reset pointer 
			aqijPtr = aqij;

			for(q=0; q < nBasis; q++)
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				// load coefficient
				Cbi = C[(b+nCore)*nBasis+i];
				Crq = C[(r+ nOcc)*nBasis+q];
				Csq = C[(s+ nOcc)*nBasis+q];
				Csj = C[(s+ nOcc)*nBasis+j];
				Crj = C[(r+ nOcc)*nBasis+j];

				// compute integral
				arbs += (*aqijPtr)*Crq*Cbi*Csj;
				asbr += (*aqijPtr)*Csq*Cbi*Crj;

				// step up pointer
				aqijPtr++;
			}

			// compute MP2 contributions
			contrib = arbs*arbs + asbr*asbr - arbs*asbr;
			contrib = contrib*2.0/(eCor[a]+eCor[b]-eVir[r]-eVir[s]);
			if(r==s) contrib = contrib/2.0;
			if(a!=b) contrib = contrib*2.0;
			sum += contrib;
		}

	}

	printf("The number of (pq|ij) integral computed:  %lld\n", nEE);

	printf("\nMP2 Energy Correction is %.6f Hartrees\n", sum);
	printf(
       "-------------------------------------------------------------\n");

	// free memory
	free(schwarz);
	free(aqij);

	return sum;
}


// mp2_uhf_direct : computes MP2 correction energy for the case of unrestricted
// orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for general case is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[a<b, r<s] { <ab||rs> / (ea+eb-er-es)  }
//
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. The integral is over spin-orbitals.
// Once the MP2 energy has been expressed explicitly in terms of spatial-orbital
// it assumes the form:
//
//   MP2 = SUM[a<b, r<s] {( (ar|bs) - (as|br) )^2 /(ea+eb-er-es) } 
//         a,b,r,s all spin up
//
//      +  SUM[a<b, r<s] {( (ar|bs) - (as|br) )^2 /(ea+eb-er-es) } 
//         a,b,r,s all spin dn
//
//      +  SUM[a b, r s] {( (ar|bs) )^2           /(ea+eb-er-es) } 
//         a,r are spin up and r,s are spin dn
//
// See research note page 45 for details
//
// Dec 22, 2014 - Teepanis Chachiyo
//    Initial implementation and testing. Adapted from mp2_rhf_direct
//
double mp2_uhf_direct(
	int nBasis,                    // number of basis function
	int nEA,                       // number of occupied spin up orbitals
	int nEB,                       // number of occupied spin dn orbitals
	const double *eA,              // eigen values
	const double *eB,              // eigen values
	const double *CA,              // molecular orbitals
	const double *CB,              // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorrA,nCorrB;         // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int p,q,i,j;               // basis function index
	double EE;                 // compute electron integral
	int nVirA,nVirB;           // number of virtual orbitals
	const double *eVirA,*eVirB;// pointer to virtual orbital eigen value
	const double *eCorA,*eCorB;// correlated orbital eigen value
	double arbs, asbr;         // current integral value
	double *schwarz;           // schwarz inequality at basis level
	double Cap,Crq,Cbi,Csj,Csq,Crj; // transform coefficients
	double Spq,Sij;                 // schwarz screening
	double contrib, sum=0.0;        // acuumulated MP2 energy
	unsigned long long nEE=0;

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVirA  = nBasis - nEA;
	nCorrA = nEA    - nCore;
	nVirB  = nBasis - nEB;
	nCorrB = nEB    - nCore;

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----        UNRESTRICTED MP2 ENERGY CORRECTION         -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:      %4d %4d\n", nCore, nCore);
	printf("The number of included occupied orbitals:  %4d %4d\n", nCorrA, nCorrB);
	printf("The number of included virtual  orbitals:  %4d %4d\n", nVirA, nVirB);
	fflush(stdout);

	// locate virtual orbital and correlated orbitals eigen value
	eVirA = eA + nEA;
	eCorA = eA + nCore;
	eVirB = eB + nEB;
	eCorB = eB + nCore;

	// generate schwarz inequality between basis function
	schwarz = create_Schwarz(nBasis, gto);

	// compute mp2 energy: spin-up contribution
	for(a=0; a < nCorrA; a++)
	for(b=a+1; b < nCorrA; b++)
	for(r=0; r < nVirA; r++)
	for(s=r+1; s < nVirA; s++){

		// reset values
		arbs = asbr = 0.0;

		for(p=0; p < nBasis; p++)
		for(q=0; q < nBasis; q++){
			Cap = CA[(a+nCore)*nBasis+p];
			Crq = CA[(r+  nEA)*nBasis+q];
			Csq = CA[(s+  nEA)*nBasis+q];
			Spq = schwarz[p*nBasis+q];
			if(fabs(Cap*Crq*Spq) > MP2_SCHWARZ_CUTOFF ||
			   fabs(Cap*Csq*Spq) > MP2_SCHWARZ_CUTOFF){

				for(i=0; i < nBasis; i++)
				for(j=0; j < nBasis; j++){
					Cbi = CA[(b+nCore)*nBasis+i];
					Csj = CA[(s+  nEA)*nBasis+j];
					Crj = CA[(r+  nEA)*nBasis+j];
					Sij = schwarz[i*nBasis+j];
					if( fabs(Cap*Crq*Cbi*Csj*Spq*Sij) > MP2_SCHWARZ_CUTOFF  ||
					    fabs(Cap*Csq*Cbi*Crj*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

						// compute two-electron integral
						EE = contr_eri(
						     gto[p].nContract,
						     gto[p].exp, gto[p].coef, gto[p].norm,
						     gto[p].x0,        gto[p].y0,  gto[p].z0,
						     gto[p].l,         gto[p].m,   gto[p].n,
	
						     gto[q].nContract,
						     gto[q].exp, gto[q].coef, gto[q].norm,
						     gto[q].x0,        gto[q].y0,  gto[q].z0,
						     gto[q].l,         gto[q].m,   gto[q].n,
	
						     gto[i].nContract,
						     gto[i].exp, gto[i].coef, gto[i].norm,
						     gto[i].x0,        gto[i].y0,  gto[i].z0,
						     gto[i].l,         gto[i].m,   gto[i].n,
	
						     gto[j].nContract,
						     gto[j].exp, gto[j].coef, gto[j].norm,
						     gto[j].x0,        gto[j].y0,  gto[j].z0,
						     gto[j].l,         gto[j].m,   gto[j].n);
	
						// counter
						nEE++;
	
						arbs += EE*Cap*Crq*Cbi*Csj;
						asbr += EE*Cap*Csq*Cbi*Crj;
					}
				}
			}
		}

		contrib = arbs - asbr;
		sum += contrib*contrib/(eCorA[a]+eCorA[b]-eVirA[r]-eVirA[s]);
	}

	// compute mp2 energy: spin-dn contribution
	for(a=0; a < nCorrB; a++)
	for(b=a+1; b < nCorrB; b++)
	for(r=0; r < nVirB; r++)
	for(s=r+1; s < nVirB; s++){

		// reset values
		arbs = asbr = 0.0;

		for(p=0; p < nBasis; p++)
		for(q=0; q < nBasis; q++){
			Cap = CB[(a+nCore)*nBasis+p];
			Crq = CB[(r+  nEB)*nBasis+q];
			Csq = CB[(s+  nEB)*nBasis+q];
			Spq = schwarz[p*nBasis+q];
			if(fabs(Cap*Crq*Spq) > MP2_SCHWARZ_CUTOFF ||
			   fabs(Cap*Csq*Spq) > MP2_SCHWARZ_CUTOFF){

				for(i=0; i < nBasis; i++)
				for(j=0; j < nBasis; j++){
					Cbi = CB[(b+nCore)*nBasis+i];
					Csj = CB[(s+  nEB)*nBasis+j];
					Crj = CB[(r+  nEB)*nBasis+j];
					Sij = schwarz[i*nBasis+j];
					if( fabs(Cap*Crq*Cbi*Csj*Spq*Sij) > MP2_SCHWARZ_CUTOFF  ||
					    fabs(Cap*Csq*Cbi*Crj*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

						// compute two-electron integral
						EE = contr_eri(
						     gto[p].nContract,
						     gto[p].exp, gto[p].coef, gto[p].norm,
						     gto[p].x0,        gto[p].y0,  gto[p].z0,
						     gto[p].l,         gto[p].m,   gto[p].n,
	
						     gto[q].nContract,
						     gto[q].exp, gto[q].coef, gto[q].norm,
						     gto[q].x0,        gto[q].y0,  gto[q].z0,
						     gto[q].l,         gto[q].m,   gto[q].n,
	
						     gto[i].nContract,
						     gto[i].exp, gto[i].coef, gto[i].norm,
						     gto[i].x0,        gto[i].y0,  gto[i].z0,
						     gto[i].l,         gto[i].m,   gto[i].n,
	
						     gto[j].nContract,
						     gto[j].exp, gto[j].coef, gto[j].norm,
						     gto[j].x0,        gto[j].y0,  gto[j].z0,
						     gto[j].l,         gto[j].m,   gto[j].n);
	
						// counter
						nEE++;
	
						arbs += EE*Cap*Crq*Cbi*Csj;
						asbr += EE*Cap*Csq*Cbi*Crj;
					}
				}
			}
		}

		contrib = arbs - asbr;
		sum += contrib*contrib/(eCorB[a]+eCorB[b]-eVirB[r]-eVirB[s]);
	}

	// compute mp2 energy: cross-spin contribution
	for(a=0; a < nCorrA; a++)
	for(b=0; b < nCorrB; b++)
	for(r=0; r < nVirA; r++)
	for(s=0; s < nVirB; s++){

		// reset values
		arbs = 0.0;

		for(p=0; p < nBasis; p++)
		for(q=0; q < nBasis; q++){
			Cap = CA[(a+nCore)*nBasis+p];
			Crq = CA[(r+  nEA)*nBasis+q];
			Spq = schwarz[p*nBasis+q];
			if(fabs(Cap*Crq*Spq) > MP2_SCHWARZ_CUTOFF){

				for(i=0; i < nBasis; i++)
				for(j=0; j < nBasis; j++){
					Cbi = CB[(b+nCore)*nBasis+i];
					Csj = CB[(s+  nEB)*nBasis+j];
					Sij = schwarz[i*nBasis+j];
					if( fabs(Cap*Crq*Cbi*Csj*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

						// compute two-electron integral
						EE = contr_eri(
						     gto[p].nContract,
						     gto[p].exp, gto[p].coef, gto[p].norm,
						     gto[p].x0,        gto[p].y0,  gto[p].z0,
						     gto[p].l,         gto[p].m,   gto[p].n,
	
						     gto[q].nContract,
						     gto[q].exp, gto[q].coef, gto[q].norm,
						     gto[q].x0,        gto[q].y0,  gto[q].z0,
						     gto[q].l,         gto[q].m,   gto[q].n,
	
						     gto[i].nContract,
						     gto[i].exp, gto[i].coef, gto[i].norm,
						     gto[i].x0,        gto[i].y0,  gto[i].z0,
						     gto[i].l,         gto[i].m,   gto[i].n,
	
						     gto[j].nContract,
						     gto[j].exp, gto[j].coef, gto[j].norm,
						     gto[j].x0,        gto[j].y0,  gto[j].z0,
						     gto[j].l,         gto[j].m,   gto[j].n);
	
						// counter
						nEE++;
	
						arbs += EE*Cap*Crq*Cbi*Csj;
					}
				}
			}
		}

		contrib = arbs;
		sum += contrib*contrib/(eCorA[a]+eCorB[b]-eVirA[r]-eVirB[s]);
	}

	printf("The number of (pq|ij) integral computed:  %lld\n", nEE);

	printf("\nMP2 Energy Correction is %.6f Hartrees\n", sum);
	printf(
       "-------------------------------------------------------------\n");

	// free memory
	free(schwarz);

	return sum;
}


// mp2_uhf_semi_direct_aqij : computes MP2 correction energy for the case of 
// unrestricted orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for general case is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[a<b, r<s] { <ab||rs> / (ea+eb-er-es)  }
//
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. The integral is over spin-orbitals.
// Once the MP2 energy has been expressed explicitly in terms of spatial-orbital
// it assumes the form:
//
//   MP2 = SUM[a<b, r<s] {( (ar|bs) - (as|br) )^2 /(ea+eb-er-es) } 
//         a,b,r,s all spin up
//
//      +  SUM[a<b, r<s] {( (ar|bs) - (as|br) )^2 /(ea+eb-er-es) } 
//         a,b,r,s all spin dn
//
//      +  SUM[a b, r s] {( (ar|bs) )^2           /(ea+eb-er-es) } 
//         a,r are spin up and r,s are spin dn
//
// See research note page 45 for details
//
// The aqij variant will stores integral in memory 
// before calculating (ar|bs) and (as|br)
//
// Dec 22, 2014 - Teepanis Chachiyo
//    Initial implementation and testing. Adapted from mp2_rhf_direct
//
double mp2_uhf_semi_direct_aqij(
	int nBasis,                    // number of basis function
	int nEA,                       // number of occupied spin up orbitals
	int nEB,                       // number of occupied spin dn orbitals
	const double *eA,              // eigen values
	const double *eB,              // eigen values
	const double *CA,              // molecular orbitals
	const double *CB,              // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorrA,nCorrB;         // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int p,q,i,j;               // basis function index
	double EE;                 // compute electron integral
	int nVirA,nVirB;           // number of virtual orbitals
	const double *eVirA,*eVirB;// pointer to virtual orbital eigen value
	const double *eCorA,*eCorB;// correlated orbital eigen value
	double arbs, asbr;         // current integral value
	double *schwarz;           // schwarz inequality at basis level
	double *aqij;              // nBasis^3 intermediate storage
	double *aqijPtr;           // pointer to current aqbj 
	double Cap,Crq,Cbi,Csj,Csq,Crj; // transform coefficients
	double Spq,Sij;                 // schwarz screening
	double contrib, sum=0.0;        // acuumulated MP2 energy
	unsigned long long nEE=0;

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVirA  = nBasis - nEA;
	nCorrA = nEA    - nCore;
	nVirB  = nBasis - nEB;
	nCorrB = nEB    - nCore;

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----        UNRESTRICTED MP2 ENERGY CORRECTION         -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:      %4d %4d\n", nCore, nCore);
	printf("The number of included occupied orbitals:  %4d %4d\n", nCorrA, nCorrB);
	printf("The number of included virtual  orbitals:  %4d %4d\n", nVirA, nVirB);
	fflush(stdout);

	// locate virtual orbital and correlated orbitals eigen value
	eVirA = eA + nEA;
	eCorA = eA + nCore;
	eVirB = eB + nEB;
	eCorB = eB + nCore;

	// generate schwarz inequality between basis function
	schwarz = create_Schwarz(nBasis, gto);

	// allocate intermediate storage
	if((aqij=calloc(nBasis*nBasis*nBasis,sizeof(double)))==NULL){
		printf("mp2_uhf_semi_direct_aqij - error cannot allocate memory\n");
		exit(-1);
	}

	//
	// compute mp2 energy: spin-up contribution
	//
	for(a=0; a < nCorrA; a++){

		// reset (aqij)
		memset(aqij,0,nBasis*nBasis*nBasis*sizeof(double));

		// generate (aqij)
		for(p=0; p < nBasis; p++){

			// load coefficients
			Cap = CA[(a+nCore)*nBasis+p];

			// reset pointer
			aqijPtr = aqij;

			// screening
			if(fabs(Cap) > MP2_SCHWARZ_CUTOFF)
			for(q=0; q < nBasis; q++)
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){

				// load schwarz inequality values
				Spq = schwarz[p*nBasis+q];
				Sij = schwarz[i*nBasis+j];

				// screening
				if(fabs(Cap*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

					// compute two-electron integral
					EE = contr_eri(
					     gto[p].nContract,
					     gto[p].exp, gto[p].coef, gto[p].norm,
					     gto[p].x0,        gto[p].y0,  gto[p].z0,
					     gto[p].l,         gto[p].m,   gto[p].n,

					     gto[q].nContract,
					     gto[q].exp, gto[q].coef, gto[q].norm,
					     gto[q].x0,        gto[q].y0,  gto[q].z0,
					     gto[q].l,         gto[q].m,   gto[q].n,

					     gto[i].nContract,
					     gto[i].exp, gto[i].coef, gto[i].norm,
					     gto[i].x0,        gto[i].y0,  gto[i].z0,
					     gto[i].l,         gto[i].m,   gto[i].n,

					     gto[j].nContract,
					     gto[j].exp, gto[j].coef, gto[j].norm,
					     gto[j].x0,        gto[j].y0,  gto[j].z0,
					     gto[j].l,         gto[j].m,   gto[j].n);

					// counter
					nEE++;

					*aqijPtr += EE*Cap;
				}

				// increment pointer index
				aqijPtr++;
			}
		}

		// use (aq|ij) to generate (ar|bs) and (as|br) for alpha-alpha case
		for(b=a+1; b < nCorrA; b++)
		for(r=0; r < nVirA; r++)
		for(s=r+1; s < nVirA; s++){

			//reset value
			arbs = asbr = 0.0;

			// reset pointer 
			aqijPtr = aqij;

			for(q=0; q < nBasis; q++)
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				// load coefficient
				Cbi = CA[(b+nCore)*nBasis+i];
				Crq = CA[(r+  nEA)*nBasis+q];
				Csq = CA[(s+  nEA)*nBasis+q];
				Csj = CA[(s+  nEA)*nBasis+j];
				Crj = CA[(r+  nEA)*nBasis+j];

				// compute integral
				arbs += (*aqijPtr)*Crq*Cbi*Csj;
				asbr += (*aqijPtr)*Csq*Cbi*Crj;

				// step up pointer
				aqijPtr++;
			}

			// compute MP2 contributions
			contrib = arbs - asbr;
			sum += contrib*contrib/(eCorA[a]+eCorA[b]-eVirA[r]-eVirA[s]);
		}


		// use (aq|ij) to generate (ar|bs) for cross-spin case
		for(b=0; b < nCorrB; b++)
		for(r=0; r < nVirA; r++)
		for(s=0; s < nVirB; s++){

			//reset value
			arbs = 0.0;

			// reset pointer 
			aqijPtr = aqij;

			for(q=0; q < nBasis; q++)
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				// load coefficient
				Cbi = CB[(b+nCore)*nBasis+i];
				Crq = CA[(r+  nEA)*nBasis+q];
				Csj = CB[(s+  nEB)*nBasis+j];

				// compute integral
				arbs += (*aqijPtr)*Crq*Cbi*Csj;

				// step up pointer
				aqijPtr++;
			}

			// compute MP2 contributions
			contrib = arbs;
			sum += contrib*contrib/(eCorA[a]+eCorB[b]-eVirA[r]-eVirB[s]);
		}

	}


	//
	// compute mp2 energy: spin-dn contribution
	//
	for(a=0; a < nCorrB; a++){

		// reset (aqij)
		memset(aqij,0,nBasis*nBasis*nBasis*sizeof(double));

		// generate (aqij)
		for(p=0; p < nBasis; p++){

			// load coefficients
			Cap = CB[(a+nCore)*nBasis+p];

			// reset pointer
			aqijPtr = aqij;

			// screening
			if(fabs(Cap) > MP2_SCHWARZ_CUTOFF)
			for(q=0; q < nBasis; q++)
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){

				// load schwarz inequality values
				Spq = schwarz[p*nBasis+q];
				Sij = schwarz[i*nBasis+j];

				// screening
				if(fabs(Cap*Spq*Sij) > MP2_SCHWARZ_CUTOFF){

					// compute two-electron integral
					EE = contr_eri(
					     gto[p].nContract,
					     gto[p].exp, gto[p].coef, gto[p].norm,
					     gto[p].x0,        gto[p].y0,  gto[p].z0,
					     gto[p].l,         gto[p].m,   gto[p].n,

					     gto[q].nContract,
					     gto[q].exp, gto[q].coef, gto[q].norm,
					     gto[q].x0,        gto[q].y0,  gto[q].z0,
					     gto[q].l,         gto[q].m,   gto[q].n,

					     gto[i].nContract,
					     gto[i].exp, gto[i].coef, gto[i].norm,
					     gto[i].x0,        gto[i].y0,  gto[i].z0,
					     gto[i].l,         gto[i].m,   gto[i].n,

					     gto[j].nContract,
					     gto[j].exp, gto[j].coef, gto[j].norm,
					     gto[j].x0,        gto[j].y0,  gto[j].z0,
					     gto[j].l,         gto[j].m,   gto[j].n);

					// counter
					nEE++;

					*aqijPtr += EE*Cap;
				}

				// increment pointer index
				aqijPtr++;
			}
		}

		// use (aq|ij) to generate (ar|bs) and (as|br) for beta-beta case
		for(b=a+1; b < nCorrB; b++)
		for(r=0; r < nVirB; r++)
		for(s=r+1; s < nVirB; s++){

			//reset value
			arbs = asbr = 0.0;

			// reset pointer 
			aqijPtr = aqij;

			for(q=0; q < nBasis; q++)
			for(i=0; i < nBasis; i++)
			for(j=0; j < nBasis; j++){
				// load coefficient
				Cbi = CB[(b+nCore)*nBasis+i];
				Crq = CB[(r+  nEB)*nBasis+q];
				Csq = CB[(s+  nEB)*nBasis+q];
				Csj = CB[(s+  nEB)*nBasis+j];
				Crj = CB[(r+  nEB)*nBasis+j];

				// compute integral
				arbs += (*aqijPtr)*Crq*Cbi*Csj;
				asbr += (*aqijPtr)*Csq*Cbi*Crj;

				// step up pointer
				aqijPtr++;
			}

			// compute MP2 contributions
			contrib = arbs - asbr;
			sum += contrib*contrib/(eCorB[a]+eCorB[b]-eVirB[r]-eVirB[s]);
		}

	}

	printf("The number of (pq|ij) integral computed:  %lld\n", nEE);

	printf("\nMP2 Energy Correction is %.6f Hartrees\n", sum);
	printf(
       "-------------------------------------------------------------\n");

	// free memory
	free(aqij);
	free(schwarz);

	return sum;
}


// mp2_uhf_aqij : computes MP2 correction energy for the case of unrestricted 
// orbitals. It returns the correction energy. Not all orbitals are 
// included into correlations. The core electrons according to "Noble Gas 
// Core" are excluded.
//
// The MP2 correlation energy for general case is defined as
// (Szabo, "Modern Quantum Chemistry", page 352, equation 6.74)
//
//    
//   MP2 = SUM[a<b, r<s] { <ab||rs> / (ea+eb-er-es)  }
//
//   
// Here, the indexes a and b are occupied orbitals; whereas the indexes
// r and r are virtual orbitals. The integral is over spin-orbitals.
// Once the MP2 energy has been expressed explicitly in terms of spatial-orbital
// it assumes the form:
//
//   MP2 = SUM[a<b, r<s] {( (ar|bs) - (as|br) )^2 /(ea+eb-er-es) } 
//         a,b,r,s all spin up
//
//      +  SUM[a<b, r<s] {( (ar|bs) - (as|br) )^2 /(ea+eb-er-es) } 
//         a,b,r,s all spin dn
//
//      +  SUM[a b, r s] {( (ar|bs) )^2           /(ea+eb-er-es) } 
//         a,r are spin up and r,s are spin dn
//
// See research note page 45 for details
//
//
// Dec 23, 2014 - Teepanis Chachiyo
//    Initial implementation and testing. Adapted from mp2_rhf_aqij
//
// Sep 17, 2014 - Teepanis Chachiyo
//    maxMem check
//
double mp2_uhf_aqij(
	int nBasis,                    // number of basis function
	int nEA,                       // number of occupied spin up orbitals
	int nEB,                       // number of occupied spin dn orbitals
	const double *eA,              // eigen values
	const double *eB,              // eigen values
	const double *CA,              // molecular orbitals
	const double *CB,              // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	int nCore;                 // number of core orbital to ignore
	int nCorrA,nCorrB;         // number of correlated orbitals
	int a,b,r,s;               // molecular orbitals
	int q,i,j;                 // basis function index
	double *aqbj, *arbj;       // intermediate integrals
	double Cbi, arbjSum;       // intermediate values
	double arbs, asbr;         // integral needed to compute MP2 energy
	double *ptr;               // pointer to current aqbj value
	int nVirA,nVirB;           // number of virtual orbitals
	const double *eVirA,*eVirB;// pointer to virtual orbital eigen value
	const double *eCorA,*eCorB;// correlated orbital eigen value
	double sum;                // accumulated mp2 correction
	double contrib;            // current contribution to sum
	double *aqijBank;          // 2e storage between orbitals and basis functions
	double **IJm;              // IJm mapping matrix
	double *Schwarz;           // schwarz inequality
	int nEE;                   // number of integral index
	int maxCorrA,maxCorrB;     // maximum correlated orbitals fits in memory
	int aCycle;                // a index for this round
	int aCount;                // counter for a index round

	// there should be a way to detect possible singularity (or 
	// degeneracy) which might happen when ea+eb-er-es = zero
	// Teepanis (Dec 23, 2014)

	// compute the number of frozen core
	nCore = getFrozenCore_RHF(mol);

	// compute number of virtual orbitals and correlated orbitals
	nVirA  = nBasis - nEA;
	nCorrA = nEA    - nCore;
	nVirB  = nBasis - nEB;
	nCorrB = nEB    - nCore;

	// report
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----        UNRESTRICTED MP2 ENERGY CORRECTION         -----\n"
	"-------------------------------------------------------------\n");

	printf("The number of excluded core orbitals:      %4d %4d\n", nCore, nCore);
	printf("The number of included occupied orbitals:  %4d %4d\n", nCorrA, nCorrB);
	printf("The number of included virtual  orbitals:  %4d %4d\n", nVirA, nVirB);
	fflush(stdout);

	// locate virtual orbital and correlated orbitals eigen value
	eVirA = eA + nEA;
	eCorA = eA + nCore;
	eVirB = eB + nEB;
	eCorB = eB + nCore;

	// memory requirement
	Schwarz = create_Schwarz(nBasis, gto);
	for(nEE=0,i=0; i < nBasis; i++)
	for(j=0; j <= i; j++)
		if(Schwarz[i*nBasis+j] >= MP2_SCHWARZ_CUTOFF) nEE++;

	// allocate computed integral
	if((aqbj=calloc(nBasis*nBasis,sizeof(double )))==NULL ||
	   (arbj=calloc(nBasis*nBasis,sizeof(double )))==NULL ||
	   (IJm =calloc(nBasis*nBasis,sizeof(double*)))==NULL){
		printf("mp2_uhf_aqij - error cannot allocate aqbj,arbj,IJm\n");
		exit(-1);
	}

	// Check the case where opt->maxMem is too low to do even one round
	// Teepanis - 17 September 2014
	if(opt->maxMem*1.0E6 < nEE*nBasis*sizeof(double)){
		printf("mp2_uhf_aqij - error maxMem must be at least %.1f MB\n",
		      nEE*nBasis*sizeof(double)*1.0E-6);
		exit(-1);
	}

	maxCorrA = (int)floor(opt->maxMem*1.0E6/nEE/nBasis/sizeof(double));
	if(maxCorrA > nCorrA) maxCorrA = nCorrA;

	maxCorrB = (int)floor(opt->maxMem*1.0E6/nEE/nBasis/sizeof(double));
	if(maxCorrB > nCorrB) maxCorrB = nCorrB;

	sum = nEE*nBasis*sizeof(double)*
	      ( (maxCorrA>maxCorrB)? maxCorrA : maxCorrB );
	mp2_report_memory(sum);

	// report the number of cycles
	printf("The number of cycles to compute (pq|ij):   %4d %4d\n",
	          nCorrA==0?0: (int)ceil((float)nCorrA/maxCorrA) 
	         ,nCorrB==0?0: (int)ceil((float)nCorrB/maxCorrB));
	fflush(stdout);

	sum = 0.0;

	//////////////////////////
	// spin-up contribution //
	//////////////////////////

	// compute MP2 energy in multiple cycles
	for(a=0,aCount=0; a < nCorrA; a+=maxCorrA,aCount++){

		// compute number of orbitals to compute for this cycles
		if((nCorrA-a) < maxCorrA) maxCorrA = (nCorrA-a);

		// generate aqbj integral
		aqijBank = mp2_gen_aqij_Quartet(nBasis, maxCorrA, nCore+a, CA, gto, Schwarz);

		// build IJm mapping matrix
		for(nEE=0,i=0; i < nBasis; i++)
		for(j=0; j <= i; j++){
			if(Schwarz[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) IJm[i*nBasis+j] = NULL;
			else{
				IJm[i*nBasis+j] = aqijBank + nEE*nBasis*maxCorrA;
				nEE++;
			}
			IJm[j*nBasis+i] = IJm[i*nBasis+j];
		}

		for(aCycle=0; aCycle < maxCorrA; aCycle++){

			//
			// alpha-alpha case
			//
			for(b=aCycle+a+1; b < nCorrA; b++){
		
				memset(aqbj,0,sizeof(double)*nBasis*nBasis);
				for(i=0; i < nBasis; i++)
				for(j=0; j < nBasis; j++)
				if((ptr=IJm[i*nBasis+j])){
					Cbi = CA[(b+nCore)*nBasis+i];
					if(fabs(Schwarz[i*nBasis+j]*Cbi) > MP2_SCHWARZ_CUTOFF)
					for(q=0; q < nBasis; q++)
						aqbj[j*nBasis+q] += Cbi*ptr[aCycle*nBasis+q];
				}
	
				memset(arbj,0,sizeof(double)*nBasis*nBasis);
				for(r=0; r < nVirA; r++)
				for(j=0; j < nBasis; j++){
					for(arbjSum=0.0,q=0; q < nBasis; q++)
						arbjSum += CA[(r+nEA)*nBasis+q]*aqbj[j*nBasis+q];
					arbj[r*nBasis+j] = arbjSum;
				}
	
				for(r=0; r < nVirA; r++)
				for(s=r+1; s < nVirA; s++){
					arbs = asbr = 0.0;
					for(j=0; j < nBasis; j++){
						arbs += CA[(s+nEA)*nBasis+j]*arbj[r*nBasis+j];
						asbr += CA[(r+nEA)*nBasis+j]*arbj[s*nBasis+j];
					}
	
					// compute MP2 contributions	
					contrib = arbs - asbr;
					sum += contrib*contrib/(eCorA[aCycle+a]+eCorA[b]-eVirA[r]-eVirA[s]);
	
				}
			}

			//
			// alpha-beta case
			//
			for(b=0; b < nCorrB; b++){
		
				memset(aqbj,0,sizeof(double)*nBasis*nBasis);
				for(i=0; i < nBasis; i++)
				for(j=0; j < nBasis; j++)
				if((ptr=IJm[i*nBasis+j])){
					Cbi = CB[(b+nCore)*nBasis+i];
					if(fabs(Schwarz[i*nBasis+j]*Cbi) > MP2_SCHWARZ_CUTOFF)
					for(q=0; q < nBasis; q++)
						aqbj[j*nBasis+q] += Cbi*ptr[aCycle*nBasis+q];
				}
	
				memset(arbj,0,sizeof(double)*nBasis*nBasis);
				for(r=0; r < nVirA; r++)
				for(j=0; j < nBasis; j++){
					for(arbjSum=0.0,q=0; q < nBasis; q++)
						arbjSum += CA[(r+nEA)*nBasis+q]*aqbj[j*nBasis+q];
					arbj[r*nBasis+j] = arbjSum;
				}
	
				for(r=0; r < nVirA; r++)
				for(s=0; s < nVirB; s++){
					arbs = 0.0;
					for(j=0; j < nBasis; j++)
						arbs += CB[(s+nEB)*nBasis+j]*arbj[r*nBasis+j];

					// compute MP2 contributions
					contrib = arbs;
					sum += contrib*contrib/(eCorA[aCycle+a]+eCorB[b]-eVirA[r]-eVirB[s]);
				}
			}

		}

		// free memory
		free(aqijBank);

		// progress report
		printf("The number of cycles already finished:     %4d\n",aCount+1);
		fflush(stdout);

	}


	//////////////////////////
	// spin-dn contribution //
	//////////////////////////

	// compute MP2 energy in multiple cycles
	for(a=0,aCount=0; a < nCorrB; a+=maxCorrB,aCount++){

		// compute number of orbitals to compute for this cycles
		if((nCorrB-a) < maxCorrB) maxCorrB = (nCorrB-a);

		// generate aqbj integral
		aqijBank = mp2_gen_aqij_Quartet(nBasis, maxCorrB, nCore+a, CB, gto, Schwarz);

		// build IJm mapping matrix
		for(nEE=0,i=0; i < nBasis; i++)
		for(j=0; j <= i; j++){
			if(Schwarz[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) IJm[i*nBasis+j] = NULL;
			else{
				IJm[i*nBasis+j] = aqijBank + nEE*nBasis*maxCorrB;
				nEE++;
			}
			IJm[j*nBasis+i] = IJm[i*nBasis+j];
		}

		for(aCycle=0; aCycle < maxCorrB; aCycle++){

			//
			// beta-beta case
			//
			for(b=aCycle+a+1; b < nCorrB; b++){
		
				memset(aqbj,0,sizeof(double)*nBasis*nBasis);
				for(i=0; i < nBasis; i++)
				for(j=0; j < nBasis; j++)
				if((ptr=IJm[i*nBasis+j])){
					Cbi = CB[(b+nCore)*nBasis+i];
					if(fabs(Schwarz[i*nBasis+j]*Cbi) > MP2_SCHWARZ_CUTOFF)
					for(q=0; q < nBasis; q++)
						aqbj[j*nBasis+q] += Cbi*ptr[aCycle*nBasis+q];
				}
	
				memset(arbj,0,sizeof(double)*nBasis*nBasis);
				for(r=0; r < nVirB; r++)
				for(j=0; j < nBasis; j++){
					for(arbjSum=0.0,q=0; q < nBasis; q++)
						arbjSum += CB[(r+nEB)*nBasis+q]*aqbj[j*nBasis+q];
					arbj[r*nBasis+j] = arbjSum;
				}
	
				for(r=0; r < nVirB; r++)
				for(s=r+1; s < nVirB; s++){
					arbs = asbr = 0.0;
					for(j=0; j < nBasis; j++){
						arbs += CB[(s+nEB)*nBasis+j]*arbj[r*nBasis+j];
						asbr += CB[(r+nEB)*nBasis+j]*arbj[s*nBasis+j];
					}
	
					// compute MP2 contributions	
					contrib = arbs - asbr;
					sum += contrib*contrib/(eCorB[aCycle+a]+eCorB[b]-eVirB[r]-eVirB[s]);
	
				}
			}

		}

		// free memory
		free(aqijBank);

		// progress report
		printf("The number of cycles already finished:     %4d\n",aCount+1);
		fflush(stdout);

	}

	printf("\n"
	       "MP2 Energy Correction is %12.6f Hartree\n"
	       "                         %12.6f eV\n", sum, HARTREE2EV*sum);
	printf(
       "-------------------------------------------------------------\n");
	fflush(stdout);

	// free memory
	mp2_gen_aqij_Quartet(0, 0, 0, NULL, NULL, NULL);
	free(aqbj);
	free(arbj);
	free(IJm);
	free(Schwarz);

	return sum;
}
