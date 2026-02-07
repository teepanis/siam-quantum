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

#ifndef MATRIX_H
#define MATRIX_H

#include "basis.h"
#include "mol.h"
#include "int.h"
#include "option.h"

// Schwarz inequality cut-off for 2e integral
#define SCHWARZ_CUTOFF 1.0E-10
#define MP2_SCHWARZ_CUTOFF 1.0E-8

double GTO_overlap(int i,                        // ith basis
                   int j,                        // jth basis
                   const struct GTOBasis_t *gto);// basis database

double GTO_moment(int i,                            // ith basis
                  int j,                            // jth basis
                  const struct GTOBasis_t *gto,     // basis database
                  int mx, int my, int mz,           // moment order
                  double xc, double yc, double zc); // moment center

double GTO_kinetic(int i,                        // ith basis
                   int j,                        // jth basis
                   const struct GTOBasis_t *gto);// basis database

double GTO_nuclei(int i,                        // ith basis
                  int j,                        // jth basis
                  const struct GTOBasis_t *gto, // basis database
                  const struct Molecule_t *mol);// molecule database

double GTO_nai(int i,                            // index of GTO_tho
               int j,                            // index of GTO_tho
               const struct GTOBasis_t *gto,     // basis function info.
               double xc, double yc, double zc); // position to find potential

void GTO_efi(int i,                               // ith basis
             int j,                               // jth basis
             const struct GTOBasis_t *gto,        // basis database
             double xc, double yc, double zc,     // point to evaluate efi
             double *ex, double *ey, double *ez); // returned value

double contr_eri(
             int lena,double *aexps,double *acoefs,double *anorms,
             double xa,double ya,double za,int la,int ma,int na,
             int lenb,double *bexps,double *bcoefs,double *bnorms,
             double xb,double yb,double zb,int lb,int mb,int nb,
             int lenc,double *cexps,double *ccoefs,double *cnorms,
             double xc,double yc,double zc,int lc,int mc,int nc,
             int lend,double *dexps,double *dcoefs,double *dnorms,
             double xd,double yd,double zd,int ld,int md,int nd);

double * create_Schwarz(
	int nBasis,                     // number of basis functions
	const struct GTOBasis_t *gto);  // basis set info

//unsigned long long get_nEE(int nBasis, double *Schwarz);

void GTO_2JK_Matrix(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G);                  // return matrix

void GTO_2JK_Matrix_NoSymm(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G);                  // return matrix

void GTO_2JK_Matrix_Symm(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G);                  // return matrix

void GTO_2JK_Matrix_Symm_Shell(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G);                  // return matrix

void GTO_2JK_Matrix_Symm_Shell(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G);                  // return matrix

void GTO_2JK_Matrix_ShellSet(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G);

int shell_prep(
	int nBasis,                    // number of basis function
	const struct GTOBasis_t *gto,  // basis function data structure
	const double *Schwarz,         // Schwarz matrix for screening
	int **RshellMap,               // returned shellMap
	int **RshellMaxL,              // returned shellMaxL
	double **RSchwarz_Shell);      // returned max Schwarz within shell

void GTO_JK_Matrix_ShellSet(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB,                    // return G for spin down
	struct option_t *opt);         // global option

/*
void GTO_JK_Matrix_Quartet(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *schwarz_basis,   // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB,                    // return G for spin down
	struct option_t *opt);         // global option
*/

void GTO_JK_Matrix_Quartet_Parallel(
	int childID,                   // child id number
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *schwarz_basis,   // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *JT,                    // returned total Coulomb matrix
	double *KA,                    // returned Exchange for spin up
	double *KB,                    // returned Exchange for spin down
	struct option_t *opt);         // global option

void GTO_JK_Matrix_PrimeSpace(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double cutoff,                 // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB);                   // return G for spin down

void GTO_JK_Matrix_PrimeSpaceShellSet(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double cutoff,                 // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB);                   // return G for spin down

double * mp2_gen_aqbj_ShellSet(
	int nBasis,                     // number of basis function
	int nOcc,                       // number of occpupied orbitals
	int nCore,                      // number of core orbitals to ignore
	const double *C,                // molecular orbitals
	const struct GTOBasis_t *gto);  // basis function structure

double * mp2_gen_aqij_ShellSet(
	int nBasis,                     // number of basis function
	int nCorr,                      // number of correlated orbitals
	int nCore,                      // number of core orbitals to ignore
	const double *C,                // molecular orbitals
	const struct GTOBasis_t *gto);  // basis function structure

double * mp2_gen_aqij_Quartet(
	int nBasis,                     // number of basis function
	int nCorr,                      // number of correlated orbitals
	int nCore,                      // number of core orbitals to ignore
	const double *C,                // molecular orbitals
	const struct GTOBasis_t *gto,   // basis function structure
	const double *schwarz_basis);   // schwarz matrix at basis level

void GTO_JK_Matrix_NoSymm(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double cutoff,                 // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB);                   // return G for spin down

int isSameShell(
	int nBasis,                    // number of basis functions
	int i, int j,                  // which basis to test
	const struct GTOBasis_t *gto); // basis set info

#define MAXBASIS    16
#define MAXCONTRACT 32
struct GTOShell_t{
	double x,y,z;                             // shell center
	int nBasis;                               // number of basis in the shell
	int nContract;                            // number contract functions
	int maxL;                                 // maximum angular index
	int min,max;                              // basis index
	int l[MAXBASIS],m[MAXBASIS],n[MAXBASIS];  // angular index
	double exps[MAXCONTRACT];                 // exponents
	double coef[MAXCONTRACT][MAXBASIS];       // coefficients
	int type;                                 // predefined type 
#define TYPES 0
#define TYPEL 1
#define TYPEP 2
#define TYPED 3
#define TYPEF 4
#define TYPEG 5
#define TYPEX 6
};
#endif
