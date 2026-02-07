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
#ifndef DFT_H
#define DFT_H

#include "basis.h"
#include "grid.h"


// set basis function value to zero if below this value
#define BGRID_CUTOFF 1.0E-6

// set rho to zero if below this value
#define RHO_CUTOFF 1.0E-10 

struct BGrid_t * genBasisGrid(
const int nBasis,               // number of basis function
const struct GTOBasis_t *gto,   // basis function data
const struct MolGrid_t *grid,   // grid point data
const double cutoff);           // cutoff value

struct BGrid_t * cleanBGrid(struct BGrid_t *bg);

void getRho(
const int nBasis,              // number of basis function
const double *PA,              // spin up density matrix
const double *PB,              // spin dn density matrix
const struct GTOBasis_t *gto,  // basis function info
const struct MolGrid_t *grid,  // grid point data
const double cutoff,           // basis function cutoff
      double *rhoa,            // returned spin up density 
      double *rhob);           // returned spin dn density

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
      double *Fz);             // returned (incremental) force

struct gridBasisCenter_t{
	int *gCenter;      // nearest basis center index
	double *bR2Max;    // maximum r2 range
	double *bR2;       // basis distance matrix
	int *bAtomID;      // basis atomID
};

struct gridBasisCenter_t * cleanGridBasisCenter(
int nBasis,
struct gridBasisCenter_t *gbCenter);

struct gridBasisCenter_t * genGridBasisCenter(
int nBasis,
const struct GTOBasis_t *gto,
const struct Molecule_t *mol,
const struct MolGrid_t *grid);

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
      double *Exc,             // returned (incremental) energy grid
      double *XA,              // returned spin up matrix
      double *XB,              // returned spin dn matrix
      double *Fx,              // returned (incremental) force
      double *Fy,              // returned (incremental) force
      double *Fz);             // returned (incremental) force

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
      double *Exc,             // returned (incremental) energy grid
      double *XA,              // returned spin up matrix
      double *XB,              // returned spin dn matrix
      double *Fx,              // returned (incremental) force
      double *Fy,              // returned (incremental) force
      double *Fz);             // returned (incremental) force

void getDFT_xSlater(
double rhoa,     // spin up electron density
double rhob,     // spin dn electron density
double *Exc,     // returned (incremental) exchange-correlation energy
double *dEdrhoa, // returned (incremental) spin up potential
double *dEdrhob);// returned (incremental) spin dn potential

extern double VWN5_A[3];
extern double VWN5_x[3];
extern double VWN5_b[3];
extern double VWN5_c[3];
void getDFT_cVWN5(
double rhoa,     // spin up electron density
double rhob,     // spin dn electron density
double *Exc,     // returned (incremental) exchange-correlation energy
double *dEdrhoa, // returned (incremental) spin-up potential
double *dEdrhob);// returned (incremental) spin-dn potential

void getDFT_cChachiyo(
double rhoa,     // spin up   electron density
double rhob,     // spin down electron density
double *Ec,      // returned (incremental)  correlation energy
double *dEdrhoa, // returned (incremental) spin up   potential
double *dEdrhob);// returned (incremental) spin down potential


#define GGAxSLATER   0  // Dirac LDA exchange
#define GGAxBECKE88  1  // Becke 1988
#define GGAxPY86     2  // Perdew-Yue 1986
#define GGAxPBE      3  // Perdew-Burke-Ernzerhof
#define GGAxCHACHIYO 4  // Chachiyo first principle exchange
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
double *Gb);     // returned (incremental) coef of spin dn scaled gradient

void getDFT_xMVS(
double rhoa,     // spin up electron density
double rhob,     // spin dn electron density
double rhoag,    // magnitude of the gradient of spin up electron density
double rhobg,    // magnitude of the gradient of spin dn electron density
double taua,     // meta-GGA tau for spin-up density
double taub,     // meta-GGA tau for spin-dn density
double *Ex);     // returned (incremental) exchange energy

void getDFT_xGGA_Chachiyo(
double rhoa,     // spin up   electron density
double rhob,     // spin down electron density
double rhoag,    // spin up   density gradient 
double rhobg,    // spin down density gradient
double *Ex,      // returned (incremental) exchange energy
double *dEdrhoa, // returned (incremental) spin up   potential
double *dEdrhob, // returned (incremental) spin down potential
double *Ga,      // returned (incremental) coef. of alpha density gradient
double *Gb);     // returned (incremental) coef. of beta  density gradient

void getDFT_cGGA_Chachiyo(
double rhoa,     // spin up   electron density
double rhob,     // spin down electron density
double rhog,     // total gradient
double *Ec,      // returned (incremental) correlation energy
double *dEdrhoa, // returned (incremental) spin up   potential
double *dEdrhob, // returned (incremental) spin down potential
double *G);      // returned (incremental) coef. of total density gradient
#endif
