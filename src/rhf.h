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

#ifndef RHF_H
#define RHF_H
#include <math.h>
#include "mol.h"
#include "basis.h"
#include "option.h"

double eval_chi(int i, struct GTOBasis_t *gto,
                double x, double y, double z);

void rhf_getDMatrix(int nBasis, int nOcc, double *C, double *P);

double rhf_rho(int nBasis,              // number of basis function
               struct GTOBasis_t *gto,  // function structure
               double *P,               // density matrix
               double x, double y, double z);

double rhf_mo(int nBasis,              // number of basis function
              struct GTOBasis_t *gto,  // function structure
              double *C,               // molecular orbital
              int n,                   // orbital index (1st orbital is zero)
              double x, double y, double z);

double rhf(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nE,                  // total number of electrons
	double *C,               // returned molecular orbitals
	double *e,               // returned eigen values
	struct option_t *opt);   // global option
#endif
