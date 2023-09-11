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

#ifndef CHECK_H
#define CHECK_H

#include "basis.h"
#include "option.h"
#include "mol.h"

double load_checkpoint_orbital(int nBasis,                  // number of basis function
                               double *CA,                  // alpha orbitals
                               double *CB,                  // beta orbitals 
                               double *eA,                  // alpha eigen values
                               double *eB,                  // beta eigen values
                               const struct option_t *opt); // global options

void guess_checkpoint_orbital(
	int nBasis,                   // number of basis function
	const struct GTOBasis_t *gto, // basis function
    const struct Molecule_t *mol, // molecular structure structure
	const struct option_t *opt,   // global options
	const double *S,              // overlap matrix 
	double *CA,                   // returned guess alpha orbital
	double *CB);                  // returned guess beta  orbital

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
                     struct option_t *opt ); // global options

struct Molecule_t *load_checkpoint_molecule(const struct option_t *opt);
int load_checkpoint_multiplicity(const struct option_t *opt);
struct GTOBasis_t *load_checkpoint_basis(const struct option_t *opt, int *nBasis);

#endif
