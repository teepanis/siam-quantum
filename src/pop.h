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
#ifndef POP_H
#define POP_H
void mulliken(
	int nBasis,               // number of basis function
	struct GTOBasis_t * gto,  // pointer to basis set information
	struct Molecule_t * mol,  // pointer to molecular structure
	int nEA,                  // number of spin up electron
	int nEB,                  // number of spin down electron
	double *CA,               // calculated molecular orbital coeff.
	double *CB,               // save but for beta spin
	double *eA,               // spin up eigen value
	double *eB,               // spin down eigen value
	struct option_t *option); // print according to options

void electrostatic(
	int nBasis,               // number of basis function
	struct GTOBasis_t * gto,  // pointer to basis set information
	struct Molecule_t * mol,  // pointer to molecular structure
	int nEA,                  // number of spin up electron
	int nEB,                  // number of spin down electron
	double *CA,               // calculated molecular orbital coeff.
	double *CB,               // save but for beta spin
	double *eA,               // spin up eigen value
	double *eB,               // spin down eigen value
	struct option_t *option); // print according to options

void electric_multipole(
	int nBasis,               // number of basis function
	struct GTOBasis_t * gto,  // pointer to basis set information
	struct Molecule_t * mol,  // pointer to molecular structure
	int nEA,                  // number of spin up electron
	int nEB,                  // number of spin down electron
	double *CA,               // calculated molecular orbital coeff.
	double *CB,               // save but for beta spin
	double *eA,               // spin up eigen value
	double *eB,               // spin down eigen value
	struct option_t *option); // print according to options
#endif
