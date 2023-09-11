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
#ifndef HESSIAN_H
#define HESSIAN_H

#include "basis.h"
#include "mol.h"
#include "option.h"

// numerical step in Bohr
#define HESSIAN_STEPSIZE        0.008
// 0.005 causes error for the lowest frequency of CH3 HF/6-31G*
//              the analytical freq is 308, but numerical one is 314
//              setting it 0.01 brings the numerical freq to 308
//
// 0.01  causes non-degeneraacy of some modes in CH4/6-31G* at
//              frequency around 1703 >> 1703 and 1702
//              setting it to 0.008, restores the degeneracy
//
// 0.008 seems to be the sweet spot
//
// Teepanis Chachiyo 27 March, 2020

 // CRC Handbook of Chemistry and Physics 85th Edition
#define HARTREE2J 4.35974417E-18  // conversion hartree to joule
#define AMU2KG    1.66053886E-27  // conversion u to kilogram
#define AMU2ME    1.82288848E+3   // conversion u to atomic unit
#define ME_KG     9.10938262E-31  // electron mass in kg
#define MHZ2CM    3.33564E-5      // conversion MHz to cm^-1
#define CM2HARTREE 4.556335E-6    // conversion cm^-1 to hartree
#define CM2KCALMOL 2.85914E-3     // conversion cm^-1 to kcal/mol

double Z2_amu_mass(int Z, int *isotope);

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
	struct option_t *opt);   // options

#endif
