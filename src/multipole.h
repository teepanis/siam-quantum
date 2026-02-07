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

#ifndef MULTIPOLE_H
#define MULTIPOLE_H

#include "basis.h"

#define MULTIPOLE_CHARGE_CUTOFF 1.0E-18
#define MULTIPOLE_RADII2_CUTOFF 1.0E-00

// multipole expansion of a product between two basis functions
struct Multipole_t{
	double q;                    // total charge
	double x,y,z;                // center
};

struct Multipole_t *genMultipole(
	int nBasis,                     // the number of basis functions
	const struct GTOBasis_t *gto);  // pointer to basis function information

double coulombMultipole(
	int nBasis,                       // number of basis function
	int p, int q,                     // index to first multipole
	int i, int j,                     // index to second multipole
	const struct Multipole_t *mpole); // pointer to multipole info

#endif // MULTIPOLE_H
