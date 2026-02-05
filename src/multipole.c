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

#include "multipole.h"
#include "basis.h"
#include "matrix.h"


//
// genMultipole : generate multipole moments for a product of basis function.
// This can be used to approximate coulomb integral.
//
// 23 Dec, 2012 - Teepanis Chachiyo
//       Initial implementation and testing
//
struct Multipole_t *genMultipole(
	int nBasis,                     // the number of basis functions
	const struct GTOBasis_t *gto){  // pointer to basis function information

	struct Multipole_t *mpole;  // pointer to multipole data structure
	double q;                   // total charge
	double px,py,pz;            // dipole moment
	double x,y,z;               // center
	int i,j;                    // basis function index

	// allocate memory
	mpole = (struct Multipole_t *)
	        calloc(nBasis*nBasis, sizeof(struct Multipole_t));
	if(mpole==NULL){
		printf("genMultipole: Error - cannot allocate memory\n");
		exit(-1);
	}

	// loop thru all basis functions
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){

		// compute total charge
		q = GTO_overlap(i,j,gto);

		// truncate total charge
		if(fabs(q) < MULTIPOLE_CHARGE_CUTOFF) q = 0.0;

		// compute charge center
		x = (gto[i].x0 + gto[j].x0)/2.0;
		y = (gto[i].y0 + gto[j].y0)/2.0;
		z = (gto[i].z0 + gto[j].z0)/2.0;

		// compute dipole moment
		px = GTO_moment(i,j, gto, 1,0,0, x,y,z);
		py = GTO_moment(i,j, gto, 0,1,0, x,y,z);
		pz = GTO_moment(i,j, gto, 0,0,1, x,y,z);

		// re-evaluate new center so that dipole becomes zero
		if(q!=0.0){
			x = x + px/q;
			y = y + py/q;
			z = z + pz/q;
		}

		// store values
		mpole[i*nBasis+j].q   = q;
		mpole[i*nBasis+j].x   = x;
		mpole[i*nBasis+j].y   = y;
		mpole[i*nBasis+j].z   = z;
	}

	return mpole;
}


//
// coulombMultipole : compute Coulomb interaction energy between two multipole.
// It returns zero if failed to do so because q = 0.
//
// 24 Dec, 2012 - Teepanis Chachiyo
//    Initital implementation and testing
//
double coulombMultipole(
	int nBasis,                       // number of basis function
	int p, int q,                     // index to first multipole
	int i, int j,                     // index to second multipole
	const struct Multipole_t *mpole){ // pointer to multipole info

	double K;                         // inverse of sepration distance
	const struct Multipole_t *m1,*m2; // pointer to multipole structure

	m1 = &(mpole[p*nBasis+q]);
	m2 = &(mpole[i*nBasis+j]);

	// check if charge is zero
	if(m1->q == 0.0 || m2->q == 0.0) return 0.0;

	// check if separation is large enough
	K   = (m1->x - m2->x)*(m1->x - m2->x) + 
	      (m1->y - m2->y)*(m1->y - m2->y) +
	      (m1->z - m2->z)*(m1->z - m2->z);
	if(K < MULTIPOLE_RADII2_CUTOFF) return 0.0;
	K  = sqrt(1/K);

	// monopole-monopole
	return (m1->q * m2->q * K);
}
