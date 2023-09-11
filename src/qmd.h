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
#ifndef QMD_H
#define QMD_H

#include "basis.h"
#include "mol.h"
#include "option.h"

/////////////////////////////////////////////////////////
////////INTERNAL VARIABLES MUST BE IN MOLECULAR DYNMAICS UNIT
//
// DISANCE  NANO-METER
// TIME     PICO-SECOND
// VELOCITY NM/PS
// MASS     DALTON
// FORCE    KJ/MOL/NANO-METER
// TEMPERATURE KELVIN
// BOLTZMAN 0.0083144621 KJ/MOL/TEMP
/////////////////////////////////////////////////////////

#define HARTREE2KJMOL 2625.49962      // Hartere to kJ/mol conversion
#define MD_kB            0.0083144621 // Boltzmann constant in MD unit
#define MD_LIGHTSPEED 29992.458       // speed of light in MD unit

#define AU2MD_FORCE    (HARTREE2KJMOL*ANGSTROM2BOHR/0.1)
#define AU2MD_ENERGY   (HARTREE2KJMOL)
#define AU2MD_DISTANCE (0.1/ANGSTROM2BOHR)

#define AU2VOLTPERNM   514.220652  // electric field AU --> V/nm

double Z2mass(int Z);

void qmd(
	int dbSize,                         // number of record in basis set db
	const struct GTOBasisSet_t *basisDB,// pointer to basis set database
	struct Molecule_t *mol,             // returned molecular coordinate
	struct option_t *opt);              // options

#endif
