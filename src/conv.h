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

#ifndef CONV_H
#define CONV_H
double conv_damping(int nBasis,     // number of basis
                    double alpha,   // drag coefficient between 0 and 1
                    double *PA,     // alpha density matrix
                    double *PB);    // beta density matrix

double conv_diis2(int nBasis,      // number of basis
                  double alpha,    // drag coefficient between 0 and 1
                  double *PA,      // alpha density matrix
                  double *PB);     // beta density matrix

double conv_diis3(int nBasis,      // number of basis
                  double alpha,    // drag coefficient between 0 and 1
                  double *PA,      // alpha density matrix
                  double *PB);     // beta density matrix

double conv_diis4(int nBasis,      // number of basis
                  double alpha,    // drag coefficient between 0 and 1
                  double *PA,      // output alpha density matrix
                  double *PB);     // output beta  density matrix

#endif
