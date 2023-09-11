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

#ifndef FGAMMA_H
#define FGAMMA_H
double fgamma_s(int m, double x);
double fgamma_steed(int m, double x);
double fgamma_lenz(int m, double x);
double fgamma_0(double x);
double fgamma(int m, double x);
void fgamma_set(int max, double x, double *fset);
#endif
