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

#ifndef LIN_H
#define LIN_H
int gen_sym_eigen(
	int nDim,                           // array dimension
	const double *A, const double *S,   // input A and S
	double *e, double *C);              // output value and vector

int gen_sym_SomeEigen(
	int nDim,                           // array dimension
	int nEigen,                         // number of eigen value requested
	const double *A, const double *S,   // input A and S
	double *e, double *C);              // output value and vector

int linear_solver(int nRow,         // number of row
                  int nCol,         // number column
                  const double *A,  // input matrix A
                  double *B);       // input matrix B, output matrix X
#endif
