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

#ifndef UTIL_H
#define UTIL_H

#define BOHR2ANGSTROM 0.529177249000
#define ANGSTROM2BOHR 1.889725988579

#define HARTREE2EV    27.2113957

#define DEBYE2AU 0.393430307
#define AU2DEBYE 2.54174623

int time_str(int max, char *str);
int findf(FILE *fd, int n, char *str, ...);
#define SYMB_SHORTNAME 0
#define SYMB_LONGNAME  1
int sym2Z(const char *sym, int type);
void Z2Sym(int Z, char *sym);
void Z2SymShort(int Z, char *sym);
double pow_int(double x,int i);
void stepProgressBar(int maxbar, int maxn);
#endif
