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

#include <math.h>
#include "fgamma.h"

// fgamma_s : compute Fm(x) using series expansion
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
//		Initial implementation and testing
//
#define CONV 1.0E-15
#define ITMAX 100
double fgamma_s(int m, double x){
	double a;
	double sum=0.0;
	int i;

	m = m+m+1;
	a = 1.0/m;
	x = 2.0*x;
	for(i=1; i < ITMAX; i++){
		sum += a;
		a    = a*x/(m+i+i);
		if(a<CONV) break;
	}
	return exp(-x/2.0)*sum;
}
#undef CONV
#undef ITMAX

// fgamma_steed : compute Fm(x) using continued fraction
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation and testing
//
#define ITMAX 100
#define CONV  1.0e-15
double fgamma_steed(int m, double x){
	int i;
	double a,bi,ai,D,df,f;
	  
	// compute continued fraction
	// using Steed's method
	a = m+0.5;
	bi = x+1.0-a;
	ai = 0.0;
	D=f=df=1.0/bi;
	for(i=1; i <= ITMAX; i++){
		ai += a-i-i+1;
		bi += 2.0;
		D   = 1.0/(bi+ai*D);
		df *= (bi*D-1.0);
		f  += df;
		if(fabs(df/f) < CONV) break;
	}
	// compute 1*3*5...(2m-1)/(2x)^m
	D  = 1.0;
	a  = 1.0/(x+x);
	for(i=1; i <= (m+m-1); i+=2)
		D *= i*a;

#define SQRT_PI_OVER_2 (double)0.8862269254528
	D *= SQRT_PI_OVER_2/sqrt(x);
#undef  SQRT_PI_OVER_2
	return D - 0.5*exp(-x)*f;
}
#undef ITMAX
#undef CONV

// fgamma_lenz : compute Fm(x) using continued fraction
// the continued fraction algorithm was adapted from 
// "Numerical Recipies in C++."
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation & testing
//
#define ITMAX 100
#define CONV  1.0e-15
#define FPMIN 1.0E-30
double fgamma_lenz(int m, double x){
	int i;
	double a,an,b,c,d,del,h;
	  
	// compute continued fraction
	// using Lenz's method
	a = m+0.5;
	b = x+1.0-a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	for(i=1; i <= ITMAX; i++){
		an  = (a-i)*i;
		b  += 2;
		d   = an*d+b;
		if(fabs(d)<FPMIN) d = FPMIN;
		c   = b+an/c;
		if(fabs(c)<FPMIN) c = FPMIN;
		d   = 1.0/d;
		del = d*c;
		h  *= del;
		if(fabs(del-1.0) < CONV) break;
	}
	// compute 1*3*5...(2m-1)/(2x)^m
	d  = 1.0;
	c  = 1.0/(x+x);
	for(i=1; i <= (m+m-1); i+=2)
		d *= i*c;

#define SQRT_PI_OVER_2 (double)0.8862269254528
	d *= SQRT_PI_OVER_2/sqrt(x);
#undef  SQRT_PI_OVER_2
	return d - 0.5*exp(-x)*h;
}
#undef ITMAX
#undef CONV
#undef FPMIN

// fgamma_0 :  compute F0(x) using erf and sqrt relation.
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation & testing
//
// May 19, 2011 - Teepanis Chachiyo
//      handle case x approaches zero
double fgamma_0(double x){

	double t;

	// case x approaches zero
	if(x < 0.0005)
		return    1.0 - x/3.0 + x*x/10.0;

	t = sqrt(x);
	return (double)0.8862269254528*erf(t)/t;
}


// fgamma :  compute Fm(x) using power serie or continued 
// fraction depending on the range of x.
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation & testing
//
double fgamma(int m, double x){
	if(x < (m+0.5))
		return fgamma_s(m, x);
	else
		return fgamma_steed(m, x);
}


// fgamma_set : compute a range of Fm(x) from m=0 upto
// (including) m=max using recurrence relation.
//
// 2008 - Teepanis Chachiyo & Narong Putuddong
// 		Initial implementation  & testing
//
void fgamma_set(int max, double x, double *fset){
	double a,b,c;
	int i,j;

	// case x approaches zero
	if(x < 0.0005){
		for(i=0; i <= max; i++)
		fset[i] =     1.0/(i+i    +1.0)
		            -   x/(i+i    +3.0)
		            + x*x/(i+i+i+i+10.0);
		return;
	}

	// compute F0(x) first
	fset[0] = fgamma_0(x);

	if(max < 1) return;

	// forward
	a = 1.0/(x+x);
	c = a*exp(-x);
	for(i=1; i <= max; i++){
		if((b=(i+i-1)*a)>1.0) break;
		fset[i] = b*fset[i-1] - c;
	}

	if(max < i) return;

	// backward
	fset[max] = fgamma(max,x);
	a = 1.0/a;
	c = c*a;
	for(j=max-1; j >= i; j--)
		fset[j] = (a*fset[j+1] + c)/(j+j+1);
}

