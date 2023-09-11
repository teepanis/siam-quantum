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
#include <stdlib.h>
#include <string.h>
#include "int.h"
#include "fgamma.h"


///////////////////////////////////////////////////////////////////
// Credits: The codes below started with the implementation      //
// of PyQuante by Richard Muller. Since then, many improvments   //
// has been made in order to make it faster. The current version //
// below is based on THO. paper. Nevertheless, I shall give      //
// Richard Muller tremendous credits here for his inspiration.   //
///////////////////////////////////////////////////////////////////

// f_coef : computes array of coefficient x^j in the expansion 
// (x+a)^l (x+b)^m. The usual binomial expression has been 
// recasted into a recurrence relation instead.
//
// 2008 - Teepanis Chachiyo 
// 		Initial implementation & testing
//
#define MAX_POWER 32
static void f_coef(double *f, int l, int m, double a, double b){
	double Ca[MAX_POWER], Cb[MAX_POWER];
	int i,j,k;

	Ca[l] = 1.0;
	for(j=l-1;j>=0;j--)
		Ca[j] = a*(j+1)/(l-j)*Ca[j+1];
	Cb[m] = 1.0;
	for(k=m-1;k>=0;k--)
		Cb[k] = b*(k+1)/(m-k)*Cb[k+1];
	for(i=0;i<=(l+m);i++)
		f[i] = 0.0;
	for(j=0;j<l+1;j++)
		for(k=0;k<m+1;k++){
			i     = j+k;
			f[i] += Ca[j]*Cb[k];
		}
	return;
}
#undef MAXPOWER

// overlap_1d : compute overlap integral in one dimension.
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
#define MAXL 8
static double overlap_1d(int l1, int l2, 
                         double PAx, double PBx,
                         double eta){
	int i;
	double f[2*MAXL];
	double sum=0.0;
	double p;

	// compute coefficient
	f_coef(f, l1, l2, PAx, PBx);

	// loop and using look-up tables
	eta = eta/2.0;
	p   = 1.0;
	for(i=0; i < (l1+l2)/2+1; i++){
		p        *= fabs(i+i-1);
		sum      += f[i+i]*p;
		p        *= eta;
	}
	return sum;
}
#undef MAXL

// overlap :  compute overlap matrix using equation 2.12 from 
// THO paper. The integral is of the form (a|b), where "a" and "b" 
// are GTO functions. 
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+\
                                  (y1-y2)*(y1-y2)+\
                                  (z1-z2)*(z1-z2))
double overlap(double alpha1, int l1, int m1, int n1,
               double xa, double ya, double za,
               double alpha2, int l2, int m2, int n2,
               double xb, double yb, double zb){

	double rab2,eta,xp,yp,zp,pre,wx,wy,wz;

	rab2 = DIST2(xa,ya,za,xb,yb,zb);
	eta  = 1.0/(alpha1+alpha2);
	xp   = (alpha1*xa+alpha2*xb)*eta;
	yp   = (alpha1*ya+alpha2*yb)*eta;
	zp   = (alpha1*za+alpha2*zb)*eta;

	wx = overlap_1d(l1,l2,xp-xa,xp-xb,eta);
	wy = overlap_1d(m1,m2,yp-ya,yp-yb,eta);
	wz = overlap_1d(n1,n2,zp-za,zp-zb,eta);

	// prefactor
	pre = M_PI*eta*sqrt(M_PI*eta)*exp(-alpha1*alpha2*rab2*eta);

	return pre*wx*wy*wz;
}


// f3_coef : computes array of coefficient x^j in the expansion 
// (x+a)^l (x+b)^m (x+c)^n. The usual binomial expression has been 
// recasted into a recurrence relation instead.
//
// 2011 - Aniwat Kesorn and Teepanis Chachiyo
//      Adapted from f_coef
//
// Oct 13, 2012 - Teepanis Chachiyo and Aniwat Kesorn
//      Added to Siam Quantum
//
#define MAX_POWER 32
static void f3_coef(double *f, int l, int m, int n, double a, double b, double c){
	double Ca[MAX_POWER], Cb[MAX_POWER], Cc[MAX_POWER];
	int iA,iB,iC;
	int I;

	// coefficient for (x+a)^l
	Ca[l] = 1.0;
	for(iA=l-1;iA>=0;iA--)
		Ca[iA] = a*(iA+1)/(l-iA)*Ca[iA+1];

	// coefficient for (x+b)^m
	Cb[m] = 1.0;
	for(iB=m-1;iB>=0;iB--)
		Cb[iB] = b*(iB+1)/(m-iB)*Cb[iB+1];

	// coefficient for (x+c)^n
	Cc[n] = 1.0;
	for(iC=n-1;iC>=0;iC--)
		Cc[iC] = c*(iC+1)/(n-iC)*Cc[iC+1];

	// reset f[I]
	for(I=0;I<=(l+m+n);I++) f[I] = 0.0;

	// multiply (x+a)^l (x+b)^m (x+c)^n together
	for(iA=0;iA<l+1;iA++)
	for(iB=0;iB<m+1;iB++)
	for(iC=0;iC<n+1;iC++){
		I     = iA+iB+iC;
		f[I] += Ca[iA]*Cb[iB]*Cc[iC];
	}
	return;
}
#undef MAXPOWER

// moment_1d : compute moment integral in one dimesion
//
// 2011 - Aniwat Kesorn and Teepanis Chachiyo
//     Adapted from overlap_1d subroutine
//
// Oct 13, 2012 - Teepanis Chachiyo and Aniwat Kesorn
//     Added to Siam Quantum
//
#define MAXL 8
#define MAXMOMENT 8
static double moment_1d(int l1, int l2, int mx,
                        double PAx, double PBx, double PCx,
                        double eta){
	int i;
	double f[2*MAXL+MAXMOMENT];
	double sum=0.0;
	double p;

	// compute coefficient
	f3_coef(f, l1, l2, mx, PAx, PBx, PCx);

	// loop and using look-up tables
	eta = eta/2.0;
	p   = 1.0;
	for(i=0; i < (l1+l2+mx)/2+1; i++){
		p        *= fabs(i+i-1);
		sum      += f[i+i]*p;
		p        *= eta;
	}
	return sum;
}
#undef MAXMOMENT
#undef MAXL

// moment :  compute moment matrix element adapted from equation 2.12 from 
// THO paper. This subroutine evaluate the integral
//
// / 3              mx      my      mz
// |d r gto(1) (x-xc)  (y-yc)  (z-zc)  gto(2)
// /
//
// where gto(1) and gto(2) are simple gaussian function with the exponent
// alpha1 and alpha2 respectively.
//
// 2011 - Aniwat Kesorn
//     Adapted from overlap subroutine
//
// Oct 13, 2012 - Teepanis Chachiyo and Aniwat Kesorn
//     Add to Siam Quantum
//
double moment(
	double alpha1, int l1, int m1, int n1, // gto(1) exponent and angular index
	double xa, double ya, double za,       // gto(1) center
	double alpha2, int l2, int m2, int n2, // gto(2) exponent and angular index
	double xb, double yb, double zb,       // gto(2) center
	int mx, int my, int mz,                // moment angular index
	double xc, double yc, double zc){      // moment center

	double rab2,eta,xp,yp,zp,pre,wx,wy,wz;

	rab2 = DIST2(xa,ya,za,xb,yb,zb);
	eta  = 1.0/(alpha1+alpha2);
	xp   = (alpha1*xa+alpha2*xb)*eta;
	yp   = (alpha1*ya+alpha2*yb)*eta;
	zp   = (alpha1*za+alpha2*zb)*eta;

	wx = moment_1d(l1,l2,mx,xp-xa,xp-xb,xp-xc,eta);
	wy = moment_1d(m1,m2,my,yp-ya,yp-yb,yp-yc,eta);
	wz = moment_1d(n1,n2,mz,zp-za,zp-zb,zp-zc,eta);

	// prefactor
	pre = M_PI*eta*sqrt(M_PI*eta)*exp(-alpha1*alpha2*rab2*eta);

	return pre*wx*wy*wz;
}


// kinetic : compute kinetic operator integral of the form
// (a|-1/2 LAPLACIAN|b), where "a" and "b" are
// GTO. Again, the actual reference is the equation
// 2.13 where the kinetic integral is expanded in 
// terms of several overlap integrals. The programming
// structure and variables are from PyQuante. No
// optimization here because it has been done in 
// the overlap integral already.
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
// March 15, 2010 - Teepanis Chachiyo
//      Bug fix in case l2<2 or m2<2 or n2<2 because
//      in such case the overlap will segfault
//
double kinetic(double alpha1, int l1, int m1, int n1,
               double xa, double ya, double za,
               double alpha2, int l2, int m2, int n2,
               double xb, double yb, double zb){

	double term0,term1,term2;

	term0 = alpha2*(2*(l2+m2+n2)+3)
	              *overlap(alpha1,l1,m1,n1,xa,ya,za,
	                       alpha2,l2,m2,n2,xb,yb,zb);

	term1 = -2*alpha2*alpha2
	          *(  overlap(alpha1,l1,m1,n1,xa,ya,za,
	                      alpha2,l2+2,m2,n2,xb,yb,zb)
	            + overlap(alpha1,l1,m1,n1,xa,ya,za,
	                      alpha2,l2,m2+2,n2,xb,yb,zb)
	            + overlap(alpha1,l1,m1,n1,xa,ya,za,
	                      alpha2,l2,m2,n2+2,xb,yb,zb)
	           );

	term2 = -0.5*(  ((l2<2)?0:l2*(l2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
	                                           alpha2,l2-2,m2,n2,xb,yb,zb))
	              + ((m2<2)?0:m2*(m2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
	                                           alpha2,l2,m2-2,n2,xb,yb,zb))
	              + ((n2<2)?0:n2*(n2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
	                                           alpha2,l2,m2,n2-2,xb,yb,zb))
	             );

  return term0+term1+term2;
}


// naiA_1d_tho : compute auxilary constant according to 
// equation 2.18 in reference THO. This is known as constant
// A in the paper. The explicit expression is given in the
// paper.
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
#define MAXL 8
static double FACT[16] = {
  1.000000000000E+00 ,
  1.000000000000E+00 ,
  2.000000000000E+00 ,
  6.000000000000E+00 ,
  2.400000000000E+01 ,
  1.200000000000E+02 ,
  7.200000000000E+02 ,
  5.040000000000E+03 ,
  4.032000000000E+04 ,
  3.628800000000E+05 ,
  3.628800000000E+06 ,
  3.991680000000E+07 ,
  4.790016000000E+08 ,
  6.227020800000E+09 ,
  8.717829120000E+10 ,
  1.307674368000E+12 };

void naiA_1d(double *A,  int l1,     int l2,
             double PAx, double PBx, double CPx,
             double eta){
	int i,r,u,I;            // loop counter
	double f[2*MAXL];       // expansion coefficient
	double CPx_pow[2*MAXL]; // (eta/4)^I
	double eta_pow[2*MAXL]; // CPx^I / I!

	// compute coefficient
	f_coef(f, l1, l2, PAx, PBx);

	// set initial values
	eta        = eta/4.0;	
	A[0]       = 0.0;
	eta_pow[0] = 1.0;
	CPx_pow[0] = 1.0;

	for(I=1; I < l1+l2+1; I++){
		A[I] = 0.0;
		eta_pow[I] = eta_pow[I-1]*eta;
		CPx_pow[I] = CPx_pow[I-1]*CPx/I;
	}

	for(i=0; i < (l1+l2)+1;   i++)
	for(r=0; r < i/2+1;       r++)
	for(u=0; u < (i-r-r)/2+1; u++){
		I     = i-r-r-u;
		A[I] += f[i]
		       *((i+u)&1?-1:1)
		       *FACT[i]/(FACT[r]*FACT[u])
		       *CPx_pow[I-u]
		       *eta_pow[r+u];
	}

	return;
}
#undef MAXL

// nai_tho : compute nuclear attraction integral of the form
// -(a|1/rC|b). Notice the negative sign. This method is based
// on THO reference in equation 2.17-2.18.
//
// Dec 16, 2015 - Teepanis Chachiyo
// 		Add pre-screening (<PRIMITIVE_CUTOFF)
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
#define MAXL 8
double nai(double x1, double y1, double z1, double norm1,
           int l1, int m1, int n1, double alpha1,
           double x2, double y2, double z2, double norm2,
           int l2, int m2, int n2, double alpha2,
           double x3, double y3, double z3){

	double eta,xp,yp,zp,sum,rab2,rcp2;
	double Ax[2*MAXL],Ay[2*MAXL],Az[2*MAXL];
	double F[2*MAXL];
	double pre;
	int kX, kY, kZ;

	eta = 1.0/(alpha1+alpha2);

	xp   = (alpha1*x1+alpha2*x2)*eta;
	yp   = (alpha1*y1+alpha2*y2)*eta;
	zp   = (alpha1*z1+alpha2*z2)*eta;

	rab2 = DIST2(x1,y1,z1,x2,y2,z2);
	rcp2 = DIST2(x3,y3,z3,xp,yp,zp);

	// cut-off
	pre = norm1*norm2*2*M_PI*eta*exp(-alpha1*alpha2*rab2*eta);
	if(pre < PRIMITIVE_CUTOFF) return 0.0;

	// integrate all 1d directions
	naiA_1d(Ax,l1,l2,xp-x1,xp-x2,xp-x3,eta);
	naiA_1d(Ay,m1,m2,yp-y1,yp-y2,yp-y3,eta);
	naiA_1d(Az,n1,n2,zp-z1,zp-z2,zp-z3,eta);

	// compute all fgamma
	fgamma_set(l1+l2+m1+m2+n1+n2,rcp2/eta,F);

	// looping
	sum = 0.0;
	for(kX=0; kX < l1+l2+1; kX++)
	for(kY=0; kY < m1+m2+1; kY++)
	for(kZ=0; kZ < n1+n2+1; kZ++)
		sum += Ax[kX]*Ay[kY]*Az[kZ]*F[kX+kY+kZ];

	return -pre*sum;
}
#undef MAXL

#define MAXL 8
static void efi_1d(double *A,  double *Ap, int l1,     int l2,
                   double PAx, double PBx, double CPx,
                   double eta){
	register int i,r,u,I;   // loop counter
	double f[2*MAXL];       // expansion coefficient
	double CPx_pow[2*MAXL]; // (eta/4)^I
	double eta_pow[2*MAXL]; // CPx^I / I!

	// compute coefficient
	f_coef(f, l1, l2, PAx, PBx);

	// set initial values
	eta        = eta/4.0;	
	A[0]       = 0.0;
	Ap[0]      = 0.0;
	eta_pow[0] = 1.0;
	CPx_pow[0] = 1.0;

	for(I=1; I < l1+l2+1; I++){
		A[I]  = 0.0;
		Ap[I] = 0.0;
		eta_pow[I] = eta_pow[I-1]*eta;
		CPx_pow[I] = CPx_pow[I-1]*CPx;
	}

	for(i=0; i < (l1+l2)+1;   i++)
	for(r=0; r < i/2+1;       r++)
	for(u=0; u < (i-r-r)/2+1; u++){
		I     = i-r-r-u;
		A[I] += f[i]
		       *((i+u)&1?-1:1)
		       *FACT[i]/(FACT[r]*FACT[u]*FACT[I-u])
		       *CPx_pow[I-u]
		       *eta_pow[r+u];

		Ap[I] += f[i]
		       *((i+u+1)&1?-1:1)
		       *FACT[i]/(FACT[r]*FACT[u]*FACT[I-u])
		       *(I-u)
		       *CPx_pow[I-u-1]
		       *eta_pow[r+u];
	}

	return;
}
#undef MAXL

// compute electric field integral of the form
// (a|ddxC(1/rC)|b), (a|ddyC(1/rC)|b), (a|ddzC(1/rC)|b)
//
// 2010 - Theerapon Khumla
//   Initial implementation
// 
#define MAXL 8
void efi(double x1, double y1, double z1, double norm1,
         int l1, int m1, int n1, double alpha1,
         double x2, double y2, double z2, double norm2,
         int l2, int m2, int n2, double alpha2,
         double x3, double y3, double z3,
         double *ex, double *ey, double *ez){

	double eta,xp,yp,zp,sum,rab2,rcp2;
	double  Ax[2*MAXL], Ay[2*MAXL], Az[2*MAXL];
	double Axp[2*MAXL],Ayp[2*MAXL],Azp[2*MAXL];
	double F[2*MAXL];
	int kX, kY, kZ;

	eta = 1.0/(alpha1+alpha2);

	xp   = (alpha1*x1+alpha2*x2)*eta;
	yp   = (alpha1*y1+alpha2*y2)*eta;
	zp   = (alpha1*z1+alpha2*z2)*eta;

	rab2 = DIST2(x1,y1,z1,x2,y2,z2);
	rcp2 = DIST2(x3,y3,z3,xp,yp,zp);

	// integrate all 1d directions
	efi_1d(Ax,Axp,l1,l2,xp-x1,xp-x2,xp-x3,eta);
	efi_1d(Ay,Ayp,m1,m2,yp-y1,yp-y2,yp-y3,eta);
	efi_1d(Az,Azp,n1,n2,zp-z1,zp-z2,zp-z3,eta);

	// compute all fgamma
	fgamma_set(l1+l2+m1+m2+n1+n2+1,rcp2/eta,F);

	*ex = 0.0; *ey = 0.0; *ez = 0.0;
	for(kX=0; kX < l1+l2+1; kX++)
	for(kY=0; kY < m1+m2+1; kY++)
	for(kZ=0; kZ < n1+n2+1; kZ++){
		*ex += Axp[kX]*Ay[kY]*Az[kZ]*F[kX+kY+kZ];
		*ey += Ax[kX]*Ayp[kY]*Az[kZ]*F[kX+kY+kZ];
		*ez += Ax[kX]*Ay[kY]*Azp[kZ]*F[kX+kY+kZ];
	}

	// looping
	sum = 0.0;
	for(kX=0; kX < l1+l2+1; kX++)
	for(kY=0; kY < m1+m2+1; kY++)
	for(kZ=0; kZ < n1+n2+1; kZ++)
		sum += Ax[kX]*Ay[kY]*Az[kZ]*F[kX+kY+kZ+1];

	// add contribution
	*ex -= sum*2.0*(x3-xp)/eta;
	*ey -= sum*2.0*(y3-yp)/eta;
	*ez -= sum*2.0*(z3-zp)/eta;

	// multiply prefactor
	sum = -norm1*norm2*2.0*M_PI*eta*exp(-alpha1*alpha2*rab2*eta);
	*ex = *ex * sum;
	*ey = *ey * sum;
	*ez = *ez * sum;
}
#undef MAXL

#define MAXL 8
static inline void get_H(double *H, const double *fi, int l1, int l2, double eta){

	int L,i,r;               // loop counter
	double eta_pow[2*MAXL];  // (eta/4)^I

/* /////////////////////////////////
	// special case: s type
	if(l1+l2==0){
		H[0] = fi[0];
		return;
	}
    
	// rescale
	eta = eta/4.0;
	
	// special case: p type
	if(l1+l2==1){
		H[0] = fi[0];
		H[1] = fi[1]*eta;
		return;
	}
	
	// special case: d type
	if(l1+l2==2){
		H[0]  = fi[0];
		H[1]  = fi[1]*eta;
		H[2]  = fi[2]*eta;
		H[0] += H[2]+H[2];
		H[2] *= eta;
		return;
	}
//////////////////////////////// */

	switch(l1+l2){

	case 0:
		H[0] = 1.0;
		return; 
		break;

	case 1:
		eta = eta/4.0;
		H[0] = fi[0];
		H[1] = fi[1]*eta;
		return;
	break;

	case 2:
		eta = eta/4.0;
		H[0]  = fi[0];
		H[1]  = fi[1]*eta;
		H[2]  = fi[2]*eta;
		H[0] += H[2]+H[2];
		H[2] *= eta;
		return;
	break;

	case 3:
		eta   = eta/4.0;
		H[0]  = fi[0];
		H[1]  = fi[1]*eta;
		H[2]  = fi[2]*eta;
		H[0] += H[2]+H[2];
		H[2] *= eta;
		H[3]  = fi[3]*eta*eta;
		H[1] += 6.0*H[3];
		H[3] *= eta; 
	break;

	case 4:
		eta   = eta/4.0;
		H[0]  = fi[0];
		H[1]  = fi[1]*eta;
		H[2]  = fi[2]*eta;
		H[0] += H[2]+H[2];
		H[2] *= eta;
		H[3]  = fi[3]*eta*eta;
		H[1] += 6.0*H[3];
		H[3] *= eta;
		H[4]  = fi[4]*eta*eta;
		H[0] += 12.0*H[4];
		H[2] += 12.0*H[4]*eta;
		H[4] *= eta*eta;
	break;

	default:
		eta = eta/4.0;
	break;
	}

	H[0]           = 0.0;
	eta_pow[0]     = 1.0;
	for(L=1; L < l1+l2+1; L++){
		H[L]       = 0.0;
		eta_pow[L] = eta_pow[L-1]*eta;
	}

	for(i=0; i < l1+l2+1; i++)
	for(r=0; r < i/2+1;   r++){
		L = i-r-r;
		H[L] = H[L] +  fi[i]
		                *FACT[i]/(FACT[r]*FACT[L])
		                *eta_pow[i-r];
	}
	return;
}

// eri : compute two electron integral of the form (ab|1/r12|cd).
//
// May 18, 2011 - Teepanis Chachiyo
//                Clean the special cases
//
// July 4, 2010 - Teepanis Chachiyo
//                Add special case to enhance speed
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
void eriB_1d(double *B, int l1, int l2, int l3, int l4,
             double Px, double Ax, double Bx,
             double Qx, double Cx, double Dx,
             double eta1, double eta2,
             double zeta){

	int L,M,u;               // looping index
	int I;                   // B array indices
	double fi1[2*MAXL];      // binomial array
	double fi2[2*MAXL];      // binomial array
	double PA,PB,QC,QD,QP;   // auxilary
	double HL[2*MAXL];       // l1 l2 part
	double HM[2*MAXL];       // l3 l4 part
	double QP_pow[4*MAXL];   // QP^I / I!
	double zeta_pow[4*MAXL]; // zeta^I

	// special cases
#define ERI1111(Ax,Bx,Cx,Dx,Px,Qx,eta1,eta2,zeta)                      \
		PA = Px-Ax;                                                    \
		PB = Px-Bx;                                                    \
		QC = Qx-Cx;                                                    \
		QD = Qx-Dx;                                                    \
		QP = Qx-Px;                                                    \
		eta1 *= 0.5;                                                   \
		eta2 *= 0.5;                                                   \
		zeta *= 0.5;                                                   \
		B[0]  = (eta1+PA*PB)*(eta2+QC*QD);                             \
		B[1]  = eta1*(PA+PB)*(QP*(eta2+QC*QD) + eta2*(QC+QD));         \
		B[1] -= eta1*eta1*(eta2+QC*QD);                                \
		B[1] -= eta2*(eta1+PA*PB)*(eta2 + QP*(QC+QD));                 \
		B[1] *= zeta;                                                  \
		B[2]  = eta2*eta2*(eta1+PA*PB) + eta1*eta1*(eta2+QC*QD);       \
		B[2] -= eta1*eta2*(PA+PB)*(QC+QD);                             \
		B[2] *= QP*QP;                                                 \
		B[2] += 3.0*QP*eta1*eta2*(eta1*(QC+QD) - eta2*(PA+PB));        \
		B[2] += 3.0*eta1*eta1*eta2*eta2;                               \
		B[2] *= zeta*zeta;                                             \
		B[3]  = QP*eta2*(PA+PB) - QP*eta1*(QC+QD) - 6.0*eta1*eta2;     \
		B[3] *= QP*QP*zeta*zeta*zeta*eta1*eta2;                        \
		B[4]  = QP*QP*eta1*eta2*zeta*zeta;                             \
		B[4]  = B[4]*B[4];

#define ERI2010(Ax,Cx,Px,Qx,eta1,eta2,zeta)                            \
		PA = Px-Ax;                                                    \
		QC = Qx-Cx;                                                    \
		QP = Qx-Px;                                                    \
		eta1 *= 0.25;                                                  \
		eta2 *= 0.25;                                                  \
		B[0]  = QC*(PA*PA+eta1+eta1);                                  \
		B[1]  = 2.0*eta1*eta2*(PA+PA-QP) + 2.0*eta1*QC*(PA*QP-eta1);   \
		B[1] -= PA*PA*QP*eta2;                                         \
		B[1] *= zeta;                                                  \
		B[2]  = eta1*eta2*6.0 - 2.0*PA*QP*eta2 + QC*QP*eta1;           \
		B[2] *= QP*eta1*zeta*zeta;                                     \
		QP   *= QP*QP;                                                 \
		zeta *= zeta*zeta;                                             \
		B[3]  = -QP*zeta*eta1*eta1*eta2;

#define ERI2100(Ax,Bx,Px,Qx,eta1,zeta)                                 \
		PA = Px-Ax;                                                    \
		PB = Px-Bx;                                                    \
		QP = Qx-Px;                                                    \
		eta1 *= 0.25;                                                  \
		B[0]  = 2.0*eta1*(PA+PA+PB) + PA*PA*PB;                        \
		B[1]  = eta1*(6.0*QP-4.0*PA-PB-PB) + PA*QP*(PA+PB+PB);         \
		B[1] *= zeta*eta1;                                             \
		B[2]  = QP*(PA+PA+PB) - 6.0*eta1;                              \
		B[2] *= QP*zeta*zeta*eta1*eta1;                                \
		B[3]  = QP*eta1*zeta;                                          \
		B[3] *= B[3]*B[3];

#define ERI0111(Bx,Cx,Dx,Px,Qx,eta1,eta2,zeta)                   \
		PB = Px-Bx;                                              \
		QC = Qx-Cx;                                              \
		QD = Qx-Dx;                                              \
		QP = Qx-Px;                                              \
		zeta  *= 0.25;                                           \
		eta1  *= 0.50;                                           \
		eta2  *= 0.50;                                           \
		B[0]  = PB*(eta2 + QC*QD);                               \
		B[1]  = eta1*eta2*(QC+QD+QP);                            \
		B[1] += eta1*QC*QD*QP;                                   \
		B[1] -= eta2*PB*QP*(QC+QD);                              \
		B[1] -= eta2*eta2*PB;                                    \
		B[1] *= 2.0*zeta;                                        \
		B[2]  = eta1*QP*(QC+QD) + 3.0*eta1*eta2 - PB*QP*eta2;    \
		B[2] *= -4.0*QP*eta2*zeta*zeta;                          \
		QP   *= QP*QP;                                           \
		zeta *= zeta*zeta;                                       \
		B[3]  = 8.0*QP*eta1*eta2*eta2*zeta;

#define ERI2000(Ax,Px,Qx,eta1,zeta)             \
		PA = Px-Ax;                             \
		QP = Qx-Px;                             \
		B[0] =  0.5*eta1 + PA*PA;               \
		B[1] = -0.125*zeta*eta1*(eta1-4*QP*PA); \
		B[2] =  0.25*QP*zeta*eta1;              \
		B[2] =  B[2]*B[2];

#define ERI1100(Ax,Bx,Px,Qx,eta1,zeta)                \
		PA = Px-Ax;                                   \
		PB = Px-Bx;                                   \
		QP = Qx-Px;                                   \
		B[0] =  0.5*eta1+PA*PB;                       \
		B[1] = -0.125*zeta*eta1*(eta1-2*QP*(PA+PB));  \
		B[2] =  0.25*QP*eta1*zeta;                    \
		B[2] =  B[2]*B[2];

#define ERI1010(Ax,Cx,Px,Qx,eta1,eta2,zeta)                    \
		PA = Px-Ax;                                            \
		QC = Qx-Cx;                                            \
		QP = Qx-Px;                                            \
		B[0] = PA*QC;                                          \
		B[1] = 0.125*zeta*(eta1*eta2+2*QP*(QC*eta1-PA*eta2));  \
		B[2] = 0.25*QP*zeta;                                   \
		B[2] = -eta1*eta2*B[2]*B[2];

#define ERI1000(Ax,Px,Qx,eta1,zeta)     \
		B[0] = Px-Ax;                   \
		B[1] = 0.25*(Qx-Px)*eta1*zeta;

	/*
	I=l1+l2+l3+l4;
	switch(I){
	case 2: if(l1==2){         ERI2000(Ax,Px,Qx,eta1,zeta);    return; //(20|00)
			}else if(l2==2){   ERI2000(Bx,Px,Qx,eta1,zeta);    return; //(02|00)
			}else if(l3==2){   ERI2000(Cx,Qx,Px,eta2,zeta);    return; //(00|20)
			}else if(l4==2){   ERI2000(Dx,Qx,Px,eta2,zeta);    return; //(00|02)
			}else if(!(l3+l4)){ERI1100(Ax,Bx,Px,Qx,eta1,zeta); return; //(11|00)
			}else if(!(l1+l2)){ERI1100(Cx,Dx,Qx,Px,eta2,zeta); return; //(00|11)
			}else if(!(l2+l4)){ERI1010(Ax,Cx,Px,Qx,eta1,eta2,zeta); return; //(10|10)
			}else if(!(l1+l4)){ERI1010(Bx,Cx,Px,Qx,eta1,eta2,zeta); return; //(01|10)
			}else if(!(l2+l3)){ERI1010(Ax,Dx,Px,Qx,eta1,eta2,zeta); return; //(10|01)
			}else if(!(l1+l3)){ERI1010(Bx,Dx,Px,Qx,eta1,eta2,zeta); return; //(01|01)
			}
	break;

	case 1: if(l1==1){          ERI1000(Ax,Px,Qx,eta1,zeta); return; //(10|00)
			}else if(l2==1){    ERI1000(Bx,Px,Qx,eta1,zeta); return; //(01|00)
			}else if(l3==1){    ERI1000(Cx,Qx,Px,eta2,zeta); return; //(00|10)
			}else{              ERI1000(Dx,Qx,Px,eta2,zeta); return; //(00|01)
			}
	break;

	case 0: B[0] = 1.0; return;
	break;
	}
	*/

	switch(64*l1+16*l2+4*l3+l4){
	case 64*0+16*0+4*0+0: B[0] = 1.0; return; break;
	
	case 64*1+16*0+4*0+0: ERI1000(Ax,Px,Qx,eta1,zeta); return; break;
	case 64*0+16*1+4*0+0: ERI1000(Bx,Px,Qx,eta1,zeta); return; break;
	case 64*0+16*0+4*1+0: ERI1000(Cx,Qx,Px,eta2,zeta); return; break;
	case 64*0+16*0+4*0+1: ERI1000(Dx,Qx,Px,eta2,zeta); return; break;

	case 64*2+16*0+4*0+0: ERI2000(Ax,Px,Qx,eta1,zeta); return; break;
	case 64*0+16*2+4*0+0: ERI2000(Bx,Px,Qx,eta1,zeta); return; break;
	case 64*0+16*0+4*2+0: ERI2000(Cx,Qx,Px,eta2,zeta); return; break;
	case 64*0+16*0+4*0+2: ERI2000(Dx,Qx,Px,eta2,zeta); return; break;

	case 64*1+16*1+4*0+0: ERI1100(Ax,Bx,Px,Qx,eta1,zeta);      return; break;
	case 64*0+16*0+4*1+1: ERI1100(Cx,Dx,Qx,Px,eta2,zeta);      return; break;
	case 64*1+16*0+4*1+0: ERI1010(Ax,Cx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*0+16*1+4*1+0: ERI1010(Bx,Cx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*1+16*0+4*0+1: ERI1010(Ax,Dx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*0+16*1+4*0+1: ERI1010(Bx,Dx,Px,Qx,eta1,eta2,zeta); return; break;

	case 64*0+16*1+4*1+1: ERI0111(Bx,Cx,Dx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*1+16*0+4*1+1: ERI0111(Ax,Cx,Dx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*1+16*1+4*0+1: ERI0111(Dx,Ax,Bx,Qx,Px,eta2,eta1,zeta); return; break;
	case 64*1+16*1+4*1+0: ERI0111(Cx,Ax,Bx,Qx,Px,eta2,eta1,zeta); return; break;

	case 64*2+16*1+4*0+0: ERI2100(Ax,Bx,Px,Qx,eta1,zeta); return; break;
	case 64*1+16*2+4*0+0: ERI2100(Bx,Ax,Px,Qx,eta1,zeta); return; break;
	case 64*0+16*0+4*2+1: ERI2100(Cx,Dx,Qx,Px,eta2,zeta); return; break;
	case 64*0+16*0+4*1+2: ERI2100(Dx,Cx,Qx,Px,eta2,zeta); return; break;

	case 64*2+16*0+4*1+0: ERI2010(Ax,Cx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*2+16*0+4*0+1: ERI2010(Ax,Dx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*0+16*2+4*1+0: ERI2010(Bx,Cx,Px,Qx,eta1,eta2,zeta); return; break;
	case 64*0+16*2+4*0+1: ERI2010(Bx,Dx,Px,Qx,eta1,eta2,zeta); return; break;

	case 64*1+16*0+4*2+0: ERI2010(Cx,Ax,Qx,Px,eta2,eta1,zeta); return; break;
	case 64*1+16*0+4*0+2: ERI2010(Dx,Ax,Qx,Px,eta2,eta1,zeta); return; break;
	case 64*0+16*1+4*2+0: ERI2010(Cx,Bx,Qx,Px,eta2,eta1,zeta); return; break;
	case 64*0+16*1+4*0+2: ERI2010(Dx,Bx,Qx,Px,eta2,eta1,zeta); return; break;

	case 64*1+16*1+4*1+1: ERI1111(Ax,Bx,Cx,Dx,Px,Qx,eta1,eta2,zeta); return; break;
	}

	PA = Px-Ax;
	PB = Px-Bx;
	QC = Qx-Cx;
	QD = Qx-Dx;
	QP = Qx-Px;

	// prepare binomial coefficients
	f_coef(fi1,l1,l2,PA,PB);
	f_coef(fi2,l3,l4,QC,QD);

	// prepare H matrix
	get_H(HL, fi1, l1, l2, eta1);
	get_H(HM, fi2, l3, l4, eta2);

	// compute CI
	B[0]            = 0.0;
	QP_pow[0]       = 1.0;
	zeta_pow[0]     = 1.0;
	for(I=1; I < l1+l2+l3+l4+1; I++){
		B[I]        = 0.0;
		QP_pow[I]   = QP_pow[I-1]*QP/I;
		zeta_pow[I] = zeta_pow[I-1]*zeta;
	}
	for(L=0; L < l1+l2+1;   L++)
	for(M=0; M < l3+l4+1;   M++)
	for(u=0; u < (L+M)/2+1; u++){
		I    = L+M-u;
		B[I] = B[I] +  HL[L]*HM[M]
		              *((u+M)&1?-1:1)
		              *FACT[L+M]/FACT[u]
		              *QP_pow[L+M-u-u]
		              *zeta_pow[I];
	} 
	return;
}

// eri : compute two electron integral of the form (ab|1/r12|cd).
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation and testing
//
double eri(double xa, double ya, double za, double norma,
           int la, int ma, int na, double alphaa,
           double xb, double yb, double zb, double normb,
           int lb, int mb, int nb, double alphab,
           double xc, double yc, double zc, double normc,
           int lc, int mc, int nc, double alphac,
           double xd, double yd, double zd, double normd,
           int ld, int md, int nd, double alphad){ 

	double Bx[4*MAXL], By[4*MAXL], Bz[4*MAXL], F[4*MAXL];
	double rab2,rcd2,rpq2,xp,yp,zp,xq,yq,zq,sum;
	double eta1, eta2, zeta;
	int kX,kY,kZ;
	double t;

	eta1 = 1.0/(alphaa+alphab);
	eta2 = 1.0/(alphac+alphad);
	zeta = 4.0*(alphaa+alphab)*(alphac+alphad)
	          /(alphaa+alphab+alphac+alphad);

	xp   = (alphaa*xa+alphab*xb)*eta1;
	yp   = (alphaa*ya+alphab*yb)*eta1;
	zp   = (alphaa*za+alphab*zb)*eta1;
	xq   = (alphac*xc+alphad*xd)*eta2;
	yq   = (alphac*yc+alphad*yd)*eta2;
	zq   = (alphac*zc+alphad*zd)*eta2;

	rab2 = DIST2(xa,ya,za,xb,yb,zb);
	rcd2 = DIST2(xc,yc,zc,xd,yd,zd);
	rpq2 = DIST2(xp,yp,zp,xq,yq,zq);

	// compute fgamma
	t     = 0.25*rpq2*zeta;
	fgamma_set(la+lb+lc+ld+ma+mb+mc+md+na+nb+nc+nd, t, F);

	// integrate in 3 directions
	eriB_1d(Bx,la,lb,lc,ld,xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);
	eriB_1d(By,ma,mb,mc,md,yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);
	eriB_1d(Bz,na,nb,nc,nd,zp,za,zb,zq,zc,zd,eta1,eta2,zeta);

	sum = 0.0;
	for(kX=0; kX < la+lb+lc+ld+1; kX++)
	for(kY=0; kY < ma+mb+mc+md+1; kY++)
	for(kZ=0; kZ < na+nb+nc+nd+1; kZ++)
		sum += Bx[kX]*By[kY]*Bz[kZ]*F[kX+kY+kZ];

#define TWO_PI_POW2_5 34.986836655249725
	return  TWO_PI_POW2_5
	        *eta1*eta2
	        /sqrt(alphaa+alphab+alphac+alphad)
	        *exp(-alphaa*alphab*rab2*eta1
	             -alphac*alphad*rcd2*eta2)
	        *sum
	        *norma*normb*normc*normd;
#undef  TWO_PI_POW2_5
}


// primERI : compute two electron integral of the form (ab|1/r12|cd).
//
// Jan 19, 2013 - Teepanis Chachiyo
// 		Initial implementation and testing
//
double primeERI(
	double xa, double ya, double za, int la, int ma, int na, double alphaa,
	double xb, double yb, double zb, int lb, int mb, int nb, double alphab,
	double xc, double yc, double zc, int lc, int mc, int nc, double alphac,
	double xd, double yd, double zd, int ld, int md, int nd, double alphad){ 

	double Bx[4*MAXL], By[4*MAXL], Bz[4*MAXL], F[4*MAXL];
	double rpq2,xp,yp,zp,xq,yq,zq,sum;
	double eta1, eta2, zeta;
	int kX,kY,kZ;
	double t;

	eta1 = 1.0/(alphaa+alphab);
	eta2 = 1.0/(alphac+alphad);
	zeta = 4.0*(alphaa+alphab)*(alphac+alphad)
	          /(alphaa+alphab+alphac+alphad);

	xp   = (alphaa*xa+alphab*xb)*eta1;
	yp   = (alphaa*ya+alphab*yb)*eta1;
	zp   = (alphaa*za+alphab*zb)*eta1;
	xq   = (alphac*xc+alphad*xd)*eta2;
	yq   = (alphac*yc+alphad*yd)*eta2;
	zq   = (alphac*zc+alphad*zd)*eta2;

	rpq2 = DIST2(xp,yp,zp,xq,yq,zq);

	// compute fgamma
	t     = 0.25*rpq2*zeta;
	fgamma_set(la+lb+lc+ld+ma+mb+mc+md+na+nb+nc+nd, t, F);

	// integrate in 3 directions
	eriB_1d(Bx,la,lb,lc,ld,xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);
	eriB_1d(By,ma,mb,mc,md,yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);
	eriB_1d(Bz,na,nb,nc,nd,zp,za,zb,zq,zc,zd,eta1,eta2,zeta);

	sum=0.0;			
	for(kX=0; kX < la+lb+lc+ld+1; kX++)
	for(kY=0; kY < ma+mb+mc+md+1; kY++)
	for(kZ=0; kZ < na+nb+nc+nd+1; kZ++)
		sum += Bx[kX]*By[kY]*Bz[kZ]*F[kX+kY+kZ];

	return  sum/sqrt(alphaa+alphab+alphac+alphad);
}
#undef MAXL

#define MAXL      4
#define MAXL2    16
#define MAXL3    64
#define MAXL4   256
#define TWOMAXL   8
#define FOURMAXL 16
// eri_list : compute multiple two electron integral of the form (ab|1/r12|cd),
//            where all gaussian function must have the same center and
//            exponent, but can have different angular part.
//
//            (*la, *ma, *na, *lb, *mb, *nb, *lc, *mc, *nc, *ld, *md, *nd) are
//            arrays containing the list of all integral with nItem elements.
// 
//            all evaluated integrals are stored in "results"
//
// March 07, 2010 - Teepanis Chachiyo
// 		Initial implementation and testing
//
void  eri_list(double xa, double ya, double za,
              const int *la, const int *ma, const int *na, double alphaa,
              double xb, double yb, double zb,
              const int *lb, const int *mb, const int *nb, double alphab,
              double xc, double yc, double zc,
              const int *lc, const int *mc, const int *nc, double alphac,
              double xd, double yd, double zd,
              const int *ld, const int *md, const int *nd, double alphad,
              int maxSumL, int nItem, double *results){ 

	char BxArray[MAXL4];                       // storage index for Bx
	char ByArray[MAXL4];                       // storage index for By
	char BzArray[MAXL4];                       // storage index for Bz
	double Bx[FOURMAXL*MAXL4];                 // large storage for Bx
	double By[FOURMAXL*MAXL4];                 // large storage for By
	double Bz[FOURMAXL*MAXL4];                 // large storage for Bz
	double *tBx, *tBy, *tBz;                   // auxilary pointers
	double F[4*MAXL];                          // fgamma values
	double rab2,rcd2,rpq2,xp,yp,zp,xq,yq,zq;   // auxilary 
	double eta1, eta2, zeta;                   // auxilary distance
	int kX,kY,kZ;                              // cartesian sum index
	double sum;                                // summation of BxByBzF
	int i,n;                                   // generic index
	double t;                                  // argument of fgamma
	int mX,mY,mZ;                              // max angular moment 
	
	eta1 = 1.0/(alphaa+alphab);
	eta2 = 1.0/(alphac+alphad);
	zeta = 4.0*(alphaa+alphab)*(alphac+alphad)
	          /(alphaa+alphab+alphac+alphad);

	xp   = (alphaa*xa+alphab*xb)*eta1;
	yp   = (alphaa*ya+alphab*yb)*eta1;
	zp   = (alphaa*za+alphab*zb)*eta1;
	xq   = (alphac*xc+alphad*xd)*eta2;
	yq   = (alphac*yc+alphad*yd)*eta2;
	zq   = (alphac*zc+alphad*zd)*eta2;

	rab2 = DIST2(xa,ya,za,xb,yb,zb);
	rcd2 = DIST2(xc,yc,zc,xd,yd,zd);
	rpq2 = DIST2(xp,yp,zp,xq,yq,zq);

	// compute fgamma
	t     = 0.25*rpq2*zeta;	
	fgamma_set(maxSumL, t, F);

	// compute overall factor
#define TWO_PI_POW2_5 34.986836655249725
	t = TWO_PI_POW2_5
	    *eta1*eta2
	    /sqrt(alphaa+alphab+alphac+alphad)
	    *exp(-alphaa*alphab*rab2*eta1
	         -alphac*alphad*rcd2*eta2);
#undef  TWO_PI_POW2_5

	// initialize BxArray, ByArray, BzArray
	memset(BxArray, 0, MAXL4);
	memset(ByArray, 0, MAXL4);
	memset(BzArray, 0, MAXL4);

	// calculate Bx, By, Bz as needed
	for(n=0; n < nItem; n++){

		i = la[n]*MAXL3+lb[n]*MAXL2+lc[n]*MAXL+ld[n];
		tBx = Bx + FOURMAXL*i;
		if(!BxArray[i]){
			eriB_1d(tBx,
			        la[n],lb[n],lc[n],ld[n],
			        xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);
			BxArray[i] = 1;
		}

		i = ma[n]*MAXL3+mb[n]*MAXL2+mc[n]*MAXL+md[n];
		tBy = By + FOURMAXL*i;
		if(!ByArray[i]){
			eriB_1d(tBy,
			        ma[n],mb[n],mc[n],md[n],
			        yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);
			ByArray[i] = 1;
		}


		i = na[n]*MAXL3+nb[n]*MAXL2+nc[n]*MAXL+nd[n];
		tBz = Bz + FOURMAXL*i;
		if(!BzArray[i]){
			eriB_1d(tBz,
			        na[n],nb[n],nc[n],nd[n],
			        zp,za,zb,zq,zc,zd,eta1,eta2,zeta);
			BzArray[i] = 1;
		}

		sum = 0.0;
//		for(kX=0; kX < la[n]+lb[n]+lc[n]+ld[n]+1; kX++)
//		for(kY=0; kY < ma[n]+mb[n]+mc[n]+md[n]+1; kY++)
//		for(kZ=0; kZ < na[n]+nb[n]+nc[n]+nd[n]+1; kZ++)
		mX = la[n]+lb[n]+lc[n]+ld[n]+1;
		mY = ma[n]+mb[n]+mc[n]+md[n]+1;
		mZ = na[n]+nb[n]+nc[n]+nd[n]+1;
		for(kX=0; kX < mX; kX++)
		for(kY=0; kY < mY; kY++)
		for(kZ=0; kZ < mZ; kZ++)
			sum += tBx[kX]*tBy[kY]*tBz[kZ]*F[kX+kY+kZ];

		results[n] = sum*t;
	}
}
#undef MAXL
#undef MAXL2
#undef MAXL3
#undef MAXL4
#undef FOURMAXL

// genSetBxyzF : generate a set of Bx,By,Bz,F according to Taketa scheme
// of evaluating 2e integral. This is efficient because we are doing them
// in batch
//
// May 19, 2011 - Teepanis Chachiyo
//   Use HRR recursive scheme for Bx,By,Bz
//
// July 10, 2010 - Teepanis Chachiyo
//   Initial implementation and testing
//
double  genSetBxyzF(double xa,  double ya,  double za, int maxa, double alphaa,
                    double xb,  double yb,  double zb, int maxb, double alphab,
                    double xc,  double yc,  double zc, int maxc, double alphac,
                    double xd,  double yd,  double zd, int maxd, double alphad,
                    double *Bx, double *By, double *Bz, double *F,int maxL){ 

	double rab2,rcd2,rpq2,xp,yp,zp,xq,yq,zq;   // auxilary 
	double eta1, eta2, zeta;                   // auxilary distance
	int a,b,c,d;                               // angular index
	int i,j,k,n;                               // generic index
	int na,nb,nc;                              // index locator
	double t;                                  // argument of fgamma

	eta1 = 1.0/(alphaa+alphab);
	eta2 = 1.0/(alphac+alphad);
	zeta = 4.0*(alphaa+alphab)*(alphac+alphad)
	          /(alphaa+alphab+alphac+alphad);

	xp   = (alphaa*xa+alphab*xb)*eta1;
	yp   = (alphaa*ya+alphab*yb)*eta1;
	zp   = (alphaa*za+alphab*zb)*eta1;
	xq   = (alphac*xc+alphad*xd)*eta2;
	yq   = (alphac*yc+alphad*yd)*eta2;
	zq   = (alphac*zc+alphad*zd)*eta2;

	rab2 = DIST2(xa,ya,za,xb,yb,zb);
	rcd2 = DIST2(xc,yc,zc,xd,yd,zd);
	rpq2 = DIST2(xp,yp,zp,xq,yq,zq);

	// compute fgamma
	t     = 0.25*rpq2*zeta;
	fgamma_set(maxa+maxb+maxc+maxd, t, F);

	// compute overall factor
#define TWO_PI_POW2_5 34.986836655249725
	t = TWO_PI_POW2_5
	    *eta1*eta2
	    /sqrt(alphaa+alphab+alphac+alphad)
	    *exp(-alphaa*alphab*rab2*eta1
	         -alphac*alphad*rcd2*eta2);
#undef  TWO_PI_POW2_5
	if(t < PRIMITIVE_CUTOFF) return 0.0;

	// generate Bx,By,Bz
	//i = 0;
	//for(a=0; a < maxa+1; a++)
	//for(b=0; b < maxb+1; b++)
	//for(c=0; c < maxc+1; c++)
	//for(d=0; d < maxd+1; d++){
	//	eriB_1d(Bx+i, a, b, c, d, xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);	
	//	eriB_1d(By+i, a, b, c, d, yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);
	//	eriB_1d(Bz+i, a, b, c, d, zp,za,zb,zq,zc,zd,eta1,eta2,zeta);
	//	i += 4*(maxL+1);
	//}

	// generate Bx,By,Bz using recursion
	nc = (maxd+1); nb = nc*(maxc+1); na = nb*(maxb+1);
	for(a=0; a < maxa+1; a++)
	for(b=0; b < maxb+1; b++)
		if(a!=0  && b!=maxb){
			for(c=0; c < maxc+1; c++)
			for(d=0; d < maxd+1; d++){
				i =4*(maxL+1)*(     a*na +     b*nb + c*nc + d);
				j =4*(maxL+1)*( (a-1)*na + (b+1)*nb + c*nc + d);
				k =4*(maxL+1)*( (a-1)*na +     b*nb + c*nc + d);
				n = a+b+c+d;
				Bx[i+n] = Bx[j+n];
				By[i+n] = By[j+n];
				Bz[i+n] = Bz[j+n];
				for(n=n-1; n>=0; n--){
					Bx[i+n] = Bx[j+n] - (xa-xb)*Bx[k+n];
					By[i+n] = By[j+n] - (ya-yb)*By[k+n];
					Bz[i+n] = Bz[j+n] - (za-zb)*Bz[k+n];
				}
			}
		}else

		for(c=0; c < maxc+1; c++)
		for(d=0; d < maxd+1; d++)
			if(c!=0  && d!=maxd){
				i =4*(maxL+1)*( a*na + b*nb +     c*nc +  d   );
				j =4*(maxL+1)*( a*na + b*nb + (c-1)*nc +  d+1 );
				k =4*(maxL+1)*( a*na + b*nb + (c-1)*nc +  d   );
				n = a+b+c+d;
				Bx[i+n] = Bx[j+n];
				By[i+n] = By[j+n];
				Bz[i+n] = Bz[j+n];
				for(n=n-1; n>=0; n--){
					Bx[i+n] = Bx[j+n] - (xc-xd)*Bx[k+n];
					By[i+n] = By[j+n] - (yc-yd)*By[k+n];
					Bz[i+n] = Bz[j+n] - (zc-zd)*Bz[k+n];
				}
			}else{
				i =4*(maxL+1)*(a*na + b*nb + c*nc + d);
				eriB_1d(Bx+i, a, b, c, d, xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);	
				eriB_1d(By+i, a, b, c, d, yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);
				eriB_1d(Bz+i, a, b, c, d, zp,za,zb,zq,zc,zd,eta1,eta2,zeta);
			}

	return t;
}


// genSetBxyzSxyzF : generates a set of Bx,By,Bz,F. In addition, it produces
// the Sx,Sy,Sz summation array
//
// Jan 23, 2013 - Teepanis Chachiyo
//   Initial implementation and testing
//
double  genSetBxyzSxyzF(
	double xa,  double ya,  double za, int maxa, double alphaa,
	double xb,  double yb,  double zb, int maxb, double alphab,
	double xc,  double yc,  double zc, int maxc, double alphac,
	double xd,  double yd,  double zd, int maxd, double alphad,
	double *Bx, double *By, double *Bz,
	double *Sx, double *Sy, double *Sz,
	double *F,int maxL){ 

	double rab2,rcd2,rpq2,xp,yp,zp,xq,yq,zq;   // auxilary 
	double eta1, eta2, zeta;                   // auxilary distance
	int a,b,c,d;                               // angular index
	int i,j,k,n;                               // generic index
	int na,nb,nc;                              // index locator
	double t;                                  // argument of fgamma
	register double xSum,ySum,zSum;

	eta1 = 1.0/(alphaa+alphab);
	eta2 = 1.0/(alphac+alphad);
	zeta = 4.0*(alphaa+alphab)*(alphac+alphad)
	          /(alphaa+alphab+alphac+alphad);

	xp   = (alphaa*xa+alphab*xb)*eta1;
	yp   = (alphaa*ya+alphab*yb)*eta1;
	zp   = (alphaa*za+alphab*zb)*eta1;
	xq   = (alphac*xc+alphad*xd)*eta2;
	yq   = (alphac*yc+alphad*yd)*eta2;
	zq   = (alphac*zc+alphad*zd)*eta2;

	rab2 = DIST2(xa,ya,za,xb,yb,zb);
	rcd2 = DIST2(xc,yc,zc,xd,yd,zd);
	rpq2 = DIST2(xp,yp,zp,xq,yq,zq);

	// compute fgamma
	t     = 0.25*rpq2*zeta;
	fgamma_set(maxa+maxb+maxc+maxd, t, F);

	// compute overall factor
#define TWO_PI_POW2_5 34.986836655249725
	t = TWO_PI_POW2_5
	    *eta1*eta2
	    /sqrt(alphaa+alphab+alphac+alphad)
	    *exp(-alphaa*alphab*rab2*eta1
	         -alphac*alphad*rcd2*eta2);
#undef  TWO_PI_POW2_5
	if(t < PRIMITIVE_CUTOFF) return 0.0;

	// generate Bx,By,Bz using recursion
	nc = (maxd+1); nb = nc*(maxc+1); na = nb*(maxb+1);
	for(a=0; a < maxa+1; a++)
	for(b=0; b < maxb+1; b++)
		if(a!=0  && b!=maxb){
			for(c=0; c < maxc+1; c++)
			for(d=0; d < maxd+1; d++){
				i =4*(maxL+1)*(     a*na +     b*nb + c*nc + d);
				j =4*(maxL+1)*( (a-1)*na + (b+1)*nb + c*nc + d);
				k =4*(maxL+1)*( (a-1)*na +     b*nb + c*nc + d);
				n = a+b+c+d;
				Bx[i+n] = Bx[j+n];
				By[i+n] = By[j+n];
				Bz[i+n] = Bz[j+n];
				for(n=n-1; n>=0; n--){
					Bx[i+n] = Bx[j+n] - (xa-xb)*Bx[k+n];
					By[i+n] = By[j+n] - (ya-yb)*By[k+n];
					Bz[i+n] = Bz[j+n] - (za-zb)*Bz[k+n];
				}

				n = a+b+c+d;
				for(k=maxa+maxb+maxc+maxd-n; k>=0; k--){
					xSum = ySum = zSum = 0.0;
					for(j=n;j>=0;j--){
						xSum+=Bx[i+j]*F[k+j];
						ySum+=By[i+j]*F[k+j];
						zSum+=Bz[i+j]*F[k+j];
					}
					Sx[i+k] = xSum;
					Sy[i+k] = ySum;
					Sz[i+k] = zSum;
				}

			}
		}else

		for(c=0; c < maxc+1; c++)
		for(d=0; d < maxd+1; d++)
			if(c!=0  && d!=maxd){
				i =4*(maxL+1)*( a*na + b*nb +     c*nc +  d   );
				j =4*(maxL+1)*( a*na + b*nb + (c-1)*nc +  d+1 );
				k =4*(maxL+1)*( a*na + b*nb + (c-1)*nc +  d   );
				n = a+b+c+d;
				Bx[i+n] = Bx[j+n];
				By[i+n] = By[j+n];
				Bz[i+n] = Bz[j+n];
				for(n=n-1; n>=0; n--){
					Bx[i+n] = Bx[j+n] - (xc-xd)*Bx[k+n];
					By[i+n] = By[j+n] - (yc-yd)*By[k+n];
					Bz[i+n] = Bz[j+n] - (zc-zd)*Bz[k+n];
				}

				n = a+b+c+d;
				for(k=maxa+maxb+maxc+maxd-n; k>=0; k--){
					xSum = ySum = zSum = 0.0;
					for(j=n;j>=0;j--){
						xSum+=Bx[i+j]*F[k+j];
						ySum+=By[i+j]*F[k+j];
						zSum+=Bz[i+j]*F[k+j];
					}
					Sx[i+k] = xSum;
					Sy[i+k] = ySum;
					Sz[i+k] = zSum;
				}

			}else{
				i =4*(maxL+1)*(a*na + b*nb + c*nc + d);
				eriB_1d(Bx+i, a, b, c, d, xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);	
				eriB_1d(By+i, a, b, c, d, yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);
				eriB_1d(Bz+i, a, b, c, d, zp,za,zb,zq,zc,zd,eta1,eta2,zeta);

				n = a+b+c+d;
				for(k=maxa+maxb+maxc+maxd-n; k>=0; k--){
					xSum = ySum = zSum = 0.0;
					for(j=n;j>=0;j--){
						xSum+=Bx[i+j]*F[k+j];
						ySum+=By[i+j]*F[k+j];
						zSum+=Bz[i+j]*F[k+j];
					}
					Sx[i+k] = xSum;
					Sy[i+k] = ySum;
					Sz[i+k] = zSum;
				}

			}

	return t;
}


#define MAXL 4
// genPrimeSetBxyzSxyzF : generates a set of Bx,By,Bz,F. In addition, it produces
// the Sx,Sy,Sz summation array. This is tailored for primitive space
//
// Jan 22, 2013 - Teepanis Chachiyo
//   Initial implementation and testing
//
void  genPrimeSetBxyzSxyzF(
	double xa,  double ya,  double za, int maxa, double alphaa,
	double xb,  double yb,  double zb, int maxb, double alphab,
	double xc,  double yc,  double zc, int maxc, double alphac,
	double xd,  double yd,  double zd, int maxd, double alphad,
	double *Bx, double *By, double *Bz,
	double *Sx, double *Sy, double *Sz, double *F){ 

	double rpq2,xp,yp,zp,xq,yq,zq;   // auxilary 
	double eta1, eta2, zeta;         // auxilary distance
	int a,b,c,d;                     // angular index
	int i,j,k,n;                     // generic index
	int na,nb,nc;                    // index locator
	double t;                        // argument of fgamma
	double xSum,ySum,zSum;

	eta1 = 1.0/(alphaa+alphab);
	eta2 = 1.0/(alphac+alphad);
	zeta = 4.0*(alphaa+alphab)*(alphac+alphad)
	          /(alphaa+alphab+alphac+alphad);

	xp   = (alphaa*xa+alphab*xb)*eta1;
	yp   = (alphaa*ya+alphab*yb)*eta1;
	zp   = (alphaa*za+alphab*zb)*eta1;
	xq   = (alphac*xc+alphad*xd)*eta2;
	yq   = (alphac*yc+alphad*yd)*eta2;
	zq   = (alphac*zc+alphad*zd)*eta2;
	rpq2 = DIST2(xp,yp,zp,xq,yq,zq);

	// compute fgamma
	t     = 0.25*rpq2*zeta;
	fgamma_set(maxa+maxb+maxc+maxd, t, F);

	// generate Bx,By,Bz using recursion
	nc = (maxd+1); nb = nc*(maxc+1); na = nb*(maxb+1);
	for(a=0; a < maxa+1; a++)
	for(b=0; b < maxb+1; b++)
		if(a!=0  && b!=maxb){
			for(c=0; c < maxc+1; c++)
			for(d=0; d < maxd+1; d++){
				i =4*(MAXL+1)*(     a*na +     b*nb + c*nc + d);
				j =4*(MAXL+1)*( (a-1)*na + (b+1)*nb + c*nc + d);
				k =4*(MAXL+1)*( (a-1)*na +     b*nb + c*nc + d);
				n = a+b+c+d;
				Bx[i+n] = Bx[j+n];
				By[i+n] = By[j+n];
				Bz[i+n] = Bz[j+n];
				for(n=n-1; n>=0; n--){
					Bx[i+n] = Bx[j+n] - (xa-xb)*Bx[k+n];
					By[i+n] = By[j+n] - (ya-yb)*By[k+n];
					Bz[i+n] = Bz[j+n] - (za-zb)*Bz[k+n];
				}

				n = a+b+c+d;
				for(k=maxa+maxb+maxc+maxd-n; k>=0; k--){
					xSum = ySum = zSum = 0.0;
					for(j=n;j>=0;j--){
						xSum+=Bx[i+j]*F[k+j];
						ySum+=By[i+j]*F[k+j];
						zSum+=Bz[i+j]*F[k+j];
					}
					Sx[i+k] = xSum;
					Sy[i+k] = ySum;
					Sz[i+k] = zSum;
				}

			}
		}else

		for(c=0; c < maxc+1; c++)
		for(d=0; d < maxd+1; d++)
			if(c!=0  && d!=maxd){
				i =4*(MAXL+1)*( a*na + b*nb +     c*nc +  d   );
				j =4*(MAXL+1)*( a*na + b*nb + (c-1)*nc +  d+1 );
				k =4*(MAXL+1)*( a*na + b*nb + (c-1)*nc +  d   );
				n = a+b+c+d;
				Bx[i+n] = Bx[j+n];
				By[i+n] = By[j+n];
				Bz[i+n] = Bz[j+n];
				for(n=n-1; n>=0; n--){
					Bx[i+n] = Bx[j+n] - (xc-xd)*Bx[k+n];
					By[i+n] = By[j+n] - (yc-yd)*By[k+n];
					Bz[i+n] = Bz[j+n] - (zc-zd)*Bz[k+n];
				}

				n = a+b+c+d;
				for(k=maxa+maxb+maxc+maxd-n; k>=0; k--){
					xSum = ySum = zSum = 0.0;
					for(j=n;j>=0;j--){
						xSum+=Bx[i+j]*F[k+j];
						ySum+=By[i+j]*F[k+j];
						zSum+=Bz[i+j]*F[k+j];
					}
					Sx[i+k] = xSum;
					Sy[i+k] = ySum;
					Sz[i+k] = zSum;
				}

			}else{
				i =4*(MAXL+1)*(a*na + b*nb + c*nc + d);
				eriB_1d(Bx+i, a, b, c, d, xp,xa,xb,xq,xc,xd,eta1,eta2,zeta);	
				eriB_1d(By+i, a, b, c, d, yp,ya,yb,yq,yc,yd,eta1,eta2,zeta);
				eriB_1d(Bz+i, a, b, c, d, zp,za,zb,zq,zc,zd,eta1,eta2,zeta);

				n = a+b+c+d;
				for(k=maxa+maxb+maxc+maxd-n; k>=0; k--){
					xSum = ySum = zSum = 0.0;
					for(j=n;j>=0;j--){
						xSum+=Bx[i+j]*F[k+j];
						ySum+=By[i+j]*F[k+j];
						zSum+=Bz[i+j]*F[k+j];
					}
					Sx[i+k] = xSum;
					Sy[i+k] = ySum;
					Sz[i+k] = zSum;
				}


			}
}
#undef MAXL

// contr_ssss : compute 2electron integral for ss|ss type
//
// July 3, 2010 - Teepanis Chachiyo
//                Implementation but need testing
//
// May 19, 2011 - Teepanis Chachiyo
//                Take rab2 and rcd2 out of the loop and use fgamma_0 
//                instead of fgamma_set
//
double contr_ssss(
        int lena,const double *aexps,const double *acoefs,const double *anorms,
        double xa,double ya,double za,
        int lenb,const double *bexps,const double *bcoefs,const double *bnorms,
        double xb,double yb,double zb,
        int lenc,const double *cexps,const double *ccoefs,const double *cnorms,
        double xc,double yc,double zc,
        int lend,const double *dexps,const double *dcoefs,const double *dnorms,
        double xd,double yd,double zd){
			 
	double val = 0.;
	int i,j,k,l;
	double EE, F; 

	double rab2,rcd2,rpq2,xp,yp,zp,xq,yq,zq;
	double eta1, eta2, zeta;
	double t;

	rab2 = DIST2(xa,ya,za,xb,yb,zb);
	rcd2 = DIST2(xc,yc,zc,xd,yd,zd);

	// proceed from highest exponent value
	for (i=0; i<lena; i++) for (j=0; j<lenb; j++){

			eta1 = 1.0/(aexps[i]+bexps[j]);
			xp   = (aexps[i]*xa+bexps[j]*xb)*eta1;
			yp   = (aexps[i]*ya+bexps[j]*yb)*eta1;
			zp   = (aexps[i]*za+bexps[j]*zb)*eta1;
			
	for (k=0; k<lenc; k++) for (l=0; l<lend; l++){

			eta2 = 1.0/(cexps[k]+dexps[l]);
			xq   = (cexps[k]*xc+dexps[l]*xd)*eta2;
			yq   = (cexps[k]*yc+dexps[l]*yd)*eta2;
			zq   = (cexps[k]*zc+dexps[l]*zd)*eta2;
			
			rpq2 = DIST2(xp,yp,zp,xq,yq,zq);	
			zeta = 4.0*(aexps[i]+bexps[j])*(cexps[k]+dexps[l])
	                  /(aexps[i]+bexps[j]+cexps[k]+dexps[l]);

			// compute fgamma
			t     = 0.25*rpq2*zeta;
			F     = fgamma_0(t);

		// compute element
#define TWO_PI_POW2_5 34.986836655249725
			EE = TWO_PI_POW2_5
			     *F
	             *eta1*eta2
	             /sqrt(aexps[i]+bexps[j]+cexps[k]+dexps[l])
	             *exp(-aexps[i]*bexps[j]*rab2*eta1
	                  -cexps[k]*dexps[l]*rcd2*eta2)
	             *acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]
	             *anorms[i]*bnorms[j]*cnorms[k]*dnorms[l];
#undef  TWO_PI_POW2_5
		val += EE;
	}
	}

	return val;
}

