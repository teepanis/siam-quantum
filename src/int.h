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

#ifndef INT_H
#define INT_H
#define PRIMITIVE_CUTOFF 1.0E-15
double overlap(double alpha1, int l1, int m1, int n1,
               double xa, double ya, double za,
               double alpha2, int l2, int m2, int n2,
               double xb, double yb, double zb);

double moment(
	double alpha1, int l1, int m1, int n1, // gto(1) exponent and angular index
	double xa, double ya, double za,       // gto(1) center
	double alpha2, int l2, int m2, int n2, // gto(2) exponent and angular index
	double xb, double yb, double zb,       // gto(2) center
	int mx, int my, int mz,                // moment angular index
	double xc, double yc, double zc);      // moment center

double kinetic(double alpha1, int l1, int m1, int n1,
               double xa, double ya, double za,
               double alpha2, int l2, int m2, int n2,
               double xb, double yb, double zb);

double nai(double x1, double y1, double z1, double norm1,
           int l1, int m1, int n1, double alpha1,
           double x2, double y2, double z2, double norm2,
           int l2, int m2, int n2, double alpha2,
           double x3, double y3, double z3);

void naiA_1d(double *A,  int l1,     int l2,
             double PAx, double PBx, double CPx,
             double eta);

void efi(double x1, double y1, double z1, double norm1,
         int l1, int m1, int n1, double alpha1,
         double x2, double y2, double z2, double norm2,
         int l2, int m2, int n2, double alpha2,
         double x3, double y3, double z3,
         double *ex, double *ey, double *ez);

void eriB_1d(double *B, int l1, int l2, int l3, int l4,
                    double Px, double Ax, double Bx,
                    double Qx, double Cx, double Dx,
                    double eta1, double eta2,
                    double zeta);

double eri(double xa, double ya, double za, double norma,
           int la, int ma, int na, double alphaa,
           double xb, double yb, double zb, double normb,
           int lb, int mb, int nb, double alphab,
           double xc, double yc, double zc, double normc,
           int lc, int mc, int nc, double alphac,
           double xd, double yd, double zd, double normd,
           int ld, int md, int nd, double alphad);

double primeERI(
	double xa, double ya, double za, int la, int ma, int na, double alphaa,
	double xb, double yb, double zb, int lb, int mb, int nb, double alphab,
	double xc, double yc, double zc, int lc, int mc, int nc, double alphac,
	double xd, double yd, double zd, int ld, int md, int nd, double alphad);

void  eri_list(double xa, double ya, double za,
              const int *la, const int *ma, const int *na, double alphaa,
              double xb, double yb, double zb,
              const int *lb, const int *mb, const int *nb, double alphab,
              double xc, double yc, double zc,
              const int *lc, const int *mc, const int *nc, double alphac,
              double xd, double yd, double zd,
              const int *ld, const int *md, const int *nd, double alphad,
              int maxSumL, int nItem, double *results);

double  genSetBxyzF(double xa,  double ya,  double za, int maxa, double alphaa,
                    double xb,  double yb,  double zb, int maxb, double alphab,
                    double xc,  double yc,  double zc, int maxc, double alphac,
                    double xd,  double yd,  double zd, int maxd, double alphad,
                    double *Bx, double *By, double *Bz, double *F,int maxL);

double  genSetBxyzSxyzF(
	double xa,  double ya,  double za, int maxa, double alphaa,
	double xb,  double yb,  double zb, int maxb, double alphab,
	double xc,  double yc,  double zc, int maxc, double alphac,
	double xd,  double yd,  double zd, int maxd, double alphad,
	double *Bx, double *By, double *Bz,
	double *Sx, double *Sy, double *Sz,
	double *F,int maxL);

void  genPrimeSetBxyzSxyzF(
	double xa,  double ya,  double za, int maxa, double alphaa,
	double xb,  double yb,  double zb, int maxb, double alphab,
	double xc,  double yc,  double zc, int maxc, double alphac,
	double xd,  double yd,  double zd, int maxd, double alphad,
	double *Bx, double *By, double *Bz,
	double *Sx, double *Sy, double *Sz, double *F);

double contr_ssss(
        int lena,const double *aexps,const double *acoefs,const double *anorms,
        double xa,double ya,double za,
        int lenb,const double *bexps,const double *bcoefs,const double *bnorms,
        double xb,double yb,double zb,
        int lenc,const double *cexps,const double *ccoefs,const double *cnorms,
        double xc,double yc,double zc,
        int lend,const double *dexps,const double *dcoefs,const double *dnorms,
        double xd,double yd,double zd);
#endif
