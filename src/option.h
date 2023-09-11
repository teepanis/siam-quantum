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

#ifndef OPTION_H
#define OPTION_H
struct option_t{
int molCharge;        // total molecular charge
int multiplicity;     // molecule spin multiplicity = 2*s + 1

int outVolumeGrid;    // number of points per Angstrom
int outVolumeType;    // type of volume information to compute and save
#define VOLUME_NONE          0
#define VOLUME_DENSITY_TOTAL 1
#define VOLUME_DENSITY_SPIN  2
#define VOLUME_MO_ALPHA      3
#define VOLUME_MO_BETA       4
#define VOLUME_POTENTIAL     5
#define VOLUME_GRADRHO_OVER_RHO43   6
#define VOLUME_GRAD_RS       7

int outWhichMO;       // molecular orbital index to compute volume information
int outFormat;        // file format to save the volume information to
#define VOLUME_CUBE   1
#define VOLUME_XSF    2

double outVolumeCut;  // cutoff accurary for computing volume information
int outGAUSSIAN;      // flag to export output in Gaussian log format

double SCFConv;       // SCF convergence threshold
double SCFDrag;       // SCF drag coefficient
double SCFCutoff;     // SCF integral schwarz cut-off
int SCFMax;           // maximum number of scf cycle
#define CONVMETHOD_DAMPING   0
#define CONVMETHOD_DIIS2     1
#define CONVMETHOD_DIIS3     2
#define CONVMETHOD_DIIS4     3
int convMethod;       // scf convergence method
int SCFAccuracy;      // how to handle integral accuracy during scf
#define SCFACCURACY_3STEP 0
#define SCFACCURACY_1STEP 1

#define SCFACCURACY_3STEP_A 1.0E-3
#define SCFACCURACY_3STEP_B 1.0E-6
#define SCFACCURACY_3STEP_C 1.0E-9

int SCFGuess;         // guess method
#define SCFGUESS_DIAG  0  // use identity density matrix
#define SCFGUESS_CORE  1  // diagonalize core harmiltonian
#define SCFGUESS_CHECK 2  // read from check point file
#define SCFGUESS_CACB  3  // use the molecular coefficients provided

int maxMem;            // maximum memory in Mbytes to use as storage 
int loadDMatrix;       // flag to load density matrix
int saveDMatrix;       // flag to save density matrix
int saveCheck;         // flag to save checkpoint file
int saveCheckAll;      // flag to save checkpoint file every scf cycle
int loadCheck;         // flag to read molecular orbital only without SCF
int opt;               // float to request geometry optimization
int optMax;            // maximum number of geometry optimization cycle

int restricted;        // flag to use the same MO for alpha and beta
int method;            // method of calculation
#define METHOD_HF         0 // Hartree-Fock (default)
#define METHOD_MP2        1 // Moller-Pleasset perturbation
#define METHOD_SVWN5      2 // Slater exchange VWN5 correlation
#define METHOD_HALF       3 // 0.5*HF + 0.5*Slater + VWN5
#define METHOD_SCHACHIYO  4 // Slater exchange Chachiyo uniform correlation
#define METHOD_HF2xDFT    5 // DFT exchange calculated from Hartree-Fock Densities
int whichHF2xDFT;           // which exchange functional to calculate
#define HF2xDFT_xSLATER          0 // Dirac exchange
#define HF2xDFT_xBECKE88         1 // Becke 1988 exchange
#define HF2xDFT_xPY86            2 // J.P.Perdew and W.Yue 1986 exchange
#define HF2xDFT_xPBE             3 // Perdew-Burke-Ernzerhof    exchange
#define HF2xDFT_xCHACHIYO        4 // T.Chachiyo and H.Chachiyo exchange
#define HF2xDFT_xMVS             5 // meta-GGA Made Very Simple exchange

#define METHOD_xBECKE88   6 // Only Becke88  exchange without correlation
#define METHOD_CHACHIYO   7 // Chachiyo exchange and Chachiyo correlation
#ifdef LIBXC
#define METHOD_LIBXC      8 // Using exchange and correlation from LibXC
int func1Id;                // functional 1 id in xc_funcs.h (0 if not active)
int func2Id;                // functional 2 id in xc_funcs.h (0 of not active)
double hyb_hf_coef;         // hartree-fock coefficent for hybrid functional
#endif

int whichGrid;           // which grid partitioning method
#define GRID_BALLS 1     // Becke partitioning + uniform radius
#define GRID_BECKE 2     // Becke partitioning + bragg atomic radius
int nradial;             // number of radial point
int nlabedev;            // labedev order
int nradialInit;         // number of radial point for early SCF cycles
int nlabedevInit;        // labedev order for early SCF cycles

//
// METHOD_LDA is equivalent to the SVWN5 in Gaussian 09
//
// METHOD_HALF is equivalent to the Gaussian 09 route section:
// SVWN5/<Basis> IOp(3/76=1000005000) IOp(3/77=0000005000) IOp(3/78=0000010000)

int force;             // flag to compute force acting on nuclei
int hessian;           // flag to compute hessian matrix
char DMatrixFile[256]; // file name to process density matrix
char CheckFile[256];   // file name to process checkpoint file
int nCPU;              // number of cpu to perform parallel calculations
char prefixStr[256];   // prefix string for job identification
int MECP;              // flag to request MECP
int mecpMax;           // maximum number of MECP cycle
int mecpMA;            // multiplicity of the first state
int mecpMB;            // multiplicity of the second state
char DMatrixFileA[256];// file name to process density matrix of the first state
char DMatrixFileB[256];// file name to process density matrix of the second state
char CheckFileA[256];  // file name to process checkpoint file
char CheckFileB[256];  // file name to process checkpoint file 
char gaussEXE[256];    // Gaussian program execution string
char gaussINA[256];    // Gaussian input file for the first state
char gaussINB[256];    // Gaussian input file for the second state

int excite;               // flag to compute excited state
#define EXCITED_STATES 5  // number of excited states by default

double Ex,Ey,Ez;          // electric field vector in AU

int qmd;                  // quantum molecular dynamics flag
double initTemp;          // generate initial velocity with the temperature
double keepTemp;          // maintain temperature 
char initVel[256];        // file name to read initial velocity (nm/ps)
char traj[256];           // file name to save trajectory [angstrom]
double qmdFreq;           // frequency of electric field in THz
double tStep;             // time step in pico-sec
int maxQMDStep;           // maximum number of qmd steps
};

void parse_option(struct option_t *opt, int argc, char *argv[]);
void option_help();
#endif

