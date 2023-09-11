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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

// time_str : puts the formatted time into the string
// str. It returns the number of characters written.
//
// Feb 20, 2008 - Teepanis Chachiyo
//     Initial implementation.
#define TIME_FORMAT "%b %d, %Y - %H:%M:%S"
int time_str(int max, char *str){
	time_t t;
	struct tm *tmp;

	t = time(NULL);
	tmp = localtime(&t);
	if(tmp == NULL){
		perror("localtime");
		exit(EXIT_FAILURE);
	}

	return strftime(str, max, TIME_FORMAT, tmp);
}


// Summary:    Search for a set of "n" string str.
//             If one of the string given is "*", then it will
//             match anything (i.e. wildcard).
//
// Parameters: fd  - file pointer to search
//             n   - number of search strings
//             str - list of keywords
//
// Return:     EOF if the end of file
//
// Author:     Teepanis Chachiyo
//
// Created:    08/02/2005
//
// Modified:
//
int findf(FILE *fd, int n, char *str, ...){

  int  s;
  int  ret;
  char buf[1024];

  for(s = 0; (ret=fscanf(fd, "%s", buf)) == 1; ){

    if( strncmp(buf, *(&str+s), 1024) == 0 ||
        strcmp("*", *(&str+s)) == 0) s++; else s = 0;

    if( s == n ) break;
  }

  return ret;
}

struct AtomName_t{
	int   Z;
	char  *shortName;
	char  *longName;
};

//
// Sep 5, 2010 - Teepanis Chachiyo
//      Fix type from Phosphorus to Phosphorous
//
// Nov 4, 2009 - Theerapon Kumla
//		Fix typo from Selinium to Selenium
//		Add atomic number 55 - 103
//
#define MAX_PERIODIC_ATOM 103
static struct AtomName_t PeriodicName[MAX_PERIODIC_ATOM] =
{
              {1, "H", "HYDROGEN"},  \
              {2, "HE","HELIUM"},    \
              {3, "LI","LITHIUM"},   \
              {4, "BE","BERYLLIUM"}, \
              {5, "B", "BORON"},     \
              {6, "C", "CARBON"},    \
              {7, "N", "NITROGEN"},  \
              {8, "O", "OXYGEN"},    \
              {9, "F", "FLUORINE"},  \
              {10,"NE","NEON"},      \
              {11,"NA","SODIUM"},    \
              {12,"MG","MAGNESIUM"}, \
              {13,"AL","ALUMINUM"},  \
              {14,"SI","SILICON"},   \
              {15,"P", "PHOSPHOROUS"},\
              {16,"S", "SULFUR"},    \
              {17,"CL","CHLORINE"},  \
              {18,"AR","ARGON"},     \
              {19,"K", "POTASSIUM"}, \
              {20,"CA","CALCIUM"},   \
              {21,"SC","SCANDIUM"},  \
              {22,"TI","TITANIUM"},  \
              {23,"V", "VANADIUM"},  \
              {24,"CR","CHROMIUM"},  \
              {25,"MN","MANGANESE"}, \
              {26,"FE","IRON"},      \
              {27,"CO","COBALT"},    \
              {28,"NI","NICKEL"},    \
              {29,"CU","COPPER"},    \
              {30,"ZN","ZINC"},      \
              {31,"GA","GALLIUM"},   \
              {32,"GE","GERMANIUM"}, \
              {33,"AS","ARSENIC"},   \
              {34,"SE","SELENIUM"},  \
              {35,"BR","BROMINE"},   \
              {36,"KR","KRYPTON"},   \
              {37,"RB","RUBIDIUM"},  \
              {38,"SR","STRONTIUM"}, \
              {39,"Y", "YTTRIUM"},   \
              {40,"ZR","ZIRCONIUM"}, \
              {41,"NB","NIOBIUM"},   \
              {42,"MO","MOLYBDENUM"},\
              {43,"TC","TECHNETIUM"},\
              {44,"RU","RUTHENIUM"}, \
              {45,"RH","RHODIUM"},   \
              {46,"PD","PALLADIUM"}, \
              {47,"AG","SILVER"},    \
              {48,"CD","CADMIUM"},   \
              {49,"IN","INDIUM"},    \
              {50,"SN","TIN"},       \
              {51,"SB","ANTIMONY"},  \
              {52,"TE","TELLURIUM"}, \
              {53,"I", "IODINE"},    \
              {54,"XE","XENON"},     \
              {55,"CS","CESIUM"},    \
              {56,"BA","BARIUM"},    \
              {57,"LA","LANTHANUM"}, \
              {58,"CE","CERIUM"},    \
              {59,"PR","PRASEODYMIUM"}, \
              {60,"ND","NEODYMIUM"},    \
              {61,"PM","PROMETHIUM"},   \
              {62,"SM","SAMARIUM"},     \
              {63,"EU","EUROPIUM"},     \
              {64,"GD","GADOLIUM"},     \
              {65,"TB","TERBIUM"},      \
              {66,"DY","DYSPROSIUM"},   \
              {67,"HO","HOLMIUM"},      \
              {68,"ER","ERBIUM"},       \
              {69,"TM","THULIUM"},      \
              {70,"YB","YTTERBIUM"},    \
              {71,"LU","LUTETIUM"},     \
              {72,"HF","HAFNIUM"},      \
              {73,"TA","TANTALUM"},     \
              {74,"W","TUNGSTEN"},      \
              {75,"RE","RHENIUM"},      \
              {76,"OS","OSMIUM"},       \
              {77,"IR","IRIDIUM"},      \
              {78,"PT","PLATINUM"},     \
              {79,"AU","GOLD"},         \
              {80,"HG","MERCURY"},      \
              {81,"TL","THALLIUM"},     \
              {82,"PB","LEAD"},         \
              {83,"BI","BISMUTH"},      \
              {84,"PO","POLONIUM"},     \
              {85,"AT","ASTATINE"},     \
              {86,"RN","RADON"},        \
              {87,"FR","FRANCIUM"},     \
              {88,"RA","RADIUM"},       \
              {89,"AC","ACTINIUM"},     \
              {90,"TH","THORIUM"},      \
              {91,"PA","PROTACTINIUM"}, \
              {92,"U","URANIUM"},       \
              {93,"NP","NEPTUNIUM"},    \
              {94,"PU","PLUTONIUM"},    \
              {95,"AM","AMERICIUM"},    \
              {96,"CM","CURIUM"},       \
              {97,"BK","BERKELIUM"},    \
              {98,"CF","CALIFORNIUM"},  \
              {99,"ES","EINSTEINIUM"},  \
              {100,"FM","FERMIUM"},     \
              {101,"MD","MENDELEVIUM"}, \
              {102,"NO","NOBELIUM"},    \
              {103,"LR","LAWRENCIUM"}   
};

// sym2Z : convert atom name which could be in both short 
// or long format into atomic number. Any name that does not
// match string will be Z=0
//
// Feb 26, 2008 - Teepanis Chachiyo
//     Initial implementation
//
#ifndef SYMB_SHORTNAME
#define SYMB_SHORTNAME 0
#endif
#ifndef SYMB_LONGNAME
#define SYMB_LONGNAME  1
#endif
int sym2Z(const char *sym, int type){
	char str[1024];
	int i;

	// make sym upper case
	for(i=0; i < strlen(sym); i++)
		str[i] = toupper(sym[i]);
	// null termination
	str[i] = '\0';

	for(i=0; i < MAX_PERIODIC_ATOM; i++){
		switch(type){
		case SYMB_SHORTNAME:
			if(strcmp(str,PeriodicName[i].shortName)==0)
				return PeriodicName[i].Z;
			break;
		case SYMB_LONGNAME:
			if(strcmp(str,PeriodicName[i].longName)==0)
				return PeriodicName[i].Z;
			break;
		}
	}
	return 0;
}


// Summary:    Convert from atomic number to standard atom
//             name defined in PeriodicName[..] above
//
// Parameters: Z   - atomic number
//             sym - pointer to recieve atom name
//
// Author:     Teepanis Chachiyo
//
// Created:    2005
//
// Modified:
//
//	Nov 4, 2009 - Theerapon Kumla
//		catch error if name not found
//
void Z2Sym(int Z, char *sym){
	if(0 < Z && Z <= MAX_PERIODIC_ATOM)
		strcpy(sym, PeriodicName[Z-1].longName);
	else{
		printf("Z2sym - Error : Z = %d\n",Z);
		exit(EXIT_FAILURE);
	}
}


// Summary:    Convert from atomic number to standard atom
//             name defined in PeriodicName[..] above
//
// Parameters: Z   - atomic number
//             sym - pointer to recieve atom name
//
// Author:     Teepanis Chachiyo
//
// Created:    2010
//
// Modified:
//
//
void Z2SymShort(int Z, char *sym){
	if(0 < Z && Z <= MAX_PERIODIC_ATOM)
		strcpy(sym, PeriodicName[Z-1].shortName);
	else{
		printf("Z2sym - Error : Z = %d\n",Z);
		exit(EXIT_FAILURE);
	}
}


// pow_int : compute x^i where i is integer. This aims to
// minimize the number of multiplication as much as possible.
//
// Feb 17, 2008 - Teepanis Chachiyo
//    Original implementation.
//
double pow_int(double x,int i){
	double t;

	switch(i){
		case -4: t = x*x;    return 1/t/t;
		case -3: t = x*x*x;  return 1/t;
		case -2: t = x*x;    return 1/t;
		case -1:             return 1/x;
		case  0:             return 1.0;
		case  1:             return x;
		case  2: t = x*x;    return t;
		case  3: t = x*x;    return t*x;
		case  4: t = x*x;    return t*t;
		case  5: t = x*x;    return t*t*x;
		case  6: t = x*x*x;  return t*t;
		case  7: t = x*x*x;  return t*t*x;
		case  8: t = x*x;    return t*t*t*t;
		case  9: t = x*x*x;  return t*t*t;
		case 10: t = x*x*x;  return t*t*t*x;
		default:             return pow(x,i);
	}
}


// stepProgressBar: increase progress step by one
// This function increment the "status" index by one.
// and print out status report if needed.
// 
// maxbar : maximum number of bar 
// maxn   : maximum number of item
//
// call with maxbar=0 to reset 
//
// March 9, 2016 - Teepanis Chachiyo
//    Initial implementation and testing
//
void stepProgressBar(int maxbar, int maxn){
	static int bar=0; // previous iteration bar item
	static int n=0;   // previous iteration     item

	// trap reset call
	if(maxbar==0){ bar=0; n=0; return; }

	// increment
	n++;

	// draw first item
	if(bar==0) { printf("["); bar=1; fflush(stdout); }

	// update if needed
	if(ceil((float)maxbar/maxn*n)>bar){
		int d=ceil((float)maxbar/maxn*n)-bar;
		int i;
		for(i=0; i < d; i++)
			if(bar+i+1==maxbar) printf("]\n"); else printf("-"); 
		fflush(stdout);
		bar += d;
	}
}
