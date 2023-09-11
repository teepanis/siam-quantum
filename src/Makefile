SRC = sq.c fgamma.c int.c util.c basis.c mol.c lin.c matrix.c rhf.c uhf.c option.c xsf.c conv.c grad.c optimize.c pop.c check.c mp2.c multipole.c rpc.c mecp.c qmd.c grid.c dft.c hessian.c
OBJ = sq.o fgamma.o int.o util.o basis.o mol.o lin.o matrix.o rhf.o uhf.o option.o xsf.o conv.o grad.o optimize.o pop.o check.o mp2.o multipole.o rpc.o mecp.o qmd.o grid.o dft.o hessian.o

LIB_PATH = ./libs_macos
INC_PATH = ./libs_macos 

LIB = ${LIB_PATH}/liblapack.a ${LIB_PATH}/librefblas.a  ${LIB_PATH}/libxc.a

#CC  = icc
#CCOPT = -O3 -funroll-loops -fno-alias -m64 -Wall
#LDOPT = -O3 -m64 -Wall

## GNU Linux using gcc compiler
#CC  = gcc
#CCOPT = -O3 -m64 -Wall -msse2 -DLIBXC -I${INC_PATH}
#LDOPT = -O3 -m64 -Wall

## macOS using homebrew compiler
CC = /opt/homebrew/bin/gcc-13
CCOPT = -O3 -m64 -Wall -DLIBXC -I${INC_PATH}
LDOPT = -O3 -m64 -Wall

all: sq

%.o: %.c
	${CC} ${CCOPT} -c $<

sq: ${OBJ}
	${CC} ${LDOPT} -o sq $^ -lm -lgfortran ${LIB}
#	gfortran -o sq $^ -lm -static ${LIB}

clean:
	rm -f sq
	rm -f *.o

#libsq.so: ${SRC}
#	${CC} -Wall -shared -O3 -o libsq.so ${SRC}

#all: libsq.a sq

#libsq.a:${SRC}
#	${CC} ${CCOPT} -c *.c
#	ar cr libsq.a ${OBJ}
#	rm -f *.o	

#sq: sq.c libsq.a
#	${CC} ${CCOPT} -c sq.c
#	${CC} ${LDOPT} -o sq sq.o libsq.a -lm -lgfortran ${LIB} 
#	rm -f sq.o

#clean:
#	rm -f sq
#	rm -f libsq.a
