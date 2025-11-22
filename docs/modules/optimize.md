# optimize.h and optimize.c

## Introduction

This module performs geometry optimization of a molecule using Newton method, which requires an inverse of an Hessian. The Hessian is updated iteratively because the geometry optimization already has to be done iteratively.

Let $M$ = `nDim` be the dimension of the optimizing parameters, which in this case coincides with 3 times the number of atoms. For example,

$$
E = E(\underbrace{x_1, y_1, z_1, x_2, y_2, z_2, \cdots}_{M \; \text{parameters}})
$$

We seek the set of parameters that minimize the total energy of the molecule. All parameters are combined into a vector as follows:

$$
\vec{R} = \begin{bmatrix}x_1 \\ y_1 \\ z_1 \\ x_2 \\ y_2 \\ z_2 \\ \vdots \end{bmatrix} \quad \longrightarrow \quad E = E(\vec{R})
$$

In initial guess starting geometry has to be given, $\vec{R}_0$. Then, the vector is updated iteratively until a convergence is reached.

In Newton method, the stepping vector which will bring the parameter one step closer to the minimum is defined as:

$$
d\vec{R} = - \widetilde{H}^{-1} \vec{g} \quad \longrightarrow \quad \vec{R}_{n+1} = \vec{R}_n + d\vec{R} 
$$

Here, $\widetilde{H}$ is the Hessian matrix, which essentially is the second derivative of the total energy with respect to the change of molecular coordinate, and $\vec{g}$ is the gradient of the energy.

$$
H_{ij} \equiv \frac{\partial^2 E}{\partial R_i \partial R_j}, \quad g_i \equiv \frac{\partial E}{\partial R_i}
$$

Siam Quantum can compute the gradient of energy $\vec{g}$ _analytically_ for Hartree–Fock and DFT. But the Hessian is generally too computationally intensive to evaluate at every step of the optimization. Therefore, an iterative scheme is used to update the inverse of the Hessian as the optimization progresses toward the minimum.

## subroutine **invHessian_BFGS(...)**

```c
void invHessian_BFGS(
        int nDim,            // matrix dimension
        const double *dR,    // displacement vector
        const double *dGrad, // change of gradient vector,
        double *invHessian){ // original inverse of the Hessian
```

This subroutine computes inverse of Hessian matrix using BFGS scheme as explained in the Equation (C.25b) from Szabo & Ostland, _Modern Quantum Chemistry Introduction to Advanced Electronic Structure Theory_.

Let $\widetilde{G}_{n+1}$ be the inverse of the Hessian for the next iteration, then:

$$
\widetilde{G}_{n+1} = \widetilde{B}\widetilde{G}_n\widetilde{B}^T + \frac{d\vec{R}\,d\vec{R}^T}{\alpha}
$$

where,

$$
\widetilde{B} \equiv \widetilde{I} - \frac{d\vec{R}\,d\vec{g}^T}{\alpha}, \quad \alpha \equiv d\vec{g}^T\,d\vec{R}
$$

Note that according to the rule of matrix multiplication, $d\vec{R}\,d\vec{R}^T$ gives a matrix of $M\times M$ dimension; and $d\vec{g}^T\,d\vec{R}$ gives a real number.

## subroutine **stepVector_Newton(...)**

```c
void stepVector_Newton(
        int nDim,            // dimension of the matrix
        const double *invH,  // pointer to inverse of the Hessian matrix
        const double *grad,  // pointer to gradient vector
        double maxStepSize,  // maximum step size in each direction
        double *dR){         // returned step vector
```

This subroutine simply performs matrix multiplication to compute the stepping vector:

$$
d\vec{R} = - \widetilde{H}^{-1} \vec{g} 
$$

## subroutine **optimize(...)**

```c
void optimize(
        int dbSize,                         // number of record in basis set db
        const struct GTOBasisSet_t *basisDB,// pointer to basis set database
        struct Molecule_t *mol,             // returned molecular coordinate
        struct option_t *opt){              // options
```

This is the main subroutine for this module. Starting with the initial guess geometry, it builds the parameter space, perform single point calculation and request gradient of energy evaluation, leading to the following iterative procedure as explained in the introduction.

To jump start the calculation, the identity matrix is used as the inverse of the Hessian, as typically done.

Let $n$=`nIter` be the index for the current iteration, the algorithm proceeds as follows:

$$
\begin{aligned}
& \text{1. With } \vec{R}_n \quad \text{compute }\vec{g}_n \\
& \text{2. Update the inverse Hessian} \\
& \text{3. Compute the stepping vector }  d\vec{R} \\
& \text{4. Update the parameter space }  \vec{R}_{n+1} = \vec{R}_n + d\vec{R}
\end{aligned}
$$

Before looping back the program checks for the convergence criteria: maximum and root-mean-square of the forces and displacement in the atomic unit. They can be adjusted in `optimize.h`

```c
#define OPT_CONV_FORCEMAX   0.000100
#define OPT_CONV_FORCERMS   0.000080
#define OPT_CONV_DISPMAX    0.000800
#define OPT_CONV_DISPRMS    0.000400
```

In addition, the step-size control is put in place to prevent over-shooting. For example, this code is implemented in the subroutine **stepVector_Newton(...)**

```c
        // truncate the step size to maxStepSize
        r = 0;
        for(i=0; i < nDim; i++) r+= dR[i]*dR[i];
        r = sqrt(r);

        if(r > maxStepSize)
        for(i=0; i < nDim; i++)
                dR[i] = dR[i]*maxStepSize/r;
```

First, it computes the length of the stepping vector and check if it exceeds the threshold. If the length is longer than the threshold, Siam Quantum rescale the stepping vector $d\vec{R}$ to precisely the `MAXSTEPSIZE` value. Currently the threshold is set to 0.5 Bohr, but can be adjusted for developmental purpose in `optimize.h`.

