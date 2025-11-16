# uhf.h and uhf.c


## subroutine **uhf_getEtotal**

This is an important subroutine in the module which supports HF and DFT calculations.

```c
static double uhf_getEtotal(
        const int nBasis,
        const struct Molecule_t *mol,
        const struct option_t *opt,
        const double *PA,
        const double *PB,
        const double *H,
        const double *JT,
        const double *KA,
        const double *KB,
        const struct MolGrid_t *grid,
        const double *Exc){  
```
As implemented in Siam Quantum, the subroutine takes several matrices and evaluate the total energy. The expression is

```c
        // core and Coulomb
        for(i=0; i < nBasis; i++)
        for(j=0; j < nBasis; j++){
                E += (H[i*nBasis+j]  + 0.5*JT[i*nBasis+j])
                    *(PA[i*nBasis+j] + PB[i*nBasis+j]);
        }
```

$$E = E_{NN} + \underbrace{ \sum_{ij}P^{(T)}(H_{ij} + \frac{1}{2}J^{(T)}_{ij}) }_{\text{core + Coulomb repulsion}} + E_{xc}$$

$E_{NN}$ is the nuclei repulsion energy. In atomic unit, it is:

$$E_{NN} = \sum_{A>B} \frac{Z_A Z_B}{|\vec{r}_A - \vec{r}_B|}$$

$E_{xc}$ is the exchange-correlation energy, and is computed differently depending on the models. For example, for pure Hartree-Fock exchange,

```c
                for(i=0; i < nBasis; i++)
                for(j=0; j < nBasis; j++){
                        E += (0.5*KA[i*nBasis+j])*PA[i*nBasis+j] +
                             (0.5*KB[i*nBasis+j])*PB[i*nBasis+j];
                }
```

$$E^\text{(HF)}_{xc} = \sum_{ij} \frac{1}{2}K^{(\alpha)}_{ij}P^{(\alpha)}_{ij} + \frac{1}{2}K^{(\beta)}_{ij}P^{(\beta)}_{ij}$$

For pure DFT, where anaytical integrations are not available, the exchange-correlation energy density has to be passed down to the routine. This is precisely the variable ```const double *Exc```.

```c
                for(i=0; i < grid->nPoint; i++)
                        E += Exc[i] * grid->w[i];
```

Here, `w[i]` is the weight of the quadrature numerical integration scheme. Using mathematical notations,

$$E^\text{DFT}_{xc}[\rho] = \int d^3 r \, E_{xc}(\vec{r}) = \sum_{i} w[i] E_{xc}[i]$$

For hybrid-functional, the situation is a little more complicated. We need to know the percent mixture and then add to two above contributions accordingly.

```c
        case METHOD_LIBXC:
                for(i=0; i < nBasis; i++)
                for(j=0; j < nBasis; j++){
                        E += opt->hyb_hf_coef*(0.5*KA[i*nBasis+j])*PA[i*nBasis+j] +
                             opt->hyb_hf_coef*(0.5*KB[i*nBasis+j])*PB[i*nBasis+j];
                }

                for(i=0; i < grid->nPoint; i++)
                        E += Exc[i] * grid->w[i];
        break;
```

Here, the mixture goes by the variable `opt->hyb_hf_coef` which is loaded early when Siam Quantum started, and stored in `opt`, the option data structure.

## subroutine **uhf**

This is the main subroutine for the SCF cycle. It supports both Hartree-Fock and DFT, since both have the same mathematical structure, the only different being the exchange-corelationtion contribution.

```c
double uhf(
        int nBasis,              // number of basis functions
        struct GTOBasis_t * gto, // pointer to function structure
        struct Molecule_t * mol, // pointer to molecule structure
        int nEA,                 // total number of spin up electrons
        int nEB,                 // total number of spin down electrons
        double *CA,              // returned molecular alpha spin orbital
        double *CB,              // returned molecular beta spin orbital 
        double *eA,              // returned eigen values
        double *eB,              // returned eigen values
        struct option_t *opt){   // global option
```

### Overview

Before diving into the details, an overview of the procedure is given here. The subroutine basically solves two generalized eigenvalue problems, one for each spin type.

$$
\begin{align*}
\widetilde{F}^{(\alpha)} \vec{c}^{(\alpha)}_k = \epsilon^{(\alpha)}_k \widetilde{S} \vec{c}^{(\alpha)}_k \\
\widetilde{F}^{(\beta)} \vec{c}^{(\beta)}_k = \epsilon^{(\beta)}_k \widetilde{S} \vec{c}^{(\beta)}_k \\
\end{align*}
$$

The ability to allow unparied spin configuration is sometimes called ___spin-polarized density functional theory___, as supposed to the closed-shell system.

Generally, if the number of basis function is `nBasis`, then there are also `nBasis` of these eigenvalues and eigenvectors.

>**Performance note:** A while back we were trying to speed up the LAPACK library while performing the matrix diagonalization to solve this generalize eigenvalue problem. For a really large basis set, this takes a significantly long time compared to the matrix element evaluation. Finally, we found a simple solution, indicatingt that we have been a fool for a long time. That is, just request only the occupied eigenvectors. For large basis set, `nBasis` can be quite large compared to just the number of occupied orbitals. Hence, we don't need to request the full computation of `nBasis` eigenvectors. Doing so, we nolonger notice any performance drop during matrix diagonialzation within LAPACK library.

The origin of the generalized eigenvalues problem given above come from the eigen equation which is another type of differential equation.

$$
\begin{align*}
\big[-\frac{1}{2}\nabla^2 + \hat{v}_\text{ext} + \hat{v}_J + \hat{v}^{(\alpha)}_{xc} \big] \phi^{(\alpha)}_k(\vec{r}) = \epsilon^{(\alpha)}_k \phi^{(\alpha)}_k(\vec{r}) \\
\big[-\frac{1}{2}\nabla^2 + \hat{v}_\text{ext} + \hat{v}_J + \hat{v}^{(\beta)}_{xc} \big] \phi^{(\beta)}_k(\vec{r}) = \epsilon^{(\beta)}_k \phi^{(\beta)}_k(\vec{r}) 
\end{align*}
$$

Note that only the exchange–correlation potential operators are specific to the spin type. If we expand the molecular orbitals $\phi_k(\vec{r})$ into a linear combination of the basis function, then we will have the generalized eigenvalue problem for each spin type as shown above. The Fock matrix elements are,

$$
\begin{align*}
\widetilde{F}^{(\alpha)}_{ij} = (\chi_i| -\frac{1}{2}\nabla^2 + \hat{v}_\text{ext} + \hat{v}_J + \hat{v}^{(\alpha)}_{xc} | \chi_j) \\
\widetilde{F}^{(\beta)}_{ij} = (\chi_i| -\frac{1}{2}\nabla^2 + \hat{v}_\text{ext} + \hat{v}_J + \hat{v}^{(\beta)}_{xc} | \chi_j)
\end{align*}
$$

Because the matrix elements of both spin types are only partially different, it is a good idea to decompose the Fock matrices into several components, and build the full magtrices as needed.


Understanding the SCF-cycle in this subroutine is perhaps the most important exercise in learning about the electronic structure theory in general.

### Building 1-Electron Matrix Elements

Early in the subroutine, it calls subroutines inside `matrix.c` to compute matrix elements. For example, these lines of code computes the kinetic energy matrix elements.

```c
        /////////////////////////////////////////////
        // Building necessary matrix elements     ///
        /////////////////////////////////////////////

        // report
        printf(
        "Computing 1-electron matrix elements .");fflush(stdout);

        // get kinetic matrix
        for(i=0; i < nBasis; i++)
        for(j=0; j <=i; j++){
                // compute explicitly
                T[i*nBasis+j] = GTO_kinetic(i,j,gto);
                // symmetrize matrix
                T[j*nBasis+i] = T[i*nBasis+j];
        }
```

Be sure to look at a number of defition of subroutines inside `matrix.h` that starts with `GTO_` to see what type of integrals are available.

At the end of this step, we will have the **core Hamitonian** component of the Fock matrices which do not depend on spin type. The core Hamiltonian matrix elements are defined as

$$
H_{ij} = (\chi_i| -\frac{1}{2}\nabla^2 + \hat{v}_\text{ext} | \chi_j) = \int d^3 r \chi_i(\vec{r}) [-\frac{1}{2}\nabla^2 + \hat{v}_\text{ext}]\chi_j(\vec{r})
$$

Sometimes, we call this the 1-electron matrix element are integrated over the positions of only 1 electron.





