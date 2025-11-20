# matrix.h and matrix.c

## Definitions

Here, we lists most often used variables to represent matrix elements. The naming convention for accessing elements in the matrix is `[i*nBasis+j]`, meaning the i<sup>th</sup> row, and j<sup>th</sup> columns.

Also, the chemist's notation for electron integral is particularly useful for writing an integral in a compact form.

$$(a|\hat{O}|b) \equiv \int d^3 r \, a(\vec{r}) b(\vec{r}) \hat{O}  b(\vec{r})$$

Here, $a(\vec{r})$ and $b(\vec{r})$ are any functions, and $\hat{O}$ is an operator acting on $b(\vec{r})$. In some cases, like for electron repulsion energy, a double integral is needed. The chemist's notation reads,

$$(ab|cd) \equiv  \iint d^3 r_1 d^3 r_2 \, a(\vec{r}_1)b(\vec{r}_1)\frac{1}{r_{12}}c(\vec{r}_2)d(\vec{r}_2)$$


### Density Matrix 

```PA[i*nBasis+j] and PB[i*nBasis+j]```  is the density matrix for alpha and beta electron, denoting $P^{(\alpha)}_{ij}, P^{(\beta)}_{ij}$.

A density matrix is the sum over occupied orbitals of the product between the molecular coefficients.

$$P_{ij} = \sum_{k \in \text{occ}} C_{ki} C_{kj}$$

Its direct usage is the computation of the electron density,

$$\rho(\vec{r}) = \sum_{k \in \text{occ}}|\phi_i(\vec{r})|^2 = \sum_{ij}P_{ij}\chi_i(\vec{r})\chi_j(\vec{r})$$

Because the density matrix represent the density of electron, it is also used predominantly to compute electronic energies.

### Core Hamiltonian Matrix

```H[i*nBasis+j]``` the the matrix element of the basis function embracing over the core hamiltonian operator, denoting $H_{ij}$.

$$H_{ij} = (\chi_i|-\frac{1}{2}\nabla^2+\hat{v}_\text{ext}|\chi_j)$$

The **external** potential means potential coming for interaction **other than** the electron themselves. Generally, they are the nuclei attraction potential, but external electric field is also possible.

$$\hat{v}_\text{ext}(\vec{r}) = \sum_{A \in \text{nuclei}} - \frac{Z_A}{|\vec{r}-\vec{R}_A|}$$

Also, the naming convention for denoting the nuclei index is the upper roman letter (i.e., A, B).

### Coulomb Potential Matrix

```JT[i*nBasis+j]```  is the matrix element over Coulomb potential $\hat{v}_J$ due to electron repulsion, denoteing $J^{(T)}_{ij}$

The ```T``` signifies that **total** electron density because the interaction does not distinguish between the alpah and beta elecron, using only the total is sufficient.

$$\hat{v}_J(\vec{r_1}) = \int d^3 r_2 \frac{\rho(\vec{r}_2)}{|\vec{r}_1 - \vec{r}_2|} $$

Therefore, the matrix element is,

$$J^{(T)}_{ij} = (\chi_i|\hat{v}_J|\chi_j)$$

If we expand the density $\rho(\vec{r}_2)$ in terms of the density matrix, we have:

$$J^{(T)}_{ij} = \sum_{pq}P^{(T)}_{pq}(\chi_i\chi_j | \chi_p \chi_q)$$

Note that we are using the chemist's notation to write the integrals.

### Hartree-Fock Exchange Matrix

```KA[i*nBasis+j] and KB[i*nBasis]``` are the matrix element for the HF exchange, denoting $K^{(\alpha)}_{ij},K^{(\beta)}_{ij}$.  

The elements are different between the alpha and beta electrons, so we need the two variants for computations.

$$\begin{align*}
K^{(\alpha)}_{ij} & = -\sum_{pq}P^{(\alpha)}_{pq}(ip|jq) \\
K^{(\beta)}_{ij} & = - \sum_{pq}P^{(\beta)}_{pq}(ip|jq)
\end{align*}
$$

## subroutine **GTO_overlap(...)**

```c
double GTO_overlap(int i,                        // ith basis
                   int j,                        // jth basis
                   const struct GTOBasis_t *gto);// basis database
```

There are many subroutines that compute various matrix elements between a pair of basis functions. The overlap, as illustrated here, will be a starting point to understand its general methods and terminology.

```c
        int iCnt, jCnt;   // contracted function 
        double sum=0.0;   // integral sum

        // looping over contracted functions
        for(iCnt=0; iCnt < gto[i].nContract; iCnt++){
        for(jCnt=0; jCnt < gto[j].nContract; jCnt++){
                sum = sum + gto[i].coef[iCnt] * gto[j].coef[jCnt] *
                            gto[i].norm[iCnt] * gto[j].norm[jCnt] *
                overlap(gto[i].exp[iCnt], gto[i].l,  gto[i].m,  gto[i].n,
                                          gto[i].x0, gto[i].y0, gto[i].z0,
                        gto[j].exp[jCnt], gto[j].l,  gto[j].m,  gto[j].n,
                                          gto[j].x0, gto[j].y0, gto[j].z0);
        }
        }
        return sum;
```

As seen from the loop above, the subroutine itself does not compute the integral. Instead, subroutines in **overlap(...)** in `int.c` will perform the most basic element of the integration—between two primitive gaussian functions.

In `matrix.c`, the subroutine simply handles the _contraction_ between primitive gaussian. This pattern is used throughout the subroutines in this module.

## subroutine **GTO_JK_Matrix_NoSymm(...)**

```c
void GTO_JK_Matrix_NoSymm(
        int nBasis,                    // number of basis functions
        const double *PA,              // density matrix for spin up
        const double *PB,              // density matrix for spin down 
        const struct GTOBasis_t *gto,  // basis set info
        const double *Schwarz,         // pointer to schwarz matrix
        double cutoff,                 // cutoff to ignore
        double *GA,                    // return G for spin up
        double *GB){                   // return G for spin down
```

This is the simplest and the slowest version of the 2-electron integral. It only uses Schwarz screening but no other symmetric properties of the integral. Being the simplest, it is also the best place to start looking at how the code works.

```c
        for(p=0; p < nBasis; p++)
        for(q=0; q < nBasis; q++)
        for(i=0; i < nBasis; i++)
        for(j=0; j < nBasis; j++){

                ...

                // compute two-electron integral
                EE = contr_eri(...);

                GA[p*nBasis+q] = GA[p*nBasis+q] + PT[i*nBasis+j]*EE;
                GA[p*nBasis+i] = GA[p*nBasis+i] - PA[q*nBasis+j]*EE;

                GB[p*nBasis+q] = GB[p*nBasis+q] + PT[i*nBasis+j]*EE;
                GB[p*nBasis+i] = GB[p*nBasis+i] - PB[q*nBasis+j]*EE;

        }
```

The code above shows the main loop which goes over all 4 indices $(pq|ij)$. If no screening were used, the entire computation would have taken $N^4$ step to complete. Not to mention each `EE` is also very expensive. Notice how the code is calling **contr_eri(...)**, a subroutine also inside `matrix.c`.

## subroutine **contr_eri(...)**

```c
double contr_eri(
             int lena,double *aexps,double *acoefs,double *anorms,
             double xa,double ya,double za,int la,int ma,int na,
             int lenb,double *bexps,double *bcoefs,double *bnorms,
             double xb,double yb,double zb,int lb,int mb,int nb,
             int lenc,double *cexps,double *ccoefs,double *cnorms,
             double xc,double yc,double zc,int lc,int mc,int nc,
             int lend,double *dexps,double *dcoefs,double *dnorms,
             double xd,double yd,double zd,int ld,int md,int nd){
```

This subroutine works similar to the **GTO_overlap(...)**, except that the arguments are rather straightforward arrays for coefficients and weights. We shall look the the main loop as follows:

```c
        for (i=0; i<lena; i++)
        for (j=0; j<lenb; j++)
        for (k=0; k<lenc; k++)
        for (l=0; l<lend; l++){
                // compute element
                EE  = acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]
                      * eri(xa,ya,za,anorms[i],la,ma,na,aexps[i],
                            xb,yb,zb,bnorms[j],lb,mb,nb,bexps[j],
                            xc,yc,zc,cnorms[k],lc,mc,nc,cexps[k],
                            xd,yd,zd,dnorms[l],ld,md,nd,dexps[l]);
                val += EE;
        }
```

Notice how there are 4 indices: `(i,j,k,l)`, making a very computationally intensive task. Again, the subroutine **eri(...)** in `int.c` actually the one that does all the mathematically complicated recursion using the Taketa _et al._ method, as explained in the Siam Quantum paper (legacy version) in details.