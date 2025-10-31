# basis.c and basis.h

## Definition

Gaussian basis functions are one of the key elements in quantum chemistry. They are the basic elements in which the molecular orbitals are built upon.

In the file **basis.h**, the <em>contracted</em> basis function data structure are defined as follows.

```c
struct GTOBasis_t{
        int           nContract; // number of contracted functions
        int           l;         // x-coordinate angular index
        int           m;         // y-coordinate angular index
        int           n;         // z-coordinate angular index
        double        x0;        // x-coordinate center
        double        y0;        // y-coordinate center
        double        z0;        // z-coordinate center
        double        *coef;     // array for coeffients
        double        *exp;      // array for exponent
        double        *norm;     // array for normalization factor
};
```

It represents the following mathematical form.

$$
\chi(\vec{r}; \vec{r}_A) = (x-x_A)^\ell (y-y_A)^m (z-z_A)^n \sum_{k} N_k w_k e^{-\alpha_n |\vec{r}-\vec{r}_A|^2}
$$

Here, the angular part $(x-x_0)^\ell (y-y_0)^m (z-z_0)^n$ controls the shape of the Gaussian; while the radial part is said to be <em>contracted</em>, meaning it is a linear combination of sum of multiple radial distributions. The vector $\vec{r}_A$ controls the center of the Gaussian basis.

From the C implementation above, **nContract** is the number of contracted Gaussian function. **coef** is the array of $w_k$, **exp** the array of $\alpha_k$, and **norm* the array of $N_k$, respectively.

Note that the normalization factor $N_k$ might seem redundant at first because the product $N_k w_k$ is always present in the mathematical definition. However, in various basis functions (e.g. 3-21G or STO-3G), only $w_k$ are defined, while $N_k$ are <em>implied</em>. Hence, having both of them predetermined and stored in the data structure can be useful.

The normalization factor $N_k$ is such that

$$
\int d^3r \, [x^\ell y^m z^n N_k e^{-\alpha_k r^2}]^2 = 1
$$

## subroutine **genBasis**

```c
struct GTOBasis_t * genBasis(
        struct Molecule_t *mol,               // pointer to molecular structure
        int *nBasis,                          // return the number of basis created
        int dbSize,                           // number of record in basis set database
        const struct GTOBasisSet_t *basisDB){ // pointer to basis set database
```

This is the main subroutine for generating a set of basis functions $\{\chi_i\}$ for a given molecular structure (e.g. XYZ file). The subroutine returns **nBasis** which represents the number of basis functions in the molecular orbitals. ssdsdsads 

**basisDB** point to the basis set database (e.g. 3-21G, STO-3G). Note that the <em>basis set</em> is a set of angular $\{\ell,m,n\}$ and its contracted radial distribution $\{w_k\}$ with $\{N_k\}$ implied by the definition above, specifically designed for each type of atom. It does not contain the position of the atoms.

**mol** specifies the positions of the atoms in the molecule. 

When combining **basisDB** and **mol**, the subroutine is able to create the complete set of the basis functions in which the molecular orbitals can be expanded as follows.

$$
\phi_{k}(\vec{r}) = \sum_i {C}_{ki} \chi_i (\vec{r})
$$

Here, $k$ is the molecular orbital index, and $i$ loops through all the basis functions. $C_{ki}$ is the molecular orbital coefficient matrix. It can be said that solving the electronic structure problem is essentially finding coefficient $C_{ki}$ for which the set of molecular orbital $\{\phi_k\}$ can be properly defined, and the Slater determinant wave-function $\tilde{\Psi}(\vec{r}_1, \vec{r}_2, \cdots)$ subsequently constructed.   
