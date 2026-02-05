# Siam Quantum 2

*GitHub repository is for source code, documentation, and reporting issues.* To download executables, please visit the official website: [https://sites.google.com/site/siamquantum](https://sites.google.com/site/siamquantum)

**Siam Quantum 2** is a modular C toolbox supporting a wide range of quantum modeling capabilities, including density functional theory (DFT), second-order Møller-Plesset perturbation theory (MP2), analytic energy gradients, molecular geometry optimization, quantum molecular dynamics (QMD), and minimum energy crossing point (MECP) search.

It is officially published in __The Journal of Chemical Physics__. Please cite:

Teepanis Chachiyo, Hathaithip Chachiyo; Siam Quantum 2: An open-source C toolbox for quantum modeling and electronic structure development. J. Chem. Phys. 7 February 2026; 164 (5): 052501. https://doi.org/10.1063/5.0310183

Siam Quantum 2 is designed for research and education in mind. For research, its speed and accuracy are comparable to the well-established quantum chemistry programs. You could use it to do research and publish your results with confident that your results can be reproduced. For education, it is extensively documented in an easy-to-read math display using Markdown which can also be viewed directly on GitHub.

## 1. Executable

You can download the executables for MacOS, Linux, and Windows from the official website: [https://sites.google.com/site/siamquantum](https://sites.google.com/site/siamquantum)

## 2. Running

Only one executatble file ```sq``` is needed to run Siam Quantum. **SQ is a terminal based program.** You have to save it somewhere and run it from a terminal. For example, running *without any argument* will print all available options. Try this first to check if everything is working properly.

```
./sq
```

The prefix `./` tells the terminal that SQ is in current directory.

In addition, molecular geometries in .XYZ format and a few basis set in GAMESS format are provided in the distribution in the directory `examples/` and `basis/` respectively. You can use these files to perform simple Hartree-Fock calculation such as,


```
./sq  examples/h2o.xyz  basis/321g.txt
```

For advance usage, see the FAQ, programming guide, souce code LaTex note, and compilation guide in the `docs/` directory.

*A tutorial video* to get you started is also available on [YouTube](https://www.youtube.com/watch?v=ZyPOtuprsS4).

## 3. Compiling

LAPACK library is required to build SQ from source code. It should be easy to build LAPACK using a fortran compiler. Optionally, LIBXC library can be called, granting access to other exchange-correlation functionals such as B3LYP, and PBE.

Please look at the `docs/` directory for extensive guides on compiling Siam Quantum on Windows, MacOS, Linux platforms.

