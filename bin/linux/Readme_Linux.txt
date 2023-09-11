The following information will help you get Siam Quantum running on Linux.

1) Simple commands: run a Hartree-Fock, and geometry optimization

# ./sq examples/h2o.xyz basis/631g_star.txt
# ./sq examples/h2o.xyz basis/631g_star.txt -OPT


2) DFT calculations with B3LYP functionals and CHACHIYO GGA

# ./sq examples/h2o.xyz basis/631g_star.txt -DFT=B3LYP
# ./sq examples/h2o.xyz basis/631g_star.txt -DFT=CHACHIYO 


3) Get all avialable options

# ./sq


Please visit our website https://sites.google.com/site/siamquantum/
and look at the FAQ under Documentation section for details
on how to use the program.


Happy Computing,

Teepanis Chachiyo <teepanisc@nu.ac.th>
Hathaithip Chachiyo, <hathaithip.chachiyo@gmail.com>

Oct, 2019.
