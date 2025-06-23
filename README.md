This repository contains data and code associated with the results presented in the paper:

## "Cooper-pair splitters as circuit elements for realizing topological superconductors”
by Guilherme Delfino, Dmitry Green, Saulius Vaitiekėnas, Charles M. Marcus, Claudio Chamon
arXiv:2408.06420

### Oberview
The content of this repository is organized into three main folders, each corresponding to a part of the numerical analysis carried out in the paper:

**minimal_model/** *(julia)*


Numerical study of the effective (coarse-grained) model for the array of Cooper-pair splitters.
Analysis includes self-consistent imposition of superconducting phase constraints in each sublattice of the Kagome geometry.
Focus computing Chern number and constructing phase diagram as a function of flux and superconducting pairing.

**molecule_model/** *(julia)*

Contains exact diagonalization codes for the extended "molecule" Hamiltonian.
Includes study of band structure and computation of Chern numbers using momentum-space techniques.

**Kwant/** *(jupiter notebook)*
Transport simulations of a single splitter using the Kwant Python package.
Includes scattering matrix calculations and conductance plots, providing complementary insights into the tunneling features of the splitter device.
