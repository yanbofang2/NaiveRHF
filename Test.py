#!/usr/bin/env python3
from pyscf import gto, scf
import numpy as np

# Define the H2 molecule: two hydrogen atoms with a bond length of 1.3983972316 Bohr
# (this corresponds to about 0.74 Å).
mol = gto.M(
    atom='''
         H 0 0 0.37;
         H 0 0 -0.37
         ''',
    basis='sto-3g',
    unit='Angstrom',  # All integrals are then computed in Angstrom.
    charge=0,
    spin=0          # Closed-shell singlet.
)

print("Atom coordinates (Bohr):", mol.atom_coords())
# Get overlap integrals
S = mol.intor('int1e_ovlp')
print("Overlap matrix S:")
print(S)

# Get kinetic energy integrals
T = mol.intor('int1e_kin')
print("\nKinetic energy matrix T:")
print(T)

# Get nuclear attraction integrals
V = mol.intor('int1e_nuc')
print("\nNuclear attraction matrix V:")
print(V)

# Compute the core Hamiltonian (H = T + V)
H = T + V
print("\nCore Hamiltonian H (T+V):")
print(H)

# Get the two-electron integrals.
# Here, aosym='s8' returns the integrals in a compact (Mulliken) notation.
ERI = mol.intor('int2e', aosym='s8')
print("\nTwo-electron integrals (ERI):")
print(ERI)

# Get the nuclear repulsion energy (in Hartree)
Enuc = mol.energy_nuc()
print("\nNuclear repulsion energy:")
print(Enuc)

# Set up the RHF calculation
mf = scf.RHF(mol)
scf_energy = mf.kernel()
print("\nSCF Energy for H2:", scf_energy)

# Retrieve the final Fock matrix.
# mf.get_fock() computes the Fock matrix using the converged density.
fock_final = mf.get_fock()
eigenergy = np.linalg.eigvals(fock_final)
print(mf.mo_coeff)
print("\nFinal Fock matrix (in Hartree):")
print("eigenenergy of Fock matrix:", eigenergy)
print(fock_final)

# Retrieve the final density matrix.
dm_final = mf.make_rdm1()
print("\nFinal density matrix:")
print(dm_final)

# Optional: Explicit computation of Fock matrix from density, core Hamiltonian, and J/K terms:
hcore = mf.get_hcore()
J, K = mf.get_jk(mol, dm_final)
fock_explicit = hcore + J - 0.5 * K
print("\nExplicitly computed Fock matrix (should match the above):")
print(fock_explicit)

print("\nOrbital energies (Hartree):")
for i, energy in enumerate(mf.mo_energy, start=1):
    print("Orbital {}: {}".format(i, energy))


