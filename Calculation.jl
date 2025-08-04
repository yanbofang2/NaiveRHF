include("NaiveHF.jl")
using LinearAlgebra
using .NaiveHF
using GaussianBasis

# RHF
# --- Run the following block only if this file is run as a script ---
if abspath(PROGRAM_FILE) == @__FILE__
    
    # factor = 1.88973  # Conversion factor from Angstrom to Bohr
    r = 0.37 * 1.88973  # Bond length of H₂ in Bohr
    str_atoms = """
    H 0.0 0.0 0.37
    H 0.0 0.0 -0.37
    """

    # Choose basis set
    basis_set_name = "sto-3g"  # or another valid basis set name

    # Create BasisSet object
    bset = GaussianBasis.BasisSet(basis_set_name, str_atoms; spherical=true, lib=:libcint)
    ham = AbInitio_RHF(bset)
    # @show size(ham.h)

    # For H₂, there are 2 electrons.
    nelec = 2
    hf = HF(nelec, nothing, nothing, NaN, NaN)

    # Run the SCF procedure
    dm_final, scf_energy, mo_coeff = scf_RHF!(hf, ham)
    println("Final Density Matrix: ", dm_final)

    function Nuclear_Repulsion(charge, coordinates)
        n = size(charge)[1]
        E = 0.0
        for i in 1:n-1
            for j in i+1:n
                r0 = norm(coordinates[i, :] - coordinates[j, :])
                E += charge[i] * charge[j] / r0
            end
        end
        return E
    end
    
    coordinate = [0.0 0.0 r; 0.0 0.0 -r]
    charge = [1.0, 1.0]
    E_unc = Nuclear_Repulsion(charge, coordinate)
    println("Nuclear repulsxion energy = ", E_unc)

    println("Final RHF SCF energy for H₂: ", scf_energy + E_unc)
    
    # Print the orbital energies (in Hartree)
    println("\nOrbital energies (Hartree):")
    for (i, e) in enumerate(hf.mo_energy)
        println("Orbital ", i, ": ", e)
    end
end


# # GHF
# # --- Run the following block only if this file is run as a script ---
# if abspath(PROGRAM_FILE) == @__FILE__
    
#     factor = 1.88973  # Conversion factor from Angstrom to Bohr
#     r = 0.37 * factor  # Bond length of H₂ in Angstrom
#     str_atoms = """
#     H 0.0 0.0 -0.6992001
#     H 0.0 0.0 0.6992001
#     """

#     # Choose basis set
#     basis_set_name = "sto-3g"  # or another valid basis set name

#     # Create BasisSet object
#     bset = GaussianBasis.BasisSet(basis_set_name, str_atoms; spherical=true, lib=:libcint)
#     ham = AbInitio_UGHF(bset)
#     # @show size(ham.h)

#     # For H₂, there are 2 electrons.
#     nelec = 2
#     hf = HF(nelec, nothing, nothing, NaN, NaN)

#     # Run the SCF procedure
#     dm_final, scf_energy, mo_coeff = scf_GHF!(hf, ham)

#     function Nuclear_Repulsion(charge, coordinates)
#         n = size(charge)[1]
#         E = 0.0
#         for i in 1:n-1
#             for j in i+1:n
#                 r0 = norm(coordinates[i, :] - coordinates[j, :])
#                 E += charge[i] * charge[j] / r0
#             end
#         end
#         return E
#     end
    
#     coordinate = [0.0 0.0 -r; 0.0 0.0 r]
#     charge = [1.0, 1.0]
#     E_unc = Nuclear_Repulsion(charge, coordinate)

#     println("Final GHF SCF energy for H₂: ", scf_energy + E_unc)
# end


