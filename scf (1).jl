using LinearAlgebra
using NLsolve
using GaussianBasis
import TensorOperations: @tensor

export scf_RHF!, scf_GHF!


function buildfock_RHF!(fock::Matrix{Float64}, ham::Hamiltonian, dm::Matrix{Float64})

    v = ham.v
    @tensor fock[p, q] = 2 * v[p, q, r, s] * dm[r, s] - v[p, r, q, s] * dm[r, s]
    # Add core Hamiltonian
    fock .+= ham.h

    return fock  # RHF
end


# Function to compute an initial guess density matrix
function core_hamiltonian_guess_RHF(ham::AbInitio_RHF, nelec::Int)
    # Solve H_core C = S C e (generalized eigenvalue problem)
    mo_energy, mo_coeff = eigen(ham.h, ham.s)
    
    # Construct the initial density matrix from occupied orbitals
    mo_occ = mo_coeff[:, 1:div(nelec, 2)]  # Take occupied orbitals
    return mo_occ * mo_occ'      # Density matrix for closed-shell systems
end


# Later, we may implement `buildfock!(::Matrix{Float64}, ::Hubbard, ::Matrix{Float64})`
# to take advantage of the sparsity of the Hubbard Hamiltonian.
# This requires adding a `U` field to the `Hubbard` composite type.


function scf_RHF!(hf::HF, ham::AbInitio_RHF;
              dm0=nothing, iterations::Int=200, dmconv::Float64=1e-10,
              damping::Float64=0.95, diis_dim::Int=8, diis_start::Int=100)

    if dm0 === nothing
        if hf.mo_coeff === nothing
            dm0 = core_hamiltonian_guess_RHF(ham, hf.nelec)
        else
            mo_occ = @view hf.mo_coeff[:, 1:div(hf.nelec, 2)]
            dm0 = mo_occ * mo_occ'
        end
    end

    nso = size(ham.h, 1)
    fock = Matrix{Float64}(undef, nso, nso)

    function residual!(res::Matrix{Float64}, dm::Matrix{Float64})
        buildfock_RHF!(fock, ham, dm)
        mo_energy, mo_coeff = eigen(Hermitian(fock), Hermitian(ham.s))
        hf.mo_energy = mo_energy
        hf.mo_coeff = mo_coeff
        
        # println("matrix", mo_coeff' * ham.s * mo_coeff)
        # Compute the inverse square root of the overlap matrix S
        # ham.s is the overlap matrix.
        # eigS = eigen(ham.s)
        # S_inv_sqrt = eigS.vectors * Diagonal(1 ./ sqrt.(eigS.values)) * eigS.vectors'
    
        # # Transform the Fock matrix into an orthonormal basis:
        # F_orth = S_inv_sqrt' * fock * S_inv_sqrt
        # mo_energy_orth, mo_coeff_orth = eigen!(Hermitian(F_orth))
        # # Transform the eigenvectors back to the original basis:
        # mo_coeff = S_inv_sqrt * mo_coeff_orth
    
        # # Update HF object with the new orbital energies and coefficients
        # hf.mo_energy = mo_energy_orth  # (the eigenvalues are invariant under the transformation)
        # hf.mo_coeff = mo_coeff

        mo_occ = @view mo_coeff[:, 1:div(hf.nelec, 2)]
       
        # mul!(res, mo_occ, mo_occ')
        # @views energy = 0.5 * 2 * (res[:]' * ham.h[:] + sum(mo_energy[1:div(hf.nelec, 2)]))
        
        Dnew = mo_occ * mo_occ'     # New Density matrix for closed-shell systems

        energy = 0.5 * tr(2 * Dnew * (ham.h + fock))
        hf.delta_energy = energy - hf.energy
        hf.energy = energy
        println("\nE = ", energy, "  dE = ", hf.delta_energy)
        println("Core Hamiltonian = ", ham.h)
        println("Two-electron integrals (ERI) = ", ham.v)
        println("Overlap matrix = ", ham.s)
        println("Fock matrix = ", fock)
        println(fock .- ham.h) 
        display(mo_coeff)

        res .= Dnew .- dm
        return res
    end

    result = nlsolve(residual!, dm0; 
                     method=:anderson, iterations=iterations, ftol=dmconv, xtol=1e-10,
                     beta=1.0-damping, m=diis_dim, aa_start=diis_start,
                     show_trace=true)

    if result.f_converged
        println("\nSCF has converged.")
    else
        println("\nSCF did not converge.")
    end
    
    println("f_converged: ", result.f_converged)
    println("x_converged: ", result.x_converged)
    println("Number of iterations: ", result.iterations)
    println("Final solution: ", result.zero)
    dump(result)

    return result.zero, hf.energy, hf.mo_coeff
end



# function buildfock_GHF!(fock::Matrix{Float64}, ham::Hamiltonian, dm::Matrix{Float64})

#     v = ham.v
#     @tensor fock[p, q] = v[p, q, r, s] * dm[r, s] - v[p, r, q, s] * dm[r, s]
#     # Add core Hamiltonian
#     fock .+= ham.h

#     return fock  # RHF
# end

# # Function to compute an initial guess density matrix
# function core_hamiltonian_guess_GHF(ham::AbInitio_UGHF, nelec::Int)
#     # Solve H_core C = S C e (generalized eigenvalue problem)
#     mo_energy, mo_coeff = eigen(ham.h, ham.s)
    
#     # Construct the initial density matrix from occupied orbitals
#     mo_occ = mo_coeff[:, 1:div(nelec, 2)]  # Take occupied orbitals
#     return 2 * (mo_occ * mo_occ')       # Density matrix for closed-shell systems
# end


# # Later, we may implement `buildfock!(::Matrix{Float64}, ::Hubbard, ::Matrix{Float64})`
# # to take advantage of the sparsity of the Hubbard Hamiltonian.
# # This requires adding a `U` field to the `Hubbard` composite type.


# function scf_GHF!(hf::HF, ham::AbInitio_UGHF;
#               dm0=nothing, iterations::Int=100, dmconv::Float64=1e-8,
#               damping::Float64=0.0, diis_dim::Int=8, diis_start::Int=1)

#     if dm0 === nothing
#         if hf.mo_coeff === nothing
#             dm0 = core_hamiltonian_guess_GHF(ham, hf.nelec)
#         else
#             mo_occ = @view hf.mo_coeff[:, 1:hf.nelec]
#             dm0 = mo_occ * mo_occ'
#         end
#     end

#     nso = size(ham.h, 1)
#     fock = Matrix{Float64}(undef, nso, nso)

#     function residual!(res::Matrix{Float64}, dm::Matrix{Float64})
#         buildfock_GHF!(fock, ham, dm)
#         mo_energy, mo_coeff = eigen!(Hermitian(fock))
#         hf.mo_energy = mo_energy
#         hf.mo_coeff = mo_coeff

#         mo_occ = @view mo_coeff[:, 1:hf.nelec]
#         mul!(res, mo_occ, mo_occ')

#         @views energy = 0.5 * (res[:]' * ham.h[:] + sum(mo_energy[1:hf.nelec]))
#         hf.delta_energy = energy - hf.energy
#         hf.energy = energy
#         println("\nE = ", energy, "  dE = ", hf.delta_energy)

#         res .-= dm
#         return res
#     end

#     result = nlsolve(residual!, dm0;
#                      method=:anderson, iterations=iterations, ftol=dmconv,
#                      beta=1.0-damping, m=diis_dim, aa_start=diis_start,
#                      show_trace=true)

#     if result.f_converged
#         println("\nSCF has converged.")
#     else
#         println("\nSCF did not converge.")
#     end
    
#     return result.zero, hf.energy, hf.mo_coeff
# end



# FC = SCE F' = S-1/2 F S-1/2 F'C' =C' E  c' = s1/2 c c =s-1/2 c'

