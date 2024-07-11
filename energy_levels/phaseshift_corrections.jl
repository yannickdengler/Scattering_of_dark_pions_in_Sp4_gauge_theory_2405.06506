"""
This function calculates the finite volume corrections to the isospin-2
hase shift p⋅cot(δ). It is based on Eq.(31) and Table.I in the paper:
Bedaque, Sato, Walker-Loud hep-lat/0601033
"""
# This is essentially the OEIS sequence A005875 where c(0)=1
function c(n)
    n == 1 && return 6
    n == 2 && return 12
    n == 3 && return 8
    n == 4 && return 6
    n == 5 && return 24
end

"""
In order to visualise the dependence of n in Eq.(31), I compute the correction 
term Δ(p⋅cot(δ)) as: -mπ/√2 ∑ (c(n)⋅A(n))
"""
# truncate the sum in Eq.(31) at the (227/24) term:
function A(n,mπ,L)
    x = n*mπ*L
    return exp(-x)*(1-227/24/x)/sqrt(x)
end
# truncate sum over n at 5
function Δpcotδ(mπ,L)
    Δ = 0.0
    for n in 1:5
        Δ += -mπ*c(n)*A(n,mπ,L)/sqrt(2π)
    end
    return Δ
end

using Pkg; Pkg.activate("./energy_levels/src_jl",io=devnull)
using I2julia
using HDF5

h5io = h5open("./output/hdf5/scattering_b6.900_m-0.870.hdf5")
h5io = h5open("./output/hdf5/scattering_b7.200_m-0.780.hdf5")
h5io = h5open("./output/hdf5/scattering_b6.900_m-0.900.hdf5")  
pcotδ_div_mπ = h5io["orig_P_cot_PS_pipi_prime"][]
pcotδ_div_mπ_sample = h5io["sample_P_cot_PS_pipi_prime"][]
Ls = h5io["orig_N_Ls"][]
mπ = h5io["orig_m_pi_inf"][1]
for i in eachindex(pcotδ_div_mπ)
    Δ = Δpcotδ(mπ,Ls[i]) / mπ
    pcotδ = pcotδ_div_mπ[i] 
    @show Δ, pcotδ
end