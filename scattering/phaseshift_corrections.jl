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
        Δ += -mπ*c(n)*A(n,mπL,L)/sqrt(2π)
    end
    return Δ
end