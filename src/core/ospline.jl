################################################################################
#  Copyright 2022, Tom Van Acker (BASF), Jose Antonio de la O Serna (UANL)     #
################################################################################
# TaylorFourierTransform.jl                                                    #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TaylorFourierTransform.jl                 #
################################################################################

"""
    TaylorFourierTransform.sample_ospline(D::Int, K::Int, N::Int)

Function to obtain samples of the up-to-Dth-degree derivatives of Kth-degree 
o-spline.

Input:
- `D::Int`  | maximum degree of the derivative [-]
- `K::Int`  | degree of the o-spline [-]
- `N::Int`  | number of samples of the fundamental cycle [-]

Output
- `Φ::Matrix{<:Real}`   | samples of the up-to-Dth-degree derivatives of the 
                        | Kth-degree o-spline 
"""
function sample_ospline(D::Int, K::Int, N::Int)
    # define the normalized time range `rU`
    ΔU  = 1 / N
    rU  = 0.0:ΔU:(1.0-ΔU)
    # define the knots range `rK`
    rK  = [-K:-1..., 1:K...]

    # initialize the Φ-matrix
    Φ = ones((K+1) * N, D+1)

    # fill the Φ-matrices
    for ctr in 0:K
        # one-based counter `cnt` wrt zero-based counter `ctr`
        ctn = ctr + 1
        # range of the first dimension of Φ wrt zero-based counter `ctr`
        rΦ  = (N * ctr + 1):(N * (ctr + 1))
        # base o-spline
        u   = -(K + 1) / 2 + ctr .+ rU
        ψ   = [1, -rK[ctn]] ./ -rK[ctn]
        for nk in 1:K
            Φ[rΦ, 1] .*= (u .- rK[ctr + nk]) ./ -rK[ctr + nk]
            if nk > 1
                ψ   = _DSP.conv(ψ, [1, -rK[ctr + nk]] ./ -rK[ctr + nk])
            end
        end
        # derivative by Horner scheme 
        # NB: the Julia Polynomials pkg, takes its coefficients in reverse, 
        # hence the need for the double reserve when determining ψ.
        for nd in 2:(D+1)
            ψ   = reverse(_POL.derivative(_POL.Polynomial(reverse(ψ))).coeffs)
            Φ[rΦ, nd] .*= ψ[1]
            for nk in 2:(K-nd+2)
                Φ[rΦ, nd] .*= u 
                Φ[rΦ, nd] .+= ψ[nk]
            end
        end
    end

    # return the samples of the o-spline and its derivatives: Φ
    return Φ 
end