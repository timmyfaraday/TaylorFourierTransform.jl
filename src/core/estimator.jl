################################################################################
#  Copyright 2022, Tom Van Acker (BASF), Jose Antonio de la O Serna (UANL)     #
################################################################################
# TaylorFourierTransform.jl                                                    #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TaylorFourierTransform.jl                 #
################################################################################

"""
    TaylorFourierTransform.harmonic_estimator(prob::AbstractDTFTProblem, H::Int)

Function to obtain the up-to-Dth-degree derivative of the Hth-harmonic dynamic 
phasor.

Input:
- `prob::AbstractDTFTProblem`   | DTFT problem struct
- `H::Int`                      | harmonic number [-]

Output:
- `X::Matrix{<:Complex}`        | up-to-Dth-degree derivative of the 
                                | Hth-harmonic dynamic phasor
"""
function harmonic_estimator(prob::AbstractDTFTProblem, H::Int)
    # define the extended normalized time range `rU`
    Δu  = 1 / prob.N
    rU  = 0.0:Δu:(prob.K + 1 - Δu)

    # samples of up-to-Dth-degree derivative of the Kth-degree o-spline `Φ`
    Φ   = sample_ospline(prob.D, prob.K, prob.N)
    # update `Φ` to reflect the samples of up-to-Dth-degree derivative of the 
    # Hth-harmonic bandpass filter
    Y   = Φ .* exp.((2 * pi * im * H) .* rU) ./ prob.N .* 
            (prob.F).^collect(0:prob.D)'

    # define the appropriate range in the convolution `rC` using an offset `O`:
    # if (K+1)*N is even    → O = ((K+1)*N) / 2 + 1
    # if (K+1)*N is oneven  → O = ((K+1)*N - 1) / 2 + 1
    # NB: the Julia DSP pkg does not allow for the keyword 'same' in its conv-
    # function, as is the case in MATLAB, the range `rC` mimics that behavior
    O   = floor(Int, ((prob.K + 1) * prob.N) / 2) + 1 
    rC  = O:O+length(prob.s)-1
    
    # return the up-to-Dth-degree derivative of the H-th harmonic dynamic phasor
    return _DSP.conv(prob.s, Y)[rC, :] 
end

"""
    TaylorFourierTransform.fast_estimator(prob::AbstractDTFTProblem)

Function to obtain the Fast Taylor Fourier transform FTFT. Uses the samples of 
the up-to-Dth-degree derivatives of Kth-degree o-spline in matrix Φ.

Input:
- `prob::AbstractDTFTProblem`       | DTFT problem struct

Output:
- `X::Dict{Int,Matrix{<:Complex}}`  | dictionary contain up-to-Dth-degree 
                                    | derivative of the dynamic phasors for the
                                    | full set of harmonics
"""
function fast_estimator(prob::AbstractDTFTProblem)
    # initialize a matrix of the size (|s|,|h|,D) to store the up-to-Dth-degree 
    # derivative of the dynamic phasors for the full set of harmonics
    X   = zeros(Number, length(prob.s), prob.N, prob.D+1)

    # samples of up-to-Dth-degree derivative of the Kth-degree o-spline `Φ`
    Φ   = sample_ospline(prob.D, prob.K, prob.N)

    # determine the shift to center the 'current' signal point
    Δs  = floor(Int, prob.N * (prob.K+1) / 2)

    for (ni,ns) in enumerate(prob.s) if Δs < ni <= length(prob.s) - Δs
        # determine the selected range of the signal with 'center' ni
        rS  = (ni-Δs):(ni+Δs-1) 
        # determine the Hadamar product of the osplines and the selected signal
        hd  = Φ .* prob.s[rS] .* (-prob.F).^collect(0:prob.D)'
        # sum the Hadamar product over o-spline degrees
        shd = [sum(hd[nh:prob.N:end, :], dims=1) for nh in 1:prob.N]
        # perform the fft and store the result at ni
        X[ni,:,:] = _FFTW.fft(reduce(vcat,shd), [1]) / prob.N
    end end

    # return a dictionary of spliced arrays over the harmonic numbers
    return Dict(nh => X[:,nh+1,:] for nh in prob.h)
end