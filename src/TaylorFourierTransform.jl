################################################################################
#  Copyright 2022, Tom Van Acker (BASF), Jose Antonio de la O Serna (UANL)     #
################################################################################
# TaylorFourierTransform.jl                                                    #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TaylorFourierTransform.jl                 #
################################################################################

module TaylorFourierTransform

    # import pkg
    import DSP
    import Polynomials
    import FFTW

    # pkg constants 
    const _DSP  = DSP
    const _POL  = Polynomials
    const _FFTW = FFTW 

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    include("types/dtft.jl")

    include("core/estimator.jl")
    include("core/ospline.jl")

    include("prob/tft.jl")
    include("prob/ftft.jl")

    include("util/util.jl")

    # export
    export  tft, ftft

    export  amplitude, phase, ar_phase, frequency, rocof, phasor, ar_phasor, 
            signal, error
    export  a, ϕ, φ, f, r, ξ, ψ

end
