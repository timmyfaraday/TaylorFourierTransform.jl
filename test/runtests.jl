################################################################################
#  Copyright 2022, Tom Van Acker (BASF), Jose Antonio de la O Serna (UANL)     #
################################################################################
# TaylorFourierTransform.jl                                                    #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TaylorFourierTransform.jl                 #
################################################################################

# using pkgs
using ForwardDiff
using TaylorFourierTransform
using Test

# pkg constants
const _FD = ForwardDiff
const TFT = TaylorFourierTransform

# tolerances
atol = 1e-6

@testset "TaylorFourierTransform.jl" begin

    include("prob/tft.jl")
    include("prob/ftft.jl")

end
