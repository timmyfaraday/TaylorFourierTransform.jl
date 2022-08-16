################################################################################
#  Copyright 2022, Tom Van Acker (BASF), Jose Antonio de la O Serna (UANL)     #
################################################################################
# TaylorFourierTransform.jl                                                    #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TaylorFourierTransform.jl                 #
################################################################################

"""
TaylorFourierTransform.ftft(prob::AbstractDTFTProblem)

Input:
- `prob::AbstractDTFTProblem`   | DTFT problem struct

Output:
- `sol::AbstractDTFTSolution`   | DTFT solution struct
"""
ftft(prob::AbstractDTFTProblem) = DTFTSolution(fast_estimator(prob), prob)

"""
TaylorFourierTransform.ftft(s::Vector{<:Number}, t::Vector{<:Number}, D::Int, F::Real, K::Int)

Input:
- `s::Vector{<:Number}` | discrete signal [?]
- `t::Vector{<:Number}` | discrete time [s]
- `D::Int`              | maximum degree of the derivatives [-]
- `F::Real`             | fundamental frequency [Hz]
- `K::Int`              | degree of the o-spline [-]

Output:
- `sol::DTFTSolution`   | DTFT solution struct
"""
ftft(s::Vector{<:Number}, t::Vector{<:Number}, D::Int, F::Number, K::Int) =
    ftft(build_problem(s, t, D, F, K)) 

"""
TaylorFourierTransform.ftft(s::Vector{<:Number}, t::Vector{<:Number}, h::Vector{<:Int}, D::Int, F::Real, K::Int)
    
Input:
- `s::Vector{<:Number}` | discrete signal [?]
- `t::Vector{<:Number}` | discrete time [s]
- `h::Vector{<:Int}`    | harmonic numbers [-]
- `D::Int`              | maximum degree of the derivatives [-]
- `F::Real`             | fundamental frequency [Hz]
- `K::Int`              | degree of the o-spline [-]
    
Output:
- `sol::DTFTSolution`   | DTFT solution struct
"""
ftft(s::Vector{<:Number}, t::Vector{<:Number}, h::Vector{<:Int}, D::Int, F::Number, K::Int) =
    ftft(build_problem(s, t, h, D, F, K)) 