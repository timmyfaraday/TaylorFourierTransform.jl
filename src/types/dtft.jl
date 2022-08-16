################################################################################
#  Copyright 2022, Tom Van Acker (BASF), Jose Antonio de la O Serna (UANL)     #
################################################################################
# TaylorFourierTransform.jl                                                    #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TaylorFourierTransform.jl                 #
################################################################################

# abstract problem types
"""
    TaylorFourierTransform.AbstractDTFTProblem
    
Base for types which define discrete taylor-fourier transform problems.
"""
abstract type AbstractDTFTProblem end

# problem types
"""
    TaylorFourierTransform.DTFTProblem

Defines a discrete taylor-fourier transform problem.

Fields:
- `s::Vector{<:Number}`         | discrete signal [?]
- `t::StepRangeLen{<:Number}`   | discrete time [s]
- `h::Vector{<:Int}`            | harmonic numbers [-]
- `D::Int`                      | degree of the highest derivative [-]
- `F::Number`                   | fundamental frequency [Hz]
- `K::Int`                      | o-spline degree [-]
- `N::Int`                      | number of samples per fundamental cycle [-]
"""
struct DTFTProblem <: AbstractDTFTProblem
    """"discrete signal"""
    s::Vector{<:Number}
    """discrete time"""
    t::StepRangeLen{<:Number}
    """harmonic numbers"""
    h::Vector{<:Int}
    """degree of the highest derivative"""
    D::Int
    """fundamental frequency"""
    F::Number
    """o-spline degree"""
    K::Int
    """number of samples per fundamental cycle"""
    N::Int

    # base constructor
    DTFTProblem(s::Vector{<:Number}, t::StepRangeLen{<:Number}, h::Vector{<:Int}, 
                D::Int, F::Number, K::Int, N::Int) =
        new(s, t, h, D, F, K, N)
end

# problem builder
function build_problem(s, t, D, F, K)
    # checks
    length(s) == length(t)          || error("the length of the discrete signal and time are inconsistent")
    all(diff(t) .≈ (t[2] - t[1]))   || error("the discrete time is not equidistance")

    # reduce the discrete time to a range `rT`
    rT  = range(first(t), last(t), length=length(t))
    
    # determine the fundamental period [s]
    T   = 1 / F

    # find the number of samples in a fundamental cycle [-]
    N   = findfirst(x -> x ≈ first(t) + T, t) - 1

    # find the harmonics set based on the number of samples in a fundamental cycle
    h   = collect(0:N-2)

    return DTFTProblem(s, rT, h, D, F, K, N)
end
function build_problem(s, t, h, D, F, K)
    # checks
    length(s) == length(t)          || error("the length of the discrete signal and time are inconsistent")
    all(diff(t) .≈ (t[2] - t[1]))   || error("the discrete time is not equidistance")

    # reduce the discrete time to a range `rT`
    rT  = range(first(t), last(t), length=length(t))
    
    # determine the fundamental period [s]
    T   = 1 / F

    # find the number of samples in a fundamental cycle [-]
    N   = findfirst(x -> x ≈ first(t) + T, t) - 1

    return DTFTProblem(s, rT, h, D, F, K, N)
end
    
# solution types
"""
    TaylorFourierTransform.AbstractDTFTSolution

Base for types which define discrete taylor-fourier transform solutions.
"""
abstract type AbstractDTFTSolution end

"""
    TaylorFourierTransform.DTFTSolution

Defines a discrete taylor-fourier transform solution

Fields:
- `X::Dict{<:Int,Matrix{<:Number}}`     | dictionary of dynamic phasors ξₕ⁽ᴰ⁾(t)
- `prob::AbstractDTFTProblem`           | DTFT problem struct
"""
struct DTFTSolution <: AbstractDTFTSolution 
    """dictionary of dynamic phasors ξₕ⁽ᴰ⁾"""
    X::Dict
    """DTFT problem struct"""
    prob::AbstractDTFTProblem
end