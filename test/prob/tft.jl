################################################################################
#  Copyright 2022, Tom Van Acker (BASF), Jose Antonio de la O Serna (UANL)     #
################################################################################
# TaylorFourierTransform.jl                                                    #
# A Julia package for Taylor-Fourier Transform.                                #
# See http://github.com/timmyfaraday/TaylorFourierTransform.jl                 #
################################################################################

@testset "TFT" begin

    # fundamental frequency and angular frequency
    F   = 50.0
    ω   = 2 * pi * F

    # tft input
    D   = 2
    K   = 9

    # discrete time
    t   = 0.0:0.001:1.0
    tm  = 0.536
    idm = findfirst(x -> x == tm, t)
    
    @testset "Zeroth Harmonic Constant Signal" begin
        # input - amplitude
        A(t)    = 10.0

        # derivatives - amplitude and anti-rotating phase
        dA(t)   = _FD.derivative(A,t)
        d2A(t)  = _FD.derivative(dA,t)

        # derived input
        Ξ(t)    = A.(t)
        Ψ(t)    = A.(t)
        S(t)    = A.(t)

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [0,10], D, F, K)

        # tests
        ## amplitude
        @test isapprox(TFT.a(sol,0,0)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,0)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,0)[idm], d2A(tm), atol=atol)
        ### --- amplitude of non-present harmonic should be zero
        @test isapprox(TFT.a(sol,0,10)[idm], 0.0, atol=atol)
        ## phase
        @test isapprox(TFT.ϕ(sol,0,0)[idm], 0.0, atol=atol)
        ## anti-rotating phase
        @test isapprox(TFT.φ(sol,0,0)[idm], 0.0, atol=atol)
        ## dynamic phasor
        @test isapprox(TFT.ξ(sol,0,0)[idm], Ξ(tm), atol=atol)
        ## anti-rotating dynamic phasor
        @test isapprox(TFT.ψ(sol,0,0)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,1,0)[idm], TFT.ξ(sol,1,0)[idm], atol=atol)
        @test isapprox(TFT.ψ(sol,2,0)[idm], TFT.ξ(sol,2,0)[idm], atol=atol)
        ## signal        
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,0)[idm], S(tm), atol=atol)
        ### --- signal of non-present harmonic should be zero
        @test isapprox(TFT.signal(sol,10)[idm], 0.0, atol=atol)
        ## error 
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)

    end

    @testset "Zeroth Harmonic Linear Signal" begin
        # input - amplitude
        A(t)    = 10.0 .- t

        # derivatives - amplitude and anti-rotating phase
        dA(t)   = _FD.derivative(A,t)
        d2A(t)  = _FD.derivative(dA,t)

        # derived input
        Ξ(t)    = A.(t)
        Ψ(t)    = A.(t)
        S(t)    = A.(t)

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [0,10], D, F, K)

        # tests
        ## amplitude
        @test isapprox(TFT.a(sol,0,0)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,0)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,0)[idm], d2A(tm), atol=atol)
        ### --- amplitude of non-present harmonic should be zero
        @test isapprox(TFT.a(sol,0,10)[idm], 0.0, atol=atol)
        ## phase
        @test isapprox(TFT.ϕ(sol,0,0)[idm], 0.0, atol=atol)
        ## anti-rotating phase
        @test isapprox(TFT.φ(sol,0,0)[idm], 0.0, atol=atol)
        ## dynamic phasor
        @test isapprox(TFT.ξ(sol,0,0)[idm], Ξ(tm), atol=atol)
        ## anti-rotating dynamic phasor
        @test isapprox(TFT.ψ(sol,0,0)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,1,0)[idm], TFT.ξ(sol,1,0)[idm], atol=atol)
        @test isapprox(TFT.ψ(sol,2,0)[idm], TFT.ξ(sol,2,0)[idm], atol=atol)
        ## signal        
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,0)[idm], S(tm), atol=atol)
        ### --- signal of non-present harmonic should be zero
        @test isapprox(TFT.signal(sol,10)[idm], 0.0, atol=atol)
        ## error 
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)

    end

    @testset "Zeroth Harmonic Quadratic Signal" begin
        # input - amplitude
        A(t)    = 10.0 .- t .+ 2.0 .* t.^2

        # derivatives - amplitude and anti-rotating phase
        dA(t)   = _FD.derivative(A,t)
        d2A(t)  = _FD.derivative(dA,t)

        # derived input
        Ξ(t)    = A.(t)
        Ψ(t)    = A.(t)
        S(t)    = A.(t)

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [0,10], D, F, K)

        # tests
        ## amplitude
        @test isapprox(TFT.a(sol,0,0)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,0)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,0)[idm], d2A(tm), atol=atol)
        ### --- amplitude of non-present harmonic should be zero
        @test isapprox(TFT.a(sol,0,10)[idm], 0.0, atol=atol)
        ## phase
        @test isapprox(TFT.ϕ(sol,0,0)[idm], 0.0, atol=atol)
        ## anti-rotating phase
        @test isapprox(TFT.φ(sol,0,0)[idm], 0.0, atol=atol)
        ## dynamic phasor
        @test isapprox(TFT.ξ(sol,0,0)[idm], Ξ(tm), atol=atol)
        ## anti-rotating dynamic phasor
        @test isapprox(TFT.ψ(sol,0,0)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,1,0)[idm], TFT.ξ(sol,1,0)[idm], atol=atol)
        @test isapprox(TFT.ψ(sol,2,0)[idm], TFT.ξ(sol,2,0)[idm], atol=atol)
        ## signal        
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,0)[idm], S(tm), atol=atol)
        ### --- signal of non-present harmonic should be zero
        @test isapprox(TFT.signal(sol,10)[idm], 0.0, atol=atol)
        ## error 
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)

    end

    @testset "Fundamental Periodic Signal" begin 
        # input - amplitude and anti-rotating phase
        A(t)    = 10.0 
        Φ(t)    = 2 * pi / 3

        # derivatives - amplitude and anti-rotating phase
        dA(t)   = _FD.derivative(A,t)
        d2A(t)  = _FD.derivative(dA,t)
        dΦ(t)   = _FD.derivative(Φ,t)
        d2Φ(t)  = _FD.derivative(dΦ,t)

        # derived input
        Φr(t)   = rem.(ω .*t .+ Φ.(t) .+ pi, 2 * pi) .- pi
        Fr(t)   = F + dΦ(t) / (2 * pi)
        Rr(t)   = d2Φ(t) / (2 * pi)^2
        Ξ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)
        Ψ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t))
        S(t)    = real.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+ 
                  conj.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1,10], D, F, K)

        # tests
        ## amplitude
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,1)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,1)[idm], d2A(tm), atol=atol)
        ### --- amplitude of derivative degree higher than two not supported
        @test TFT.a(sol,3,1) === nothing
        ### --- amplitude of non-present harmonic should be zero
        @test isapprox(TFT.a(sol,0,10)[idm], 0.0, atol=atol)
        ## phase
        @test isapprox(TFT.ϕ(sol,0,1)[idm], Φr(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ### --- phase of derivative degree higher than two not supported
        @test TFT.ϕ(sol,3,1) === nothing
        ## anti-rotating phase
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ### --- ar phase of derivative degree higher than two not supported
        @test TFT.φ(sol,3,1) === nothing
        ## frequency 
        @test isapprox(TFT.f(sol,1)[idm], Fr(tm), atol=atol)
        ## rocof
        @test isapprox(TFT.r(sol,1)[idm], Rr(tm), atol=atol)
        ## dynamic phasor
        @test isapprox(TFT.ξ(sol,0,1)[idm], Ξ(tm), atol=atol)
        ## anti-rotating dynamic phasor
        @test isapprox(TFT.ψ(sol,0,1)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,1,1)[idm], TFT.ξ(sol,1,1)[idm] * exp(-im * ω * tm), atol=atol)
        @test isapprox(TFT.ψ(sol,2,1)[idm], TFT.ξ(sol,2,1)[idm] * exp(-im * ω * tm), atol=atol)
        ## signal        
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1)[idm], S(tm), atol=atol)
        ### --- signal of non-present harmonic should be zero
        @test isapprox(TFT.signal(sol,10)[idm], 0.0, atol=atol)
        ## error 
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

    @testset "Fundamental Aperiodic Linear Signal" begin 
        # input - amplitude and anti-rotating phase
        A(t)    = 10.0 .- t
        Φ(t)    = pi / 2 .* t

        # derivatives - amplitude and anti-rotating phase
        dA(t)   = _FD.derivative(A,t)
        d2A(t)  = _FD.derivative(dA,t)
        dΦ(t)   = _FD.derivative(Φ,t)
        d2Φ(t)  = _FD.derivative(dΦ,t)

        # derived input
        Φr(t)   = rem.(ω .*t .+ Φ.(t) .+ pi, 2 * pi) .- pi
        Fr(t)   = F + dΦ(t) / (2 * pi)
        Rr(t)   = d2Φ(t) / (2 * pi)^2
        Ξ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)
        Ψ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t))
        S(t)    = real.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+ 
                  conj.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1,10], D, F, K)

        # tests
        ## amplitude
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,1)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,1)[idm], d2A(tm), atol=atol)
        ### --- amplitude of derivative degree higher than two not supported
        @test TFT.a(sol,3,1) === nothing
        ### --- amplitude of non-present harmonic should be zero
        @test isapprox(TFT.a(sol,0,10)[idm], 0.0, atol=atol)
        ## phase
        @test isapprox(TFT.ϕ(sol,0,1)[idm], Φr(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ### --- phase of derivative degree higher than two not supported
        @test TFT.ϕ(sol,3,1) === nothing
        ## anti-rotating phase
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ### --- ar phase of derivative degree higher than two not supported
        @test TFT.φ(sol,3,1) === nothing
        ## frequency 
        @test isapprox(TFT.f(sol,1)[idm], Fr(tm), atol=atol)
        ## rocof
        @test isapprox(TFT.r(sol,1)[idm], Rr(tm), atol=atol)
        ## dynamic phasor
        @test isapprox(TFT.ξ(sol,0,1)[idm], Ξ(tm), atol=atol)
        ## anti-rotating dynamic phasor
        @test isapprox(TFT.ψ(sol,0,1)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,1,1)[idm], TFT.ξ(sol,1,1)[idm] * exp(-im * ω * tm), atol=atol)
        @test isapprox(TFT.ψ(sol,2,1)[idm], TFT.ξ(sol,2,1)[idm] * exp(-im * ω * tm), atol=atol)
        ## signal        
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1)[idm], S(tm), atol=atol)
        ### --- signal of non-present harmonic should be zero
        @test isapprox(TFT.signal(sol,10)[idm], 0.0, atol=atol)
        ## error 
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

    @testset "Fundamental Aperiodic Quadratic Signal" begin 
        # input - amplitude and anti-rotating phase
        A(t)    = 10.0 .- t .+ 2.0 .* t.^2
        Φ(t)    = pi / 2 .* t.^2

        # derivatives - amplitude and anti-rotating phase
        dA(t)   = _FD.derivative(A,t)
        d2A(t)  = _FD.derivative(dA,t)
        dΦ(t)   = _FD.derivative(Φ,t)
        d2Φ(t)  = _FD.derivative(dΦ,t)

        # derived input
        Φr(t)   = rem.(ω .* t .+ Φ.(t) .+ pi, 2 * pi) .- pi
        Fr(t)   = F + dΦ(t) / (2 * pi)
        Rr(t)   = d2Φ(t) / (2 * pi)^2
        Ξ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)
        Ψ(t)    = A.(t) ./ 2 .* exp.(im .* Φ.(t))
        S(t)    = real.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t) .+ 
                  conj.(A.(t) ./ 2 .* exp.(im .* Φ.(t)) .* exp.(im .* ω .* t)))

        # perform taylor-fourier transform
        sol = tft(S.(t), collect(t), [1,10], D, F, K)

        # tests
        ## amplitude
        @test isapprox(TFT.a(sol,0,1)[idm], A(tm), atol=atol)
        @test isapprox(TFT.a(sol,1,1)[idm], dA(tm), atol=atol)
        @test isapprox(TFT.a(sol,2,1)[idm], d2A(tm), atol=atol)
        ### --- amplitude of derivative degree higher than two not supported
        @test TFT.a(sol,3,1) === nothing
        ### --- amplitude of non-present harmonic should be zero
        @test isapprox(TFT.a(sol,0,10)[idm], 0.0, atol=atol)
        ## phase
        @test isapprox(TFT.ϕ(sol,0,1)[idm], Φr(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.ϕ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ### --- phase of derivative degree higher than two not supported
        @test TFT.ϕ(sol,3,1) === nothing
        ## anti-rotating phase
        @test isapprox(TFT.φ(sol,0,1)[idm], Φ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,1,1)[idm], dΦ(tm), atol=atol)
        @test isapprox(TFT.φ(sol,2,1)[idm], d2Φ(tm), atol=atol)
        ### --- ar phase of derivative degree higher than two not supported
        @test TFT.φ(sol,3,1) === nothing
        ## frequency 
        @test isapprox(TFT.f(sol,1)[idm], Fr(tm), atol=atol)
        ## rocof
        @test isapprox(TFT.r(sol,1)[idm], Rr(tm), atol=atol)
        ## dynamic phasor
        @test isapprox(TFT.ξ(sol,0,1)[idm], Ξ(tm), atol=atol)
        ## anti-rotating dynamic phasor
        @test isapprox(TFT.ψ(sol,0,1)[idm], Ψ(tm), atol=atol)
        @test isapprox(TFT.ψ(sol,1,1)[idm], TFT.ξ(sol,1,1)[idm] * exp(-im * ω * tm), atol=atol)
        @test isapprox(TFT.ψ(sol,2,1)[idm], TFT.ξ(sol,2,1)[idm] * exp(-im * ω * tm), atol=atol)
        ## signal        
        @test isapprox(TFT.signal(sol)[idm], S(tm), atol=atol)
        @test isapprox(TFT.signal(sol,1)[idm], S(tm), atol=atol)
        ### --- signal of non-present harmonic should be zero
        @test isapprox(TFT.signal(sol,10)[idm], 0.0, atol=atol)
        ## error 
        @test isapprox(TFT.error(sol)[idm], 0.0, atol=atol)
    end

end