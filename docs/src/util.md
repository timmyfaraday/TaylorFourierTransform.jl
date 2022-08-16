# Utilities

A number of functions are made available to the user to retrieve specific
components of the TaylorFourierTransform solution

## Amplitude

```@docs
TaylorFourierTransform.a(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TaylorFourierTransform.amplitude(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

!!! note
    The amplitude `aₕ⁽¹⁾(t)` denotes the first derivate of the amplitude of the 
    zeroth derivative of the dynamic phasor `ξₕ⁽⁰⁾(t)`, not the amplitude of the 
    first derivative of the dynamic phasor `ξₕ⁽¹⁾(t)`.

## Phase

```@docs
TaylorFourierTransform.ϕ(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TaylorFourierTransform.phase(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

!!! note
    The phase `ϕₕ⁽¹⁾(t)` denotes the first derivate of the phase of the 
    zeroth derivative of the dynamic phasor `ξₕ⁽⁰⁾(t)`, not the phase of the first 
    derivative of the dynamic phasor `ξₕ⁽¹⁾(t)`.

## Anti-Rotating Phase
```@docs
TaylorFourierTransform.φ(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TaylorFourierTransform.ar_phase(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

!!! note
    The anti-rotating phase `φₕ⁽¹⁾(t)` denotes the first derivate of the 
    anti-rotating phase of the zeroth derivative of the dynamic phasor 
    `ξₕ⁽⁰⁾(t)`, not the anti-rotating phase of the first derivative of the 
    dynamic phasor `ξₕ⁽¹⁾(t)`.

## Frequency
```@docs
TaylorFourierTransform.f(sol::TaylorFourierTransform.AbstractDTFTSolution, H::Int=1)
```
```@docs
TaylorFourierTransform.frequency(sol::TaylorFourierTransform.AbstractDTFTSolution, H::Int=1)
```

## Rate-Of-Change-Of-Frequency (ROCOF)
```@docs
TaylorFourierTransform.r(sol::TaylorFourierTransform.AbstractDTFTSolution, H::Int=1)
```
```@docs
TaylorFourierTransform.rocof(sol::TaylorFourierTransform.AbstractDTFTSolution, H::Int=1)
```

## Dynamic Phasor
```@docs
TaylorFourierTransform.ξ(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TaylorFourierTransform.phasor(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

## Anti-Rotating Dynamic Phasor
```@docs
TaylorFourierTransform.ψ(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```
```@docs
TaylorFourierTransform.ar_phasor(sol::TaylorFourierTransform.AbstractDTFTSolution, D::Int=0, H::Int=1)
```

## Signal
```@docs
TaylorFourierTransform.signal(sol::TaylorFourierTransform.AbstractDTFTSolution)
```
```@docs
TaylorFourierTransform.signal(sol::TaylorFourierTransform.AbstractDTFTSolution, H::Int)
```

## Error
```@docs
TaylorFourierTransform.error(sol::TaylorFourierTransform.AbstractDTFTSolution)
```