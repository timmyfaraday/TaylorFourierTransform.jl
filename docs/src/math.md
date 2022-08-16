# Mathematical Background on Taylor-Fourier Transform

## Nomenclature

| Symbol    | Description                   | Domain    |
|:----------|:------------------------------|:----------|
| `s`       | signal                        | 𝐑         |
| `a`       | amplitude                     | 𝐑⁺        |
| `ϕ`       | phase [(rad)]                 | [-π,π]    |
| `φ`       | anti-rotating phase [(rad)]   | [-π,π]    |
| `ξ`       | dynamic phasor                | 𝐂         |
| `ψ`       | anti-rotating dynamic phasor  | 𝐂         |

| Set       | Description                   | Domain    |
|:----------|:------------------------------|:----------|
| `d ∈ D`   | set of derivatives            | 𝐍         |
| `h ∈ H`   | set of harmonic number        | 𝐍         |

| Parameter | Description                   | Domain    |
|:----------|:------------------------------|:----------|
| `F`       | frequency [Hz]                | 𝐑⁺        |
| `ω`       | angular frequency [(rad)/s]   | 𝐑⁺        |

## Taylor-Fourier Transform

The Taylor-Fourier Transform function `tft()` gives the up-to-Dth derivative of 
the Hth-harmonic dynamic phasors `ξₕ⁽ᵈ⁾(t) ∈ ℂ, ∀ d ∈ {0,..,D}, h ∈ H`. 

### Zeroth Harmonic

The up-to-Dth derivative of the zeroth-harmonic dynamic phasor 
`ξ₀⁽ᵈ⁾(t) ∈ ℝ, ∀ d ∈ {0,..,D}` is given by,
```math
    \xi^{(d)}_{0}(t)        =& a^{(d)}_{0}(t)
```

### Zeroth Derivative

The zeroth derivative of the hth-harmonic dynamic phasor `ξₕ⁽⁰⁾(t)` is given by,
```math
\begin{aligned}
    \xi^{(0)}_{h}(t)        =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \phi^{(0)}_{h}(t)) \\
                            =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \varphi^{(0)}_{h}(t)) \, \exp(j \omega_{h} t) \\
                            =& \psi^{(0)}_{h}(t) \, \exp(j \omega_{h} t)

\end{aligned}
```
where `ψₕ⁽⁰⁾(t)` denotes the anti-rotating zeroth derivative of the hth-harmonic 
dynamic phasor. 

The zeroth derivative of the amplitude `aₕ⁽⁰⁾(t)`, phase `ϕₕ⁽⁰⁾(t)` and anti-rotating 
phase `φₕ⁽⁰⁾(t)` are, respectively,
```math
\begin{aligned}
    a^{(0)}_{h}(t)          =& |2 \, \xi^{(0)}_{h}(t)| \\
    \phi^{(0)}_{h}(t)       =& \angle \big( \xi^{(0)}_{h}(t) \big) \\
    \varphi^{(0)}_{h}(t)    =& \angle \big( \xi^{(0)}_{h}(t) \, \exp(-j \omega_{h} t) \big) = \angle \big( \psi^{(0)}_{h}(t) \big)
\end{aligned}
```
The zeroth derivative hth-harmonic signal `sₕ⁽⁰⁾(t)` and overall signal `s⁽⁰⁾(t)` 
are, respectively,
```math
\begin{aligned}
    s^{(0)}_{h}(t)          =& \xi^{(0)}_{h}(t) + \text{conj}\big( \xi^{(0)}_{h}(t) \big) \\
    s^{(0)}(t)              =& \sum_{h \in H} \xi^{(0)}_{h}(t) + \text{conj}\big( \xi^{(0)}_{h}(t) \big),
\end{aligned}
```
where `conj(ξₕ⁽⁰⁾(t))` denotes the complex conjugate of the zeroth derivative of the 
hth-harmonic dynamic phasor.

### First Derivative

The first derivate of the hth-harmonic dynamic phasor `ξₕ⁽¹⁾(t)` is given by,
```math
\begin{aligned}
    \xi_{h}^{(1)}(t)        =&  \frac{\mathrm{d}\xi_{h}^{(0)}(t)}{\mathrm{d}t} \\
                            =&  \frac{1}{2} \big( a^{(1)}_{h}(t) +
                                j \, \phi^{(1)}_{h}(t) \, a^{(0)}_{h}(t)) \big)
                                \, \exp(j \phi^{(0)}_{h}(t)).

\end{aligned}
```

The first derivative of the hth-harmomic anti-rotating dynamic phasor `ψₕ⁽¹⁾(t)`
is,
```math 
\begin{aligned}
    \psi_{h}^{(1)}(t)       =&  \frac{\mathrm{d}\psi_{h}^{(0)}(t)}{\mathrm{d}t} \\
                            =&  \xi_{h}^{(1)}(t) \, \exp(-j \omega_{h} t).
\end{aligned}
```

The first derivative of the hth-harmonic amplitude `aₕ⁽¹⁾(t)`, phase `ϕₕ⁽¹⁾(t)` and
anti-rotating phase `φₕ⁽¹⁾(t)` are, respectively,
```math
\begin{aligned}
    a_{h}^{(1)}(t)          =& ℜ[2 \, \xi_{h}^{(1)}(t) \, \exp(-j \phi_{h}^{(0)}(t))] \\
    \phi_{h}^{(1)}(t)       =& \frac{ℑ[2 \, \xi_{h}^{(1)}(t) \, \exp(-j \phi_{h}^{(0)}(t))]}{a^{(0)}_{h}(t)} \\
    \varphi_{h}^{(1)}(t)    =& \frac{ℑ[2 \, \psi_{h}^{(1)}(t) \, \exp(-j \varphi_{h}^{(0)}(t))]}{a^{(0)}_{h}(t)}.
\end{aligned}
```
Note that implicitly this means that `ϕₕ⁽¹⁾(t) = φₕ⁽¹⁾(t)`.

### Second Derivative

The second derivate of the hth-harmonic dynamic phasor `ξₕ⁽²⁾(t)` is given by,
```math
\begin{aligned}
    \xi_{h}^{(2)}(t)        =&  \frac{\mathrm{d}\xi_{h}^{(1)}(t)}{\mathrm{d}t} \\
                            =&  \frac{1}{2} \big( 
                                a^{(2)}_{h}(t) - 
                                \phi^{(1)}_{h}(t)^2 \, a^{(0)}_{h}(t) +
                                j \, 2 \, \phi^{(1)}_{h}(t) \, a^{(1)}_{h}(t) +
                                j \, \phi^{(2)}_{h}(t) \, a^{(0)}_{h}(t)) \big)
                                \, \exp(j \phi^{(0)}_{h}(t)).

\end{aligned}
```

The second derivative of the hth-harmomic anti-rotating dynamic phasor `ψₕ⁽²⁾(t)`
is,
```math 
\begin{aligned}
    \psi_{h}^{(2)}(t)       =&  \frac{\mathrm{d}\psi_{h}^{(1)}(t)}{\mathrm{d}t} \\
                            =&  \xi_{h}^{(2)}(t) \exp(-j \omega_{h} t).
\end{aligned}
```

The second derivative of the hth-harmonic amplitude `aₕ⁽²⁾(t)`, phase `ϕₕ⁽²⁾(t)` and
anti-rotating phase `φₕ⁽²⁾(t)` are, respectively,
```math
\begin{aligned}
    a_{h}^{(2)}(t)          =& ℜ[2 \, \xi_{h}^{(2)}(t) \, \exp(-j \phi_{h}^{(0)}(t))] + a_{h}^{(0)}(t) \, \phi_{h}^{(1)}(t)^{2} \\
    \phi_{h}^{(1)}(t)       =& \frac{ℑ[2 \, \xi_{h}^{(2)}(t) \, \exp(-j \phi_{h}^{(0)}(t))] - 2 \, a_{h}^{(1)}(t) \, \phi_{h}^{(1)}(t)}{a^{(0)}_{h}(t)} \\
    \varphi_{h}^{(1)}(t)    =& \frac{ℑ[2 \, \psi_{h}^{(2)}(t) \, \exp(-j \varphi_{h}^{(0)}(t))] - 2 \, a_{h}^{(1)}(t) \, \varphi_{h}^{(1)}(t)}{a^{(0)}_{h}(t)}.
\end{aligned}
```
Note that implicitly this means that `ϕₕ⁽²⁾(t) = φₕ⁽²⁾(t)`.