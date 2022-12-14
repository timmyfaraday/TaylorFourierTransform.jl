# Mathematical Background on Taylor-Fourier Transform

## Nomenclature

| Symbol    | Description                   | Domain    |
|:----------|:------------------------------|:----------|
| `s`       | signal                        | π         |
| `a`       | amplitude                     | πβΊ        |
| `Ο`       | phase [(rad)]                 | [-Ο,Ο]    |
| `Ο`       | anti-rotating phase [(rad)]   | [-Ο,Ο]    |
| `ΞΎ`       | dynamic phasor                | π         |
| `Ο`       | anti-rotating dynamic phasor  | π         |

| Set       | Description                   | Domain    |
|:----------|:------------------------------|:----------|
| `d β D`   | set of derivatives            | π         |
| `h β H`   | set of harmonic number        | π         |

| Parameter | Description                   | Domain    |
|:----------|:------------------------------|:----------|
| `F`       | frequency [Hz]                | πβΊ        |
| `Ο`       | angular frequency [(rad)/s]   | πβΊ        |

## Taylor-Fourier Transform

The Taylor-Fourier Transform function `tft()` gives the up-to-Dth derivative of 
the Hth-harmonic dynamic phasors `ΞΎββ½α΅βΎ(t) β β, β d β {0,..,D}, h β H`. 

### Zeroth Harmonic

The up-to-Dth derivative of the zeroth-harmonic dynamic phasor 
`ΞΎββ½α΅βΎ(t) β β, β d β {0,..,D}` is given by,
```math
\begin{aligned}
    \xi^{(d)}_{0}(t)        =& a^{(d)}_{0}(t)
\end{aligned}
```

### Zeroth Derivative

The zeroth derivative of the hth-harmonic dynamic phasor `ΞΎββ½β°βΎ(t)` is given by,
```math
\begin{aligned}
    \xi^{(0)}_{h}(t)        =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \phi^{(0)}_{h}(t)) \\
                            =& \frac{a^{(0)}_{h}(t)}{2} \, \exp(j \varphi^{(0)}_{h}(t)) \, \exp(j \omega_{h} t) \\
                            =& \psi^{(0)}_{h}(t) \, \exp(j \omega_{h} t)

\end{aligned}
```
where `Οββ½β°βΎ(t)` denotes the anti-rotating zeroth derivative of the hth-harmonic 
dynamic phasor. 

The zeroth derivative of the amplitude `aββ½β°βΎ(t)`, phase `Οββ½β°βΎ(t)` and anti-rotating 
phase `Οββ½β°βΎ(t)` are, respectively,
```math
\begin{aligned}
    a^{(0)}_{h}(t)          =& |2 \, \xi^{(0)}_{h}(t)| \\
    \phi^{(0)}_{h}(t)       =& \angle \big( \xi^{(0)}_{h}(t) \big) \\
    \varphi^{(0)}_{h}(t)    =& \angle \big( \xi^{(0)}_{h}(t) \, \exp(-j \omega_{h} t) \big) = \angle \big( \psi^{(0)}_{h}(t) \big)
\end{aligned}
```
The zeroth derivative hth-harmonic signal `sββ½β°βΎ(t)` and overall signal `sβ½β°βΎ(t)` 
are, respectively,
```math
\begin{aligned}
    s^{(0)}_{h}(t)          =& \xi^{(0)}_{h}(t) + \text{conj}\big( \xi^{(0)}_{h}(t) \big) \\
    s^{(0)}(t)              =& \sum_{h \in H} \xi^{(0)}_{h}(t) + \text{conj}\big( \xi^{(0)}_{h}(t) \big),
\end{aligned}
```
where `conj(ΞΎββ½β°βΎ(t))` denotes the complex conjugate of the zeroth derivative of the 
hth-harmonic dynamic phasor.

### First Derivative

The first derivate of the hth-harmonic dynamic phasor `ΞΎββ½ΒΉβΎ(t)` is given by,
```math
\begin{aligned}
    \xi_{h}^{(1)}(t)        =&  \frac{\mathrm{d}\xi_{h}^{(0)}(t)}{\mathrm{d}t} \\
                            =&  \frac{1}{2} \big( a^{(1)}_{h}(t) +
                                j \, \phi^{(1)}_{h}(t) \, a^{(0)}_{h}(t)) \big)
                                \, \exp(j \phi^{(0)}_{h}(t)).

\end{aligned}
```

The first derivative of the hth-harmomic anti-rotating dynamic phasor `Οββ½ΒΉβΎ(t)`
is,
```math 
\begin{aligned}
    \psi_{h}^{(1)}(t)       =&  \frac{\mathrm{d}\psi_{h}^{(0)}(t)}{\mathrm{d}t} \\
                            =&  \xi_{h}^{(1)}(t) \, \exp(-j \omega_{h} t).
\end{aligned}
```

The first derivative of the hth-harmonic amplitude `aββ½ΒΉβΎ(t)`, phase `Οββ½ΒΉβΎ(t)` and
anti-rotating phase `Οββ½ΒΉβΎ(t)` are, respectively,
```math
\begin{aligned}
    a_{h}^{(1)}(t)          =& β[2 \, \xi_{h}^{(1)}(t) \, \exp(-j \phi_{h}^{(0)}(t))] \\
    \phi_{h}^{(1)}(t)       =& \frac{β[2 \, \xi_{h}^{(1)}(t) \, \exp(-j \phi_{h}^{(0)}(t))]}{a^{(0)}_{h}(t)} \\
    \varphi_{h}^{(1)}(t)    =& \frac{β[2 \, \psi_{h}^{(1)}(t) \, \exp(-j \varphi_{h}^{(0)}(t))]}{a^{(0)}_{h}(t)}.
\end{aligned}
```
Note that implicitly this means that `Οββ½ΒΉβΎ(t) = Οββ½ΒΉβΎ(t)`.

### Second Derivative

The second derivate of the hth-harmonic dynamic phasor `ΞΎββ½Β²βΎ(t)` is given by,
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

The second derivative of the hth-harmomic anti-rotating dynamic phasor `Οββ½Β²βΎ(t)`
is,
```math 
\begin{aligned}
    \psi_{h}^{(2)}(t)       =&  \frac{\mathrm{d}\psi_{h}^{(1)}(t)}{\mathrm{d}t} \\
                            =&  \xi_{h}^{(2)}(t) \exp(-j \omega_{h} t).
\end{aligned}
```

The second derivative of the hth-harmonic amplitude `aββ½Β²βΎ(t)`, phase `Οββ½Β²βΎ(t)` and
anti-rotating phase `Οββ½Β²βΎ(t)` are, respectively,
```math
\begin{aligned}
    a_{h}^{(2)}(t)          =& β[2 \, \xi_{h}^{(2)}(t) \, \exp(-j \phi_{h}^{(0)}(t))] + a_{h}^{(0)}(t) \, \phi_{h}^{(1)}(t)^{2} \\
    \phi_{h}^{(1)}(t)       =& \frac{β[2 \, \xi_{h}^{(2)}(t) \, \exp(-j \phi_{h}^{(0)}(t))] - 2 \, a_{h}^{(1)}(t) \, \phi_{h}^{(1)}(t)}{a^{(0)}_{h}(t)} \\
    \varphi_{h}^{(1)}(t)    =& \frac{β[2 \, \psi_{h}^{(2)}(t) \, \exp(-j \varphi_{h}^{(0)}(t))] - 2 \, a_{h}^{(1)}(t) \, \varphi_{h}^{(1)}(t)}{a^{(0)}_{h}(t)}.
\end{aligned}
```
Note that implicitly this means that `Οββ½Β²βΎ(t) = Οββ½Β²βΎ(t)`.
