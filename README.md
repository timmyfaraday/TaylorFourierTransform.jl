# TaylorFourierTransform.jl

<a href="https://github.com/timmyfaraday/TaylorFourierTransform.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/TaylorFourierTransform.jl/workflows/CI/badge.svg"></img></a>
<a href="https://codecov.io/gh/timmyfaraday/TaylorFourierTransform.jl"><img src="https://img.shields.io/codecov/c/github/timmyfaraday/TaylorFourierTransform.jl?logo=Codecov"></img></a>
<a href="https://timmyfaraday.github.io/TaylorFourierTransform.jl/"><img src="https://github.com/timmyfaraday/TaylorFourierTransform.jl/workflows/Documentation/badge.svg"></img></a>

[![Active Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://github.com/timmyfaraday/TaylorFourierTransform.jl)

## Overview

TaylorFourierTransform.jl is a Julia package for Taylor-Fourier transform. Similar to the Fourier transform, a Taylor-Fourier transform enables decomposing a time-variant function into its corresponding temporal frequency components. However, contrary to the Fourier transform, the Taylor-Fourier transform does not require the input signal to be periodic, rather approximates it aperiodicity through a Taylor polynomial. Consequently, the resulting frequency components are called dynamic phasors, rather than, as for the Fourier transform, simply phasors.

## Installation

The latest stable release of TaylorFourierTransform.jl can be installed using the Julia package 
manager:

```julia
(v1.6) pkg> add TaylorFourierTransform
```
This package supports Julia v1.6 and later.

In order to test whether the package works, run:
```julia
(v1.6) pkg> test TaylorFourierTransform
```

## Acknowledgements

The primary developers are:
- Tom Van Acker - BASF ([@timmyfaraday](https://github.com/timmyfaraday)), and 
- Jose Antonio de la O Serna - UANL ([@jadlos](https://github.com/jadlos)).

## License

This code is provided under a BSD license.

