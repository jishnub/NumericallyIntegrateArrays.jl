# NumericallyIntegrateArrays.jl
[![Build Status](https://travis-ci.com/jishnub/NumericallyIntegrateArrays.jl.svg?branch=master)](https://travis-ci.com/jishnub/NumericallyIntegrateArrays.jl)

Julia functions for Simpson's rule and trapezoidal rule analogous to scipy.integrate

# Installation

Install the package using 

```julia
julia> ]
pkg> add https://github.com/jishnub/NumericallyIntegrateArrays.jl.git
```

# Usage 

The package exports the functions `trapz` and `simps` that perform integration using the trapezoidal rule and Simpson's rule respectively.

```julia
julia> x = range(0,π,length=10);

julia> f = sin.(x);

julia> trapz(f,x)
1.9796508112164835

julia> simps(f,x)
1.9995487365804028
```

Multidimensional arrays work as well, where the integration is over the first axis.

```julia
julia> f2 = sin.(x) * cos.(x)';

julia> hcat( simps(f2,x), 2cos.(x) ) # the second column is the expected result
10×2 Array{Float64,2}:
  1.99955    2.0
  1.87896    1.87939
  1.53174    1.53209
  0.999774   1.0
  0.347218   0.347296
 -0.347218  -0.347296
 -0.999774  -1.0
 -1.53174   -1.53209
 -1.87896   -1.87939
 -1.99955   -2.0
```