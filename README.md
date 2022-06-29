# wasm-math

A simple math library, written in Rust, targeting WebAssembly.

_WARNING: STILL DOESN'T COMPILE CORRECTLY TO WASM!_

## Implemented functions
Some functions are implemented from scratch, while others are just
wrappers for Rust's library functions. Here's a complete list:

### Algebra
- square root
- power
- exponential
- logarithm

### Trigonometry
- sine
- cosine
- tangent
- arcsine
- arccosine
- arctangent

### Statistics
- normal distribution
  - probability density distribution
    - standard
    - non-standard
  - cumulate density distribution
  - quantile
- gamma functions
  - complete 
  - probability density distribution
  - cumulate density distribution
  - inverse
  - quantile
- beta
  - beta function
  - incomplete beta function
  - regularized beta function
  - probability density distribution
  - cumulate density distribution
  - quantile
- Student's t distribution
  - probability density distribution
  - cumulate density distribution
  - quantile
- Chi^2 distribution
  - probability density distribution
  - cumulate density distribution
  - quantile

### Calculus
_Note: both functions are not accessible from WASM._
- numerical integrator
- numerical differentiator

### Optimizers
_Note: all of these functions are not accessible from WASM._
- numeric maximum finder
- numeric minimum finder
- numeric stationary points finder
- numeric root finder

### Matrix Math
- square and regular matrices
- matrix operations
  - sum
  - multiplication
  - inversion
  - transposition
