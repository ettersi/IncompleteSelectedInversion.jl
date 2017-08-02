# IncompleteSelectedInversion

Example usage:
```julia
# Load package
using IncompleteSelectedInversion

# Parameters
n = 10
c = n    # Cut-off level-of-fill. c = n is complete factorisation
τ = 0.0  # Dropping tolerance. τ = 0.0 is complete factorisation

# Create test matrix
A = sprand(n,n,0.1)
A += A' + 5*I

# Factorise and invert
F  = ldlt(A) 
Fc = cldlt(A,c)
Fτ = τldlt(A,τ)
B  = selinv(F)
Bc = selinv(Fc)
Bτ = selinv(Fτ)

# Check result
@show B == Bc == Bτ
# ^ prints true because we didn't actually drop anything in the incomplete factorisations
display(inv(full(A))); println()
display(full(B)); println()
```
