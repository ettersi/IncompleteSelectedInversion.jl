# IncompleteSelectedInversion

Example usage:
```julia
# Load package
using IncompleteSelectedInversion

# Parameters
n = 10
c = n # Cut-off level-of-fill. c = n is complete factorisation

# Create test matrix
A = sprand(n,n,0.1)
A += A' + 5*I
Ap,Ai,Av = unpacksparse(A)

# Factorise and invert
Fp,Fi,Fl = symbolic_ldlt(Ap,Ai,c) 
Fp,Fi,Fl = dropfillin(Fp,Fi,Fl,c) # Drop fillin from precomputed symbolic factorisation
Fv = numeric_ldlt(Ap,Ai,Av,Fp,Fi)
Bv = selinv_ldlt(Fp,Fi,Fv)
B = packsparse(Fp,Fi,Bv)

# Compare results
display(inv(full(A))); println()
display(full(B)); println()
```
