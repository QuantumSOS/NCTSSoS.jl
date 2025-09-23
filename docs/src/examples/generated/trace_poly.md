```@meta
EditURL = "../literate/trace_poly.jl"
```

# Tracial Polynomial Optimization

## Toy Example
Let's learn how to do [tracial polynomial optimization](@ref
tracial-polynomial) from a toy example.

We use [`NCTSSoS.FastPolynomials.tr`](@ref) to declare a part of a term in
tracial polynomial.

````julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial
@ncpolyvar x[1:3]

p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(Monomial)
````

````
1.0 * tr(x₃¹) * tr(x₁¹x₂¹) * 1 + 1.0 * tr(x₁¹x₂¹x₃¹) * 1
````

Polynomial Optimization declaration and solving interface is the same as regular
polynomial optimization.

````julia
spop = polyopt(p; is_projective=true, comm_gps=[x])

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)

result = cs_nctssos(spop, solver_config)

@assert isapprox(result.objective , -0.046717378455438933, atol = 1e-6)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=3)

result = cs_nctssos(spop, solver_config)

@assert isapprox(result.objective, -0.03124998978001017, atol = 1e-6)
````

````
Problem
  Name                   :                 
  Objective sense        : maximize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 81              
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 1               
  Matrix variables       : 1 (scalarized: 496)
  Integer variables      : 0               

Optimizer started.
Presolve started.
Linear dependency checker started.
Linear dependency checker terminated.
Eliminator started.
Freed constraints in eliminator : 0
Eliminator terminated.
Eliminator - tries                  : 1                 time                   : 0.00            
Lin. dep.  - tries                  : 1                 time                   : 0.00            
Lin. dep.  - primal attempts        : 1                 successes              : 1               
Lin. dep.  - dual attempts          : 0                 successes              : 0               
Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
Presolve terminated. Time: 0.00    
Optimizer  - threads                : 12              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 81              
Optimizer  - Cones                  : 1               
Optimizer  - Scalar variables       : 2                 conic                  : 2               
Optimizer  - Semi-definite variables: 1                 scalarized             : 496             
Factor     - setup time             : 0.00            
Factor     - dense det. time        : 0.00              GP order time          : 0.00            
Factor     - nonzeros before factor : 3321              after factor           : 3321            
Factor     - dense dim.             : 0                 flops                  : 5.82e+05        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  0.00  
1   5.1e-01  2.5e-01  9.0e-02  8.85e-01   -1.203138438e-01  -7.729174112e-03  2.5e-01  0.00  
2   1.7e-01  8.6e-02  1.5e-02  2.46e+00   -4.819246707e-02  -3.625632780e-02  8.6e-02  0.00  
3   3.9e-02  1.9e-02  1.5e-03  1.31e+00   -5.002577131e-02  -4.747173568e-02  1.9e-02  0.00  
4   1.5e-02  7.4e-03  3.3e-04  1.04e+00   -4.859732575e-02  -4.744618389e-02  7.4e-03  0.00  
5   4.8e-03  2.4e-03  5.9e-05  1.03e+00   -4.767388974e-02  -4.725934756e-02  2.4e-03  0.00  
6   1.1e-03  5.3e-04  5.3e-06  1.11e+00   -4.657974322e-02  -4.647086106e-02  5.3e-04  0.00  
7   2.0e-04  9.9e-05  4.0e-07  1.08e+00   -4.670949611e-02  -4.668850464e-02  9.9e-05  0.00  
8   3.7e-05  1.9e-05  3.2e-08  1.03e+00   -4.671561937e-02  -4.671164060e-02  1.9e-05  0.00  
9   1.6e-06  8.1e-07  2.9e-10  1.01e+00   -4.671741226e-02  -4.671723469e-02  8.1e-07  0.01  
10  4.1e-08  2.0e-08  1.1e-12  1.00e+00   -4.671737846e-02  -4.671737400e-02  2.0e-08  0.01  
Optimizer terminated. Time: 0.01    

Problem
  Name                   :                 
  Objective sense        : maximize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 395             
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 1               
  Matrix variables       : 1 (scalarized: 5886)
  Integer variables      : 0               

Optimizer started.
Presolve started.
Linear dependency checker started.
Linear dependency checker terminated.
Eliminator started.
Freed constraints in eliminator : 0
Eliminator terminated.
Eliminator - tries                  : 1                 time                   : 0.00            
Lin. dep.  - tries                  : 1                 time                   : 0.00            
Lin. dep.  - primal attempts        : 1                 successes              : 1               
Lin. dep.  - dual attempts          : 0                 successes              : 0               
Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
Presolve terminated. Time: 0.00    
Optimizer  - threads                : 12              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 395             
Optimizer  - Cones                  : 1               
Optimizer  - Scalar variables       : 2                 conic                  : 2               
Optimizer  - Semi-definite variables: 1                 scalarized             : 5886            
Factor     - setup time             : 0.00            
Factor     - dense det. time        : 0.00              GP order time          : 0.00            
Factor     - nonzeros before factor : 7.82e+04          after factor           : 7.82e+04        
Factor     - dense dim.             : 0                 flops                  : 6.45e+07        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   2.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  0.00  
1   6.1e-01  3.0e-01  1.2e-01  9.56e-01   -1.098826759e-01  8.384488153e-03   3.0e-01  0.01  
2   2.4e-01  1.2e-01  2.0e-02  2.74e+00   -1.722756033e-02  2.650243250e-03   1.2e-01  0.02  
3   7.4e-02  3.7e-02  2.7e-03  2.08e+00   -1.942829823e-02  -1.541146059e-02  3.7e-02  0.02  
4   2.8e-02  1.4e-02  5.5e-04  1.48e+00   -1.904050130e-02  -1.767364397e-02  1.4e-02  0.03  
5   1.1e-02  5.5e-03  1.4e-04  1.04e+00   -2.539912202e-02  -2.482691544e-02  5.5e-03  0.04  
6   3.2e-03  1.6e-03  2.2e-05  1.03e+00   -2.905843863e-02  -2.887813024e-02  1.6e-03  0.04  
7   1.1e-03  5.4e-04  4.3e-06  9.81e-01   -3.055329434e-02  -3.049111225e-02  5.4e-04  0.05  
8   3.0e-05  1.5e-05  1.8e-08  1.01e+00   -3.123181602e-02  -3.122977473e-02  1.5e-05  0.06  
9   8.7e-07  4.4e-07  8.8e-11  1.00e+00   -3.124951147e-02  -3.124945186e-02  4.4e-07  0.06  
10  2.0e-08  1.0e-08  3.0e-13  1.00e+00   -3.124998978e-02  -3.124998842e-02  1.0e-08  0.07  
Optimizer terminated. Time: 0.07    


````

The results matches within $10^{-6}$ absolute tolerance comparing to answer in
[klep2022Optimization](@cite)!

## Polynomial Bell Inequalities

Polynomial Bell inequalities provide a powerful framework for detecting quantum
entanglement and non-locality in bipartite quantum systems. These inequalities
impose constraints on the correlations that can be achieved by local hidden
variable models, and their violation serves as a signature of quantum mechanical
behavior. For maximally entangled bipartite states, such as Bell states, the
quantum correlations can exceed the classical bounds imposed by these polynomial
inequalities, demonstrating the non-local nature of quantum entanglement. The
following examples illustrate how tracial polynomial optimization can be used to
compute the maximum violation of specific Bell inequalities, revealing the
extent to which quantum mechanics transcends classical limitations.

````julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial

@ncpolyvar x[1:2] y[1:2]

p = -1.0 * tr(x[1] * y[1]) - 1.0 * tr(x[1] * y[2]) - 1.0 * tr(x[2] * y[1]) + 1.0 * tr(x[2] * y[2])

tpop = polyopt(p * one(Monomial); is_unipotent=true)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=1, ts_algo=MaximalElimination())

result = cs_nctssos(tpop, solver_config)

@assert isapprox(result.objective, -2.8284271157283083, atol = 1e-5)
````

````
Problem
  Name                   :                 
  Objective sense        : maximize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 17              
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 1               
  Matrix variables       : 2 (scalarized: 37)
  Integer variables      : 0               

Optimizer started.
Presolve started.
Linear dependency checker started.
Linear dependency checker terminated.
Eliminator started.
Freed constraints in eliminator : 0
Eliminator terminated.
Eliminator - tries                  : 1                 time                   : 0.00            
Lin. dep.  - tries                  : 1                 time                   : 0.00            
Lin. dep.  - primal attempts        : 1                 successes              : 1               
Lin. dep.  - dual attempts          : 0                 successes              : 0               
Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
Presolve terminated. Time: 0.00    
Optimizer  - threads                : 12              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 17              
Optimizer  - Cones                  : 1               
Optimizer  - Scalar variables       : 3                 conic                  : 2               
Optimizer  - Semi-definite variables: 1                 scalarized             : 36              
Factor     - setup time             : 0.00            
Factor     - dense det. time        : 0.00              GP order time          : 0.00            
Factor     - nonzeros before factor : 153               after factor           : 153             
Factor     - dense dim.             : 0                 flops                  : 4.70e+03        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   4.0e+00  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  0.00  
1   1.6e+00  4.0e-01  4.2e-01  -4.35e-02  -1.045860826e+00  -1.478035613e+00  4.0e-01  0.00  
2   4.0e-01  1.0e-01  2.2e-02  8.68e-01   -2.641906266e+00  -2.542499069e+00  1.0e-01  0.00  
3   1.1e-02  2.8e-03  7.4e-05  9.90e-01   -2.828991473e+00  -2.825523398e+00  2.8e-03  0.00  
4   3.0e-07  7.6e-08  7.2e-12  1.01e+00   -2.828427258e+00  -2.828427157e+00  7.6e-08  0.00  
5   3.8e-14  1.1e-14  1.0e-18  1.00e+00   -2.828427125e+00  -2.828427125e+00  1.1e-14  0.00  
Optimizer terminated. Time: 0.00    


````

Our computation matches with the theoretical prediction for maximally entangled
bipartite state with $10^{-6}$ absolute tolerance [klep2022Optimization](@cite)!

## Covariance of quantum correlation

As introduced in [Bell Inequalities example](@ref bell-inequalities), we may
also compute the covariance of quantum correlations while limiting the state to
maximally entangled bipartite state.

````julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial

@ncpolyvar x[1:3] y[1:3]

cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))
tpop = polyopt(p * one(Monomial); is_unipotent=true)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)

result = cs_nctssos(tpop, solver_config)

@assert isapprox(result.objective,-5.0, atol = 1e-5)
````

````
Problem
  Name                   :                 
  Objective sense        : maximize        
  Type                   : CONIC (conic optimization problem)
  Constraints            : 1010            
  Affine conic cons.     : 0               
  Disjunctive cons.      : 0               
  Cones                  : 0               
  Scalar variables       : 1               
  Matrix variables       : 1 (scalarized: 6670)
  Integer variables      : 0               

Optimizer started.
Presolve started.
Linear dependency checker started.
Linear dependency checker terminated.
Eliminator started.
Freed constraints in eliminator : 0
Eliminator terminated.
Eliminator - tries                  : 1                 time                   : 0.00            
Lin. dep.  - tries                  : 1                 time                   : 0.00            
Lin. dep.  - primal attempts        : 1                 successes              : 1               
Lin. dep.  - dual attempts          : 0                 successes              : 0               
Lin. dep.  - primal deps.           : 0                 dual deps.             : 0               
Presolve terminated. Time: 0.00    
Optimizer  - threads                : 12              
Optimizer  - solved problem         : the primal      
Optimizer  - Constraints            : 1010            
Optimizer  - Cones                  : 1               
Optimizer  - Scalar variables       : 2                 conic                  : 2               
Optimizer  - Semi-definite variables: 1                 scalarized             : 6670            
Factor     - setup time             : 0.01            
Factor     - dense det. time        : 0.00              GP order time          : 0.00            
Factor     - nonzeros before factor : 5.11e+05          after factor           : 5.11e+05        
Factor     - dense dim.             : 0                 flops                  : 4.38e+08        
ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME  
0   3.7e+01  1.0e+00  1.0e+00  0.00e+00   -0.000000000e+00  -0.000000000e+00  1.0e+00  0.01  
1   1.0e+01  2.7e-01  1.1e-01  6.23e-01   -7.830053809e-01  -6.637260067e-01  2.7e-01  0.03  
2   3.3e+00  9.0e-02  1.9e-02  1.90e+00   -2.817322144e+00  -2.794549298e+00  9.0e-02  0.05  
3   7.5e-01  2.0e-02  2.0e-03  1.32e+00   -4.570370369e+00  -4.567256467e+00  2.0e-02  0.07  
4   1.8e-02  4.8e-04  5.5e-06  1.11e+00   -4.990174678e+00  -4.990017169e+00  4.8e-04  0.09  
5   2.6e-04  7.0e-06  9.4e-09  1.00e+00   -4.999850733e+00  -4.999848272e+00  7.0e-06  0.11  
6   4.4e-07  1.2e-08  6.4e-13  1.00e+00   -4.999999751e+00  -4.999999747e+00  1.2e-08  0.13  
7   5.6e-08  1.8e-09  2.6e-14  1.00e+00   -4.999999971e+00  -4.999999970e+00  1.4e-09  0.16  
8   3.7e-10  4.0e-09  1.5e-17  1.00e+00   -5.000000000e+00  -5.000000000e+00  9.7e-12  0.21  
Optimizer terminated. Time: 0.21    


````

Again, the result matches the theoretical prediction for maximally entangled
bipartite state with $10^{-6}$ absolute tolerance [klep2022Optimization](@cite)!

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

