# Parametric Stoichiometry
Catalyst supports making having stoichiometric coefficients that involve parameters, species, or even general expressions. In this tutorial we show several examples of how to use parametric stoichiometry, and discuss several caveats to be aware of.

Let's first consider a simple reversible reaction where the number of reactants is a parameter and the number of products is the product of two parameters. Note, currently Catalyst's `@reaction_network` macro does not support stoichiometric parameters, so they need to be specified through the symbolic API interface:
```julia
using Catalyst, Latexify, DifferentialEquations, Plots
@parameters k₊,k₋,m,n
@variables t, A(t), B(t)
rxs = [Reaction(k₊,[A],[B],[m],[m*n]),
       Reaction(k₋,[B],[A])]
@named revsys = ReactionSystem(rxs,t)
latexify(revsys)
```
giving
```math

``` 