# Hodgkin-Huxley Equation
As an illustration of adding constraint equations to Catalyst
[`ReactionSystem`](@ref)s, let's see how to build a Hodgkin-Huxley neuronal
dynamics model.

!!!
    This tutorial uses features that require at least Julia version 1.8.

```@example hheqs
using Catalyst, ModelingToolkit, DifferentialEquations, Plots

@variables t

function αₙ(V)
    theta = (V + 60) / 10
    ifelse(theta == 0.0, .1, .1*theta / (1 - exp(-theta)))
end
βₙ(V) = .125 * exp(-(V + 70)/80)

function Kchannel(; name, V)
	@parameters ḡK = 36.0 EK = -82.0
    openrx = @reaction αₙ($V), n′ --> n
    closerx = @reaction βₙ($V), n --> n′
    Iₖ = ḡK * n^4 * (V - EK)
end

function αₘ(V)
    theta = (V + 45) / 10
    IfElse.ifelse(theta == 0.0, 1.0, theta/(1 - exp(-theta)))
end
βₘ(V) = 4*exp(-(V + 70)/18)

αₕ(V) = .07 * exp(-(V + 70)/20)
βₕ(V) = 1/(1 + exp(-(V + 40)/10))

```