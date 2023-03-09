# Hodgkin-Huxley Equation
As an illustration of adding constraint equations to Catalyst
[`ReactionSystem`](@ref)s, let's see how to build a Hodgkin-Huxley neuronal
dynamics model.

!!!
    This tutorial uses features that require at least Julia version 1.8.

```@example hheqs
using Catalyst, ModelingToolkit, DifferentialEquations, Plots

@variables t

# Potassium channels
function αₙ(V)
    theta = (V + 60) / 10
    ifelse(theta == 0.0, .1, .1*theta / (1 - exp(-theta)))
end
βₙ(V) = .125 * exp(-(V + 70)/80)

function Kchannel(; name, V)
    gatingrxs = @reaction_network $name begin
	    (αₙ($V),βₙ($V)), n′ <--> n
    end

    # equation for K current
	@parameters ḡK = 36.0 EK = -82.0
    @variables Iₖ(t)
    @unpack n = gatingrxs

    @named current = ODESystem([Iₖ ~ ḡK * n^4 * (V - EK)], t)
    extend(current, gatingrxs; name = name)
end

# Sodium channels
function αₘ(V)
    theta = (V + 45) / 10
    ifelse(theta == 0.0, 1.0, theta/(1 - exp(-theta)))
end
βₘ(V) = 4*exp(-(V + 70)/18)
αₕ(V) = .07 * exp(-(V + 70)/20)
βₕ(V) = 1/(1 + exp(-(V + 40)/10))

function Nachannel(; name, V)
    gatingrxs = @reaction_network $name begin
	    (αₘ($V),βₘ($V)), m′ <--> m
	    (αₕ($V),βₕ($V)), h′ <--> h
    end

    # equation for Na current
	@parameters ḡNa = 120.0 ENa = 45.0
    @variables Iₙₐ(t)
    @unpack m,h = gatingrxs

    @named current = ODESystem([Iₙₐ ~ ḡNa*m^3*h*(V-ENa)], t)
    extend(current, gatingrxs; name = name)
end

# leakage current
function leakchannel(; name, V)
    @parameters ḡL = 0.3 EL = -59.0
    @variables Iₗ(t)
    ODESystem([Iₗ ~ ḡL*(V-EL)], t; name)
end

# HH equation

```