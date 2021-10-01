using Catalyst, LinearAlgebra, OrdinaryDiffEq, Test, NonlinearSolve

# Repressilator model
@parameters t α₀ α K n δ β μ
@variables m(t) P(t) R(t)
rxs = [
       Reaction(α₀, nothing, [m]),
       Reaction(α / (1 + (R/K)^n), nothing, [m]),
       Reaction(δ, [m], nothing),
       Reaction(β, [m], [m,P]),
       Reaction(μ, [P], nothing)
      ]

specs = [m,P,R]
pars  = [α₀,α,K,n,δ,β,μ]
@named rs = ReactionSystem(rxs, t, specs, pars)

# using ODESystem components
@named sys₁ = convert(ODESystem, rs; include_zero_odes=false)
@named sys₂ = convert(ODESystem, rs; include_zero_odes=false)
@named sys₃ = convert(ODESystem, rs; include_zero_odes=false)
connections = [sys₁.R ~ sys₃.P,
               sys₂.R ~ sys₁.P,
               sys₃.R ~ sys₂.P]
@named connected = ODESystem(connections, t, [], [], systems=[sys₁,sys₂,sys₃])
oderepressilator = structural_simplify(connected)

pvals = [sys₁.α₀ => 5e-4,
         sys₁.α => .5,
         sys₁.K => 40.0,
         sys₁.n => 2,
         sys₁.δ => (log(2)/120),
         sys₁.β => (20*log(2)/120),
         sys₁.μ => (log(2)/600),
         sys₂.α₀ => 5e-4,
         sys₂.α => .5,
         sys₂.K => 40.0,
         sys₂.n => 2,
         sys₂.δ => (log(2)/120),
         sys₂.β => (20*log(2)/120),
         sys₂.μ => (log(2)/600),
         sys₃.α₀ => 5e-4,
         sys₃.α => .5,
         sys₃.K => 40.0,
         sys₃.n => 2,
         sys₃.δ => (log(2)/120),
         sys₃.β => (20*log(2)/120),
         sys₃.μ => (log(2)/600)]
u₀    = [sys₁.m => 0.0, sys₁.P => 20.0, sys₁.R => 0.0, sys₂.m => 0.0, sys₂.P => 0.0, sys₂.R => 0.0, sys₃.m => 0.0, sys₃.P => 0.0, sys₃.R => 0.0]
tspan = (0.0, 100000.0)
oprob = ODEProblem(oderepressilator, u₀, tspan, pvals)
sol = solve(oprob, Tsit5())

# hardcoded network
function repress!(f, y, p, t)
    α = p.α; α₀ = p.α₀; β = p.β; δ = p.δ; μ = p.μ; K = p.K; n = p.n
    f[1] = α / (1 + (y[6] / K)^n) - δ * y[1] + α₀
    f[2] = α / (1 + (y[4] / K)^n) - δ * y[2] + α₀
    f[3] = α / (1 + (y[5] / K)^n) - δ * y[3] + α₀
    f[4] = β * y[1] - μ * y[4]
    f[5] = β * y[2] - μ * y[5]
    f[6] = β * y[3] - μ * y[6]
    nothing
end
ps = (α₀=5e-4, α=.5, K=40.0, n=2, δ=(log(2)/120), β=(20*log(2)/120), μ=(log(2)/600))
u0 = [0.0,0.0,0.0,20.0,0.0,0.0]
oprob2 = ODEProblem(repress!, u0, tspan, ps)
sol2 = solve(oprob2, Tsit5())
tvs = 0:1:tspan[end]
@test all(isapprox.(sol(tvs,idxs=sys₁.P),sol2(tvs, idxs=4),atol=1e-4))

# using ReactionSystem components
@named sys₁ = ReactionSystem(rxs, t, specs, pars)
@named sys₂ = ReactionSystem(rxs, t, specs, pars)
@named sys₃ = ReactionSystem(rxs, t, specs, pars)
connections = [ParentScope(sys₁.R) ~ ParentScope(sys₃.P),
               ParentScope(sys₂.R) ~ ParentScope(sys₁.P),
               ParentScope(sys₃.R) ~ ParentScope(sys₂.P)]
@named csys = ODESystem(connections, t, [], [])
@named repressilator = ReactionSystem(t; systems=[csys,sys₁,sys₂,sys₃])
@named oderepressilator2 = convert(ODESystem, repressilator, include_zero_odes=false)
sys2 = structural_simplify(oderepressilator2)  # FAILS currently
oprob = ODEProblem(sys2, u₀, tspan, pvals)
sol = solve(oprob, Tsit5())
@test all(isapprox.(sol(tvs,idxs=sys₁.P),sol2(tvs, idxs=4),atol=1e-4))

# test conversion to nonlinear system 
@named nsys = NonlinearSystem(connections, [], [])
@named ssrepressilator = ReactionSystem(t; systems=[nsys,sys₁,sys₂,sys₃])
@named nlrepressilator = convert(NonlinearSystem, ssrepressilator, include_zero_odes=false)
sys2 = structural_simplify(nlrepressilator)
@test length(equations(sys2)) == 6
nlprob = NonlinearProblem(sys2, u₀, pvals)
sol = solve(nlprob, NewtonRaphson(), tol=1e-9)
@test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
@test sol[sys₁.m] ≈ sol[sys₂.m] ≈ sol[sys₃.m]
@test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]

# flattening
fsys = Catalyst.flatten(ssrepressilator)
@named nlrepressilator = convert(NonlinearSystem, fsys, include_zero_odes=false)
sys2 = structural_simplify(nlrepressilator)
@test length(equations(sys2)) == 6
nlprob = NonlinearProblem(sys2, u₀, pvals)
sol = solve(nlprob, NewtonRaphson(), tol=1e-9)
@test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
@test sol[sys₁.m] ≈ sol[sys₂.m] ≈ sol[sys₃.m]
@test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]

# test constraints
connections = [sys₁.R ~ sys₃.P,
               sys₂.R ~ sys₁.P,
               sys₃.R ~ sys₂.P]
@named csys = NonlinearSystem(connections, [sys₁.R,sys₃.P,sys₂.R,sys₁.P,sys₃.R,sys₂.P],[])
@named repressilator2 = ReactionSystem(t; constraints=csys, systems=[sys₁,sys₂,sys₃])
@named nlrepressilator = convert(NonlinearSystem, repressilator2, include_zero_odes=false)
sys2 = structural_simplify(nlrepressilator)
@test length(equations(sys2)) == 6
nlprob = NonlinearProblem(sys2, u₀, pvals)
sol = solve(nlprob, NewtonRaphson(), tol=1e-9)
@test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
@test sol[sys₁.m] ≈ sol[sys₂.m] ≈ sol[sys₃.m]
@test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]

# test can make ODESystem
@named oderepressilator = convert(ODESystem, repressilator2, include_zero_odes=false)
sys2 = structural_simplify(oderepressilator)  # FAILS currently
oprob = ODEProblem(sys2, u₀, tspan, pvals)
sol = solve(oprob, Tsit5())
@test all(isapprox.(sol(tvs,idxs=sys₁.P),sol2(tvs, idxs=4),atol=1e-4))

# test extending with NonlinearSystem
@named repressilator2 = ReactionSystem(t; systems=[sys₁,sys₂,sys₃])
repressilator2 = Catalyst.flatten(repressilator2)
repressilator2 = extend(csys, repressilator2)
@named nlrepressilator = convert(NonlinearSystem, repressilator2, include_zero_odes=false)
sys2 = structural_simplify(nlrepressilator)
@test length(equations(sys2)) == 6
nlprob = NonlinearProblem(sys2, u₀, pvals)
sol = solve(nlprob, NewtonRaphson(), tol=1e-9)
@test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
@test sol[sys₁.m] ≈ sol[sys₂.m] ≈ sol[sys₃.m]
@test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]


# TODO add conversion to SDE and JumpSystems once supported

# adding algebraic constraints
@parameters t, r₊, r₋, β
@variables A(t), B(t), C(t), D(t)
rxs1 = [Reaction(r₊, [A,B], [C])]
rxs2 = [Reaction(r₋, [C], [A,B])]
@named rs1 = ReactionSystem(rxs1, t, [A,B,C], [r₊])
@named rs2 = ReactionSystem(rxs2, t, [A,B,C], [r₋])
@named rs  = extend(rs1, rs2)
@test issetequal(states(rs), [A,B,C])
@test issetequal(parameters(rs), [r₊,r₋])
@test issetequal(equations(rs), union(rxs1,rxs2))
A2 = ModelingToolkit.ParentScope(A)
B2 = ModelingToolkit.ParentScope(B)
nseqs = [D ~ 2*A2 + β*B2]
@named ns = ODESystem(nseqs, t, [A2,B2,D], [β])
rs = compose(rs, [ns])
osys = convert(ODESystem, rs; include_zero_odes=false)
p  = [r₊ => 1.0, r₋ => 2.0, ns.β => 3.0]
u₀ = [A => 1.0, B => 2.0, C => 0.0]
oprob = ODEProblem(structural_simplify(osys), u₀, (0.0,10.0), p)
sol = solve(oprob, Tsit5())
@test isapprox(0, norm(sol[ns.D] .- 2*sol[A] - 3*sol[B]), atol=(100*eps())) 

# test API functions for composed model
@test issetequal(species(rs), [A,B,C])
@test issetequal(states(rs), [A,B,C,ns.D])
@test issetequal(reactionparams(rs), [r₊,r₋])
@test issetequal(parameters(rs), [r₊,r₋,ns.β])
@test issetequal(reactions(rs), union(rxs1,rxs2))
@test issetequal(filter(eq -> eq isa Reaction, equations(rs)), union(rxs1,rxs2))
@test issetequal(filter(eq -> eq isa Equation, equations(rs)), [ModelingToolkit.namespace_equation(nseqs[1],ns)])

# check several levels of nesting namespace and filter ok for the API functions
@parameters t, p1, p2a, p2b, p3a, p3b
@variables A1(t), A2a(t), A2b(t), A3a(t), A3b(t)
rxs1 = [Reaction(p1, [A1], nothing)]
rxs2 = [Reaction(p2a, [A2a], nothing), Reaction(p2b, [ParentScope(A1)], nothing)]
eqs2 = [ParentScope(A1) ~ ParentScope(p1)*A2b]
rxs3 = [Reaction(p3a, [A3a], nothing), Reaction(ParentScope(p2a), nothing, [ParentScope(A2a)])]
eqs3 = [ParentScope(A2a) ~ p3b*A3b]

@named rs3 = ReactionSystem(rxs3,t,[A3a,ParentScope(A2a)],[p3a,ParentScope(p2a)])
@named ns3 = NonlinearSystem(eqs3,[ParentScope(A2a),A3b],[p3b])
@named rs2 = ReactionSystem(rxs2,t,[A2a,ParentScope(A1)],[p2a,p2b], systems=[rs3,ns3])
@named ns2 = NonlinearSystem(eqs2,[ParentScope(A1),A2b],[ParentScope(p1)])
@named rs1 = ReactionSystem(rxs1,t,[A1],[p1],systems=[rs2,ns2])

@test issetequal(states(rs1), [A1,rs2.A2a,ns2.A2b,rs2.rs3.A3a,rs2.ns3.A3b])
@test issetequal(species(rs1), [A1,rs2.A2a,rs2.rs3.A3a])
@test issetequal(parameters(rs1), [p1,rs2.p2a,rs2.p2b,rs2.rs3.p3a,rs2.ns3.p3b])
@test issetequal(reactionparams(rs1), [p1, rs2.p2a, rs2.p2b, rs2.rs3.p3a])
rxs = [rx for rx in Iterators.filter(eq -> eq isa Reaction, equations(rs1))]
@test issubset(rxs, reactions(rs1)) && issubset(reactions(rs1), rxs)

# test throw error if there are ODE constraints and convert to NonlinearSystem
rn = @reaction_network rn begin
    (k1,k2), A <--> B
    end k1 k2
@parameters a,b
@variables t,A(t),C(t)
D = Differential(t)
eqs = [D(A) ~ -a*A + C, D(C) ~ -b*C + a*A]
@named osys = ODESystem(eqs,t,[A,C],[a,b])
rn2 = extend(osys,rn)
rnodes = convert(ODESystem,rn2)
@test_throws ErrorException convert(NonlinearSystem,rn2)

@variables G(t)
eqs = [D(G) ~ -G]
@named osys2 = ODESystem(eqs,t)
rn3 = compose(rn2, osys2)
@test length(equations(rn3)) == 5

# check conversions work with algebraic constraints
eqs = [0 ~ -a*A + C, 0 ~ -b*C + a*A]
@named nlsys = NonlinearSystem(eqs,[A,C],[a,b])
rn2 = extend(nlsys,rn)
rnodes = convert(ODESystem,rn2)
rnnlsys = convert(NonlinearSystem,rn2)
@named nlsys = ODESystem(eqs,t,[A,C],[a,b])
rn2 = extend(nlsys,rn)
rnodes = convert(ODESystem,rn2)
rnnlsys = convert(NonlinearSystem,rn2)