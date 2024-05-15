


######################## Conversion to ODEs/SDEs/jump, etc ##############################

"""
    oderatelaw(rx; combinatoric_ratelaw=true)

Given a [`Reaction`](@ref), return the symbolic reaction rate law used in
generated ODEs for the reaction. Note, for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `k*X(t)^2*Y(t)*Z(t)`. For a reaction of
the form

`k, 2X+3Y --> Z`

the expression that is returned will be `k * (X(t)^2/2) * (Y(t)^3/6)`.

Notes:
- Allocates
- `combinatoric_ratelaw=true` uses factorial scaling factors in calculating the
    rate law, i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. If
    `combinatoric_ratelaw=false` then the ratelaw is `k*S^2`, i.e. the scaling
    factor is ignored.
"""
function oderatelaw(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate

    # if the stoichiometric coefficients are not integers error if asking to scale rates
    !all(s -> s isa Union{Integer, Symbolic}, substoich) &&
        (combinatoric_ratelaw == true) &&
        error("Non-integer stoichiometric coefficients require the combinatoric_ratelaw=false keyword to oderatelaw, or passing combinatoric_ratelaws=false to convert or ODEProblem.")

    if !only_use_rate
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i, stoich) in enumerate(substoich)
            combinatoric_ratelaw && (coef *= factorial(stoich))
            rl *= isequal(stoich, one(stoich)) ? substrates[i] : substrates[i]^stoich
        end
        combinatoric_ratelaw && (!isequal(coef, one(coef))) && (rl /= coef)
    end
    rl
end

function assemble_oderhs(rs, ispcs; combinatoric_ratelaws = true, remove_conserved = false)
    nps = get_networkproperties(rs)
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    rhsvec = Any[0 for _ in ispcs]
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    for rx in get_rxs(rs)
        rl = oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)
        remove_conserved && (rl = substitute(rl, depspec_submap))
        for (spec, stoich) in rx.netstoich
            # dependent species don't get an ODE, so are skipped
            remove_conserved && (spec in nps.depspecs) && continue

            # constant or BC species also do not get equations
            drop_dynamics(spec) && continue

            i = species_to_idx[spec]
            if _iszero(rhsvec[i])
                if stoich isa Symbolic
                    rhsvec[i] = stoich * rl
                else
                    signedrl = (stoich > zero(stoich)) ? rl : -rl
                    rhsvec[i] = isone(abs(stoich)) ? signedrl : stoich * rl
                end
            else
                if stoich isa Symbolic
                    rhsvec[i] += stoich * rl
                else
                    Δspec = isone(abs(stoich)) ? rl : abs(stoich) * rl
                    rhsvec[i] = (stoich > zero(stoich)) ? (rhsvec[i] + Δspec) :
                                (rhsvec[i] - Δspec)
                end
            end
        end
    end

    rhsvec
end

function assemble_drift(rs, ispcs; combinatoric_ratelaws = true, as_odes = true,
                        include_zero_odes = true, remove_conserved = false)
    rhsvec = assemble_oderhs(rs, ispcs; combinatoric_ratelaws, remove_conserved)
    if as_odes
        D = Differential(get_iv(rs))
        eqs = [Equation(D(x), rhs)
               for (x, rhs) in zip(ispcs, rhsvec) if (include_zero_odes || (!_iszero(rhs)))]
    else
        eqs = [Equation(0, rhs) for rhs in rhsvec if (include_zero_odes || (!_iszero(rhs)))]
    end
    eqs
end

# this doesn't work with constraint equations currently
function assemble_diffusion(rs, sts, ispcs; combinatoric_ratelaws = true,
                            remove_conserved = false)
    # as BC species should ultimately get an equation, we include them in the noise matrix
    num_bcsts = count(isbc, get_unknowns(rs))

    # we make a matrix sized by the number of reactions
    eqs = Matrix{Any}(undef, length(sts) + num_bcsts, length(get_rxs(rs)))
    eqs .= 0
    species_to_idx = Dict((x => i for (i, x) in enumerate(ispcs)))
    nps = get_networkproperties(rs)
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    for (j, rx) in enumerate(get_rxs(rs))
        rlsqrt = sqrt(abs(oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)))
        has_noise_scaling(rx) && (rlsqrt *= get_noise_scaling(rx))
        remove_conserved && (rlsqrt = substitute(rlsqrt, depspec_submap))

        for (spec, stoich) in rx.netstoich
            # dependent species don't get an equation
            remove_conserved && (spec in nps.depspecs) && continue

            # constant or BC species also do not get equations
            drop_dynamics(spec) && continue

            i = species_to_idx[spec]
            if stoich isa Symbolic
                eqs[i, j] = stoich * rlsqrt
            else
                signedrlsqrt = (stoich > zero(stoich)) ? rlsqrt : -rlsqrt
                eqs[i, j] = isone(abs(stoich)) ? signedrlsqrt : stoich * rlsqrt
            end
        end
    end
    eqs
end

"""
    jumpratelaw(rx; combinatoric_ratelaw=true)

Given a [`Reaction`](@ref), return the symbolic reaction rate law used in
generated stochastic chemical kinetics model SSAs for the reaction. Note,
for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `k*X^2*Y*Z`. For a reaction of
the form

`k, 2X+3Y --> Z`

the expression that is returned will be `k * binomial(X,2) *
binomial(Y,3)`.

Notes:
- Allocates
- `combinatoric_ratelaw=true` uses binomials in calculating the rate law, i.e. for `2S ->
  0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaw=false` then
  the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling
  factor.
"""
function jumpratelaw(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    if !only_use_rate
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i, stoich) in enumerate(substoich)
            s = substrates[i]
            if stoich isa Symbolic
                rl *= combinatoric_ratelaw ? binomial(s, stoich) :
                      factorial(s) / factorial(s - stoich)
            else
                rl *= s
                isone(stoich) && continue
                for i in one(stoich):(stoich - one(stoich))
                    rl *= (s - i)
                end
                combinatoric_ratelaw && (coef *= factorial(stoich))
            end
        end
        combinatoric_ratelaw && !isequal(coef, one(coef)) && (rl /= coef)
    end
    rl
end

# if haveivdep=false then time dependent rates will still be classified as mass action
"""
```julia
ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                              haveivdep = nothing,
                              unknownset = Set(unknowns(rs)),
                              ivset = nothing)
```

True if a given reaction is of mass action form, i.e. `rx.rate` does not depend
on any chemical species that correspond to unknowns of the system, and does not
depend explicitly on the independent variable (usually time).

# Arguments
- `rx`, the [`Reaction`](@ref).
- `rs`, a [`ReactionSystem`](@ref) containing the reaction.
- Optional: `rxvars`, `Variable`s which are not in `rxvars` are ignored as
  possible dependencies.
- Optional: `haveivdep`, `true` if the [`Reaction`](@ref) `rate` field
  explicitly depends on any independent variable (i.e. t or for spatial systems x,y,etc).
  If not set, will be automatically calculated.
- Optional: `unknownset`, set of unknowns which if the rxvars are within mean rx is
  non-mass action.
- Optional: `ivset`, a `Set` of the independent variables of the system. If not provided and
  the system is spatial, i.e. `isspatial(rs) == true`, it will be created with all the
  spatial variables and the time variable. If the rate expression contains any element of
  `ivset`, then `ismassaction(rx,rs) == false`. Pass a custom set to control this behavior.

Notes:
- Non-integer stoichiometry is treated as non-mass action. This includes
  symbolic variables/terms or floating point numbers for stoichiometric
  coefficients.
"""
function ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                      haveivdep::Union{Nothing, Bool} = nothing,
                      unknownset = Set(get_unknowns(rs)), ivset = nothing)

    # we define non-integer (i.e. float or symbolic) stoich to be non-mass action
    ((eltype(rx.substoich) <: Integer) && (eltype(rx.prodstoich) <: Integer)) ||
        return false

    # if no dependencies must be zero order
    (length(rxvars) == 0) && return true

    if (haveivdep === nothing)
        if isspatial(rs)
            if (ivset === nothing)
                ivs = Set(get_sivs(rs))
                push!(ivs, get_iv(rs))
                ivdep = any(var -> var ∈ ivs, rxvars)
            else
                ivdep = any(var -> var ∈ ivset, rxvars)
            end
        else
            ivdep = any(var -> isequal(get_iv(rs), var), rxvars)
        end
        ivdep && return false
    else
        haveivdep && return false
    end
    rx.only_use_rate && return false
    @inbounds for var in rxvars
        # not mass action if have a non-constant, variable species in the rate expression
        (var in unknownset) && return false
    end

    return true
end

@inline function makemajump(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, netstoich = rx
    zeroorder = (length(substoich) == 0)
    reactant_stoch = Vector{Pair{Any, eltype(substoich)}}()
    @inbounds for (i, spec) in enumerate(substrates)
        # move constant species into the rate
        if isconstant(spec)
            rate *= spec
            isone(substoich[i]) && continue
            for i in 1:(substoich[i] - 1)
                rate *= spec - i
            end
        else
            push!(reactant_stoch, substrates[i] => substoich[i])
        end
    end

    if (!zeroorder) && combinatoric_ratelaw
        coef = prod(factorial, substoich)
        (!isone(coef)) && (rate /= coef)
    end

    net_stoch = filter(p -> !drop_dynamics(p[1]), netstoich)
    isempty(net_stoch) &&
        error("$rx has no net stoichiometry change once accounting for constant and boundary condition species. This is not supported.")

    MassActionJump(Num(rate), reactant_stoch, net_stoch, scale_rates = false,
                   useiszero = false)
end

# get_depgraph(rs)[i] is the list of reactions with rates depending on species changed by
# i'th reaction.
function get_depgraph(rs)
    jdeps = asgraph(rs)
    vdeps = variable_dependencies(rs)
    eqeq_dependencies(jdeps, vdeps).fadjlist
end

# recursively visit each neighbor's rooted tree and mark everything in it as vrj
function dfs_mark!(isvrjvec, visited, depgraph, i)
    visited[i] = true
    nhbrs = depgraph[i]
    for nhbr in nhbrs
        if !visited[nhbr]
            isvrjvec[nhbr] = true
            dfs_mark!(isvrjvec, visited, depgraph, nhbr)
        end
    end
    nothing
end

function assemble_jumps(rs; combinatoric_ratelaws = true)
    meqs = MassActionJump[]
    ceqs = ConstantRateJump[]
    veqs = VariableRateJump[]
    unknownset = Set(get_unknowns(rs))
    all(isspecies, unknownset) ||
        error("Conversion to JumpSystem currently requires all unknowns to be species.")
    rxvars = []

    isempty(get_rxs(rs)) &&
        error("Must give at least one reaction before constructing a JumpSystem.")

    # first we determine vrjs with an explicit time-dependent rate
    rxs = get_rxs(rs)
    isvrjvec = falses(length(rxs))
    havevrjs = false
    for (i, rx) in enumerate(rxs)
        empty!(rxvars)
        (rx.rate isa Symbolic) && get_variables!(rxvars, rx.rate)
        @inbounds for rxvar in rxvars
            if isequal(rxvar, get_iv(rs))
                isvrjvec[i] = true
                havevrjs = true
                break
            end
        end
    end

    # now we determine vrj's that depend on species modified by a previous vrj
    if havevrjs
        depgraph = get_depgraph(rs)
        visited = falses(length(isvrjvec))
        for (i, isvrj) in enumerate(isvrjvec)
            if isvrj && !visited[i]
                # dfs from the vrj node to propagate vrj classification
                dfs_mark!(isvrjvec, visited, depgraph, i)
            end
        end
    end

    for (i, rx) in enumerate(rxs)
        empty!(rxvars)
        (rx.rate isa Symbolic) && get_variables!(rxvars, rx.rate)

        isvrj = isvrjvec[i]
        if (!isvrj) && ismassaction(rx, rs; rxvars = rxvars, haveivdep = false,
                        unknownset = unknownset)
            push!(meqs, makemajump(rx; combinatoric_ratelaw = combinatoric_ratelaws))
        else
            rl = jumpratelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)
            affect = Vector{Equation}()
            for (spec, stoich) in rx.netstoich
                # don't change species that are constant or BCs
                (!drop_dynamics(spec)) && push!(affect, spec ~ spec + stoich)
            end
            if isvrj
                push!(veqs, VariableRateJump(rl, affect))
            else
                push!(ceqs, ConstantRateJump(rl, affect))
            end
        end
    end
    vcat(meqs, ceqs, veqs)
end

# merge constraint components with the ReactionSystem components
# also handles removing BC and constant species
function addconstraints!(eqs, rs::ReactionSystem, ists, ispcs; remove_conserved = false)
    # if there are BC species, put them after the independent species
    rssts = get_unknowns(rs)
    sts = any(isbc, rssts) ? vcat(ists, filter(isbc, rssts)) : ists
    ps = get_ps(rs)

    # make dependent species observables and add conservation constants as parameters
    if remove_conserved
        nps = get_networkproperties(rs)

        # add the conservation constants as parameters and set their values
        ps = vcat(ps, collect(eq.lhs for eq in nps.constantdefs))
        defs = copy(MT.defaults(rs))
        for eq in nps.constantdefs
            defs[eq.lhs] = eq.rhs
        end

        # add the dependent species as observed
        obs = copy(MT.observed(rs))
        append!(obs, nps.conservedeqs)
    else
        defs = MT.defaults(rs)
        obs = MT.observed(rs)
    end

    ceqs = Equation[eq for eq in get_eqs(rs) if eq isa Equation]
    if !isempty(ceqs)
        if remove_conserved
            @info """
                  Be careful mixing constraints and elimination of conservation laws.
                  Catalyst does not check that the conserved equations still hold for the
                  final coupled system of equations. Consider using `remove_conserved =
                  false` and instead calling ModelingToolkit.structural_simplify to simplify
                  any generated ODESystem or NonlinearSystem.
                  """
        end
        append!(eqs, ceqs)
    end

    eqs, sts, ps, obs, defs
end

# used by flattened systems that don't support constraint equations currently
function error_if_constraints(::Type{T}, sys::ReactionSystem) where {T <: MT.AbstractSystem}
    any(eq -> eq isa Equation, get_eqs(sys)) &&
        error("Can not convert to a system of type ", T,
              " when there are constraint equations.")
    nothing
end

# used by flattened systems that don't support differential equation constraint eqs
function error_if_constraint_odes(::Type{T},
                                  rs::ReactionSystem) where {T <: MT.AbstractSystem}
    any(eq -> (eq isa Equation) && MT.isdiffeq(eq), get_eqs(rs)) &&
        error("Cannot convert to system type $T when there are ODE constraint equations.")
    nothing
end

function spatial_convert_err(rs::ReactionSystem, systype)
    isspatial(rs) && error("Conversion to $systype is not supported for spatial networks.")
end

COMPLETENESS_ERROR = "A ReactionSystem must be complete before it can be converted to other system types. A ReactionSystem can be marked as complete using the `complete` function."

"""
```julia
Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
```
Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.ODESystem`.

Keyword args and default values:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
"""
function Base.convert(::Type{<:ODESystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, remove_conserved = false, checks = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, ODESystem)
    fullrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(fullrs)
    ists, ispcs = get_indep_sts(fullrs, remove_conserved)
    eqs = assemble_drift(fullrs, ispcs; combinatoric_ratelaws, remove_conserved,
                         include_zero_odes)
    eqs, us, ps, obs, defs = addconstraints!(eqs, fullrs, ists, ispcs; remove_conserved)

    ODESystem(eqs, get_iv(fullrs), us, ps;
              observed = obs,
              name,
              defaults = _merge(defaults,defs),
              checks,
              continuous_events = MT.get_continuous_events(fullrs),
              discrete_events = MT.get_discrete_events(fullrs),
              kwargs...)
end

"""
```julia
Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.NonlinearSystem`.

Keyword args and default values:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
"""
function Base.convert(::Type{<:NonlinearSystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, remove_conserved = false, checks = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      all_differentials_permitted = false, kwargs...)
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, NonlinearSystem)
    fullrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(fullrs)
    ists, ispcs = get_indep_sts(fullrs, remove_conserved)
    eqs = assemble_drift(fullrs, ispcs; combinatoric_ratelaws, remove_conserved,
                         as_odes = false, include_zero_odes)
    eqs, us, ps, obs, defs = addconstraints!(eqs, fullrs, ists, ispcs; remove_conserved)

    # Throws a warning if there are differential equations in non-standard format.
    # Next, sets all differential terms to `0`.
    all_differentials_permitted || nonlinear_convert_differentials_check(rs)
    eqs = [remove_diffs(eq.lhs) ~ remove_diffs(eq.rhs) for eq in eqs]


    NonlinearSystem(eqs, us, ps;
                    name,
                    observed = obs,
                    defaults = _merge(defaults,defs),
                    checks,
                    kwargs...)
end

# Finds and differentials in an expression, and sets these to 0.
function remove_diffs(expr)
    if Symbolics._occursin(Symbolics.is_derivative, expr)
        return Symbolics.replace(expr, diff_2_zero)
    else
        return expr
    end
end
diff_2_zero(expr) = (Symbolics.is_derivative(expr) ? 0.0 : expr)

# Ideally, when `ReactionSystem`s are converted to `NonlinearSystem`s, any coupled ODEs should be 
# on the form D(X) ~ ..., where lhs is the time derivative w.r.t. a single variable, and the rhs
# does not contain any differentials. If this is not teh case, we throw a warning to let the user
# know that they should be careful.
function nonlinear_convert_differentials_check(rs::ReactionSystem)
    for eq in filter(is_diff_equation, equations(rs))
        # For each differential equation, checks in order: 
        # If there is a differential on the right hand side.
        # If the lhs is not on the form D(...).
        # If the lhs upper level function is not a differential w.r.t. time.
        # If the contenct of the differential is not a variable (and nothing more).
        # If either of this is a case, throws the warning.
        if Symbolics._occursin(Symbolics.is_derivative, eq.rhs) ||
                !Symbolics.istree(eq.lhs) ||
                !isequal(Symbolics.operation(eq.lhs), Differential(get_iv(rs))) || 
                (length(arguments(eq.lhs)) != 1) ||
                !any(isequal(arguments(eq.lhs)[1]), nonspecies(rs))
            error("You are attempting to convert a `ReactionSystem` coupled with differential equations to a `NonlinearSystem`. However, some of these differentials are not of the form `D(x) ~ ...` where: 
                    (1) The left-hand side is a differential of a single variable with respect to the time independent variable, and
                    (2) The right-hand side does not contain any differentials.
                This is generally not permitted.
                
                If you still would like to perform this conversions, please use the `all_differentials_permitted = true` option. In this case, all differential will be set to `0`. 
                However, it is recommended to proceed with caution to ensure that the produced nonlinear equation makes sense for you intended application."
                )
        end
    end
end

"""
```julia
Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.SDESystem`.

Notes:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
- Does not currently support `ReactionSystem`s that include coupled algebraic or
  differential equations.
"""
function Base.convert(::Type{<:SDESystem}, rs::ReactionSystem;
                      name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, checks = false, remove_conserved = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, SDESystem)

    flatrs = Catalyst.flatten(rs)

    remove_conserved && conservationlaws(flatrs)
    ists, ispcs = get_indep_sts(flatrs, remove_conserved)
    eqs = assemble_drift(flatrs, ispcs; combinatoric_ratelaws, include_zero_odes,
                         remove_conserved)
    noiseeqs = assemble_diffusion(flatrs, ists, ispcs; combinatoric_ratelaws, remove_conserved)
    eqs, us, ps, obs, defs = addconstraints!(eqs, flatrs, ists, ispcs; remove_conserved)

    if any(isbc, get_unknowns(flatrs))
        @info "Boundary condition species detected. As constraint equations are not currently supported when converting to SDESystems, the resulting system will be undetermined. Consider using constant species instead."
    end

    SDESystem(eqs, noiseeqs, get_iv(flatrs), us, ps;
              observed = obs,
              name,
              defaults = defs,
              checks,
              continuous_events = MT.get_continuous_events(flatrs),
              discrete_events = MT.get_discrete_events(flatrs),
              kwargs...)
end

"""
```julia
Base.convert(::Type{<:JumpSystem},rs::ReactionSystem; combinatoric_ratelaws=true)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.JumpSystem`.

Notes:
- `combinatoric_ratelaws=true` uses binomials in calculating the rate law, i.e. for `2S ->
  0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaws=false` then
  the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling factor.
  Defaults to the value given when the `ReactionSystem` was constructed (which itself
  defaults to true).
- Does not currently support `ReactionSystem`s that include coupled algebraic or
  differential equations.
- Does not currently support continuous events as these are not supported by
  `ModelingToolkit.JumpSystems`.
"""
function Base.convert(::Type{<:JumpSystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      remove_conserved = nothing, checks = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, JumpSystem)

    (remove_conserved !== nothing) &&
        error("Catalyst does not support removing conserved species when converting to JumpSystems.")

    flatrs = Catalyst.flatten(rs)
    error_if_constraints(JumpSystem, flatrs)

    (length(MT.continuous_events(flatrs)) > 0) &&
        (@warn "continuous_events will be dropped as they are not currently supported by JumpSystems.")

    eqs = assemble_jumps(flatrs; combinatoric_ratelaws)

    # handle BC species
    sts, ispcs = get_indep_sts(flatrs)
    any(isbc, get_unknowns(flatrs)) && (sts = vcat(sts, filter(isbc, get_unknowns(flatrs))))
    ps = get_ps(flatrs)

    JumpSystem(eqs, get_iv(flatrs), sts, ps;
               observed = MT.observed(flatrs),
               name,
               defaults = _merge(defaults,MT.defaults(flatrs)),
               checks,
               discrete_events = MT.discrete_events(flatrs),
               kwargs...)
end

### Converts a reaction system to ODE or SDE problems ###

# ODEProblem from AbstractReactionNetwork
function DiffEqBase.ODEProblem(rs::ReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               check_length = false, name = nameof(rs),
                               combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                               include_zero_odes = true, remove_conserved = false,
                               checks = false, structural_simplify = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    osys = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
                   remove_conserved)

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify 
        (osys = MT.structural_simplify(osys))
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify=true` within `ODEProblem` call.")
    else
        osys = complete(osys)
    end
    
    return ODEProblem(osys, u0map, tspan, pmap, args...; check_length, kwargs...)
end

# NonlinearProblem from AbstractReactionNetwork
function DiffEqBase.NonlinearProblem(rs::ReactionSystem, u0,
                                     p = DiffEqBase.NullParameters(), args...;
                                     name = nameof(rs), include_zero_odes = true,
                                     combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                     remove_conserved = false, checks = false,
                                     check_length = false, all_differentials_permitted = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    nlsys = convert(NonlinearSystem, rs; name, combinatoric_ratelaws, include_zero_odes,
                    checks, all_differentials_permitted, remove_conserved)
    nlsys = complete(nlsys)
    return NonlinearProblem(nlsys, u0map, pmap, args...; check_length, 
                            kwargs...)
end

# SDEProblem from AbstractReactionNetwork
function DiffEqBase.SDEProblem(rs::ReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                               include_zero_odes = true, checks = false, check_length = false,
                               remove_conserved = false, structural_simplify = false, kwargs...)

    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    sde_sys = convert(SDESystem, rs; name, combinatoric_ratelaws,
                      include_zero_odes, checks, remove_conserved)

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify 
        (sde_sys = MT.structural_simplify(sde_sys))
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify=true` within `ODEProblem` call.")
    else
        sde_sys = complete(sde_sys)
    end
    
    p_matrix = zeros(length(get_unknowns(sde_sys)), numreactions(rs))
    return SDEProblem(sde_sys, u0map, tspan, pmap, args...; check_length,
                      noise_rate_prototype = p_matrix, kwargs...)
end

# DiscreteProblem from AbstractReactionNetwork
function DiffEqBase.DiscreteProblem(rs::ReactionSystem, u0, tspan::Tuple,
                                    p = DiffEqBase.NullParameters(), args...;
                                    name = nameof(rs),
                                    combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                    checks = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    jsys = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    jsys = complete(jsys)
    return DiscreteProblem(jsys, u0map, tspan, pmap, args...; kwargs...)
end

# JumpProblem from AbstractReactionNetwork
function JumpProcesses.JumpProblem(rs::ReactionSystem, prob, aggregator, args...;
                                   name = nameof(rs),
                                   combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                   checks = false, kwargs...)
    jsys = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    jsys = complete(jsys)
    return JumpProblem(jsys, prob, aggregator, args...; kwargs...)
end

# SteadyStateProblem from AbstractReactionNetwork
function DiffEqBase.SteadyStateProblem(rs::ReactionSystem, u0,
                                       p = DiffEqBase.NullParameters(), args...;
                                       check_length = false, name = nameof(rs),
                                       combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                       remove_conserved = false, include_zero_odes = true,
                                       checks = false, structural_simplify = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    osys = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
                   remove_conserved)

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify 
        (osys = MT.structural_simplify(osys))
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify=true` within `ODEProblem` call.")
    else
        osys = complete(osys)
    end

    return SteadyStateProblem(osys, u0map, pmap, args...; check_length, kwargs...)
end