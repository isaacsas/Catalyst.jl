# Functions for querying network properties.

# DEPRECIATED after v9.0
function params(network)
    Base.depwarn("`params` is depreciated, please use `ModelingToolkit.parameters` for all system and subsystem parameters, or `reactionparams` for all parameters within system and subsystem `Reaction`s.", :params, force=true)
    parameters(network)
end

function numparams(network)
    Base.depwarn("`numparams` is depreciated, please use `length(ModelingToolkit.parameters)` for the total number of parameters across all systems and subsystems, or `numreactionparams` for the number of parameters within system and subsystem `Reaction`s.", :params, force=true)
    length(parameters(network))
end

"""
    merge(network1::ReactionSystem, network2::ReactionSystem)

Create a new network merging `network1` and `network2`.

Notes:
- Duplicate reactions between the two networks are not filtered out.
- [`Reaction`](@ref)s are not deepcopied to minimize allocations, so the new
  network will share underlying data arrays.
- Subsystems are not deepcopied between the two networks and will hence be
  shared.
- Returns the merged network.
"""
function Base.merge(network1::ReactionSystem, network2::ReactionSystem)
    Base.depwarn("`merge(sys1::ReactionSystem, sys2::ReactionNetwork)` is depreciated, please use `ModelingToolkit.extend` instead.", :merge, force=true)
    network = make_empty_network()
    merge!(network, network1)
    merge!(network, network2)
    network
end

######### Accessors: #########

function filter_nonrxsys(network)
    systems   = get_systems(network)
    rxsystems = ReactionSystem[]
    for sys in systems
        (sys isa ReactionSystem) && push!(rxsystems, sys)
    end
    rxsystems
end


function species(network::ReactionSystem, sts)
    [MT.renamespace(network,st) for st in sts]
end

"""
    species(network)

Given a [`ReactionSystem`](@ref), return a vector of all species defined in the system
and any subsystems that are of type `ReactionSystem`. To get the variables in the system
and all subsystems, including non-`ReactionSystem` subsystems, uses `states(network)`.

Notes:
- If `ModelingToolkit.get_systems(network)` is non-empty will allocate.
"""
function species(network)
    sts = get_states(network)
    systems = filter_nonrxsys(network)
    isempty(systems) && return sts
    unique([sts; reduce(vcat, map(sys -> species(sys,species(sys)), systems))])
end

"""
    reactionparams(network)

Given a [`ReactionSystem`](@ref), return a vector of all parameters defined
within the system and any subsystems that are of type `ReactionSystem`. To get
the parameters in the system and all subsystems, including non-`ReactionSystem`
subsystems, use `parameters(network)`.

Notes:
- If `ModelingToolkit.get_systems(network)` is non-empty will allocate.
"""
function reactionparams(network)
    ps = get_ps(network)
    systems = filter_nonrxsys(network)
    isempty(systems) && return ps
    unique([ps; reduce(vcat, map(sys -> species(sys,reactionparams(sys)), systems))])
end

function namespace_reactions(network::ReactionSystem)
    rxs = reactions(network)
    isempty(rxs) && return Reaction[]
    map(rx -> namespace_equation(rx, network), rxs)
end

"""
    reactions(network)

Given a [`ReactionSystem`](@ref), return a vector of all `Reactions` in the system.

Notes:
- If `ModelingToolkit.get_systems(network)` is not empty, will allocate.
"""
function reactions(network)    
    rxs = get_eqs(network)
    systems = filter_nonrxsys(network)
    isempty(systems) && (return rxs)
    [rxs; reduce(vcat, namespace_reactions.(systems); init=Reaction[])]
end

"""
    speciesmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping species that
participate in `Reaction`s to their index within [`species(network)`](@ref).
"""
function speciesmap(network)
    Dict(S => i for (i,S) in enumerate(species(network)))
end

"""
    paramsmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping from all
parameters that appear within the system to their index within
`parameters(network)`.
"""
function paramsmap(network)
    Dict(p => i for (i,p) in enumerate(parameters(network)))
end


"""
    reactionparamsmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping from parameters that
appear within `Reaction`s to their index within [`reactionparams(network)`](@ref).
"""
function reactionparamsmap(network)
    Dict(p => i for (i,p) in enumerate(reactionparams(network)))
end

"""
    numspecies(network)

Return the total number of species within the given [`ReactionSystem`](@ref) and
subsystems that are `ReactionSystem`s. 

Notes
- If there are no subsystems this will be fast.
- As this calls [`species`](@ref), it can be slow and will allocate if there are
  any subsystems. 
"""
function numspecies(network)
    length(species(network))
end

"""
    numreactions(network)

Return the total number of reactions within the given [`ReactionSystem`](@ref)
and subsystems that are `ReactionSystem`s.
"""
function numreactions(network)
    nr = length(get_eqs(network))
    for sys in get_systems(network)
        (sys isa ReactionSystem) && (nr += numreactions(sys))
    end
    nr
end

"""
    numreactionparams(network)

Return the total number of parameters within the given [`ReactionSystem`](@ref)
and subsystems that are `ReactionSystem`s.

Notes
- If there are no subsystems this will be fast.
- As this calls [`reactionparams`](@ref), it can be slow and will allocate if
  there are any subsystems. 
"""
function numreactionparams(network)
    length(reactionparams(network))
end

"""
    dependents(rx, network)

Given a [`Reaction`](@ref) and a [`ReactionSystem`](@ref), return a vector of
`ModelingToolkit.Num`s corresponding to species the *reaction rate
law* depends on. E.g., for

`k*W, 2X + 3Y --> 5Z + W`

the returned vector would be `[W(t),X(t),Y(t)]`.

Notes:
- Allocates
- Does not check for dependents within any subsystems.
"""
function dependents(rx, network)
    if rx.rate isa Number
        return rx.substrates
    else
        rvars = get_variables(rx.rate, species(network))
        return union!(rvars, rx.substrates)
    end
end

"""
    dependents(rx, network)

See documentation for [`dependents`](@ref).
"""
function dependants(rx, network)
    dependents(rx, network)
end

"""
    reactionrates(network)

Given a [`ReactionSystem`](@ref), returns a vector of the symbolic reaction
rates for each reaction.
"""
function reactionrates(rn)
    [r.rate for r in reactions(rn)]
end

"""
    substoichmat(rn; sparse=false, smap=speciesmap(rn))

Returns the substrate stoichiometry matrix, ``S``, with ``S_{i j}`` the
stoichiometric coefficient of the ith substrate within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
"""
function substoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; smap=speciesmap(rn))
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i,sub) in enumerate(rx.substrates)
            push!(Js, k)
            push!(Is, smap[sub])
            push!(Vs, stoich[i])
        end
    end
    sparse(Is,Js,Vs,numspecies(rn),numreactions(rn))
end
function substoichmat(::Type{Matrix{Int}},rn::ReactionSystem; smap=speciesmap(rn))
    smat = zeros(Int, numspecies(rn), numreactions(rn))
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i,sub) in enumerate(rx.substrates)
            smat[smap[sub],k] = stoich[i]
        end
    end
    smat
end
function substoichmat(rn::ReactionSystem; sparse::Bool=false, smap=speciesmap(rn))
    isempty(get_systems(rn)) || error("substoichmat does not currently support subsystems.")
	sparse ? substoichmat(SparseMatrixCSC{Int,Int}, rn; smap=smap) : substoichmat(Matrix{Int}, rn; smap=smap)
end


"""
    prodstoichmat(rn; sparse=false, smap=speciesmap(rn))

Returns the product stoichiometry matrix, ``P``, with ``P_{i j}`` the
stoichiometric coefficient of the ith product within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
"""
function prodstoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; smap=speciesmap(rn))
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i,prod) in enumerate(rx.products)
			push!(Js, k)
			push!(Is, smap[prod])
			push!(Vs, stoich[i])
        end
    end
    sparse(Is,Js,Vs,numspecies(rn),numreactions(rn))
end
function prodstoichmat(::Type{Matrix{Int}},rn::ReactionSystem; smap=speciesmap(rn))
    pmat = zeros(Int, numspecies(rn), numreactions(rn))
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i,prod) in enumerate(rx.products)
            pmat[smap[prod],k] = stoich[i]
        end
    end
    pmat
end
function prodstoichmat(rn::ReactionSystem; sparse=false, smap=speciesmap(rn))
    isempty(get_systems(rn)) || error("prodstoichmat does not currently support subsystems.")
	sparse ? prodstoichmat(SparseMatrixCSC{Int,Int}, rn; smap=smap) : prodstoichmat(Matrix{Int}, rn; smap=smap)
end


"""
    netstoichmat(rn, sparse=false; smap=speciesmap(rn))

Returns the net stoichiometry matrix, ``N``, with ``N_{i j}`` the net
stoichiometric coefficient of the ith species within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
"""
function netstoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; smap=speciesmap(rn))
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (k,rx) in pairs(reactions(rn))
        for (spec,coef) in rx.netstoich
			push!(Js, k)
			push!(Is, smap[spec])
			push!(Vs, coef)
        end
    end
    sparse(Is,Js,Vs,numspecies(rn),numreactions(rn))
end
function netstoichmat(::Type{Matrix{Int}},rn::ReactionSystem; smap=speciesmap(rn))
    nmat = zeros(Int,numspecies(rn),numreactions(rn))
    for (k,rx) in pairs(reactions(rn))
        for (spec,coef) in rx.netstoich
            nmat[smap[spec],k] = coef
        end
    end
    nmat
end
function netstoichmat(rn::ReactionSystem; sparse=false, smap=speciesmap(rn))
    isempty(get_systems(rn)) || error("netstoichmat does not currently support subsystems.")
	sparse ? netstoichmat(SparseMatrixCSC{Int,Int}, rn; smap=smap) : netstoichmat(Matrix{Int}, rn; smap=smap)
end

######################## reaction complexes and reaction rates ###############################
"""
$(TYPEDEF)
One reaction complex element

# Fields
$(FIELDS)
"""
struct ReactionComplexElement{T}
    """The integer id of the species representing this element."""
    speciesid::Int
    """The stoichiometric coefficient of this species."""
    speciesstoich::T
end

"""
$(TYPEDEF)
One reaction complex.

# Fields
$(FIELDS)
"""
struct ReactionComplex{V<:Integer} <: AbstractVector{ReactionComplexElement{V}}
    """The integer ids of all species participating in this complex."""
    speciesids::Vector{Int}
    """The stoichiometric coefficients of all species participating in this complex."""
    speciesstoichs::Vector{V}
end

function (==)(a::ReactionComplex{V},b::ReactionComplex{V}) where {V <: Integer} 
    (a.speciesids == b.speciesids) &&
    (a.speciesstoichs == b.speciesstoichs) 
end
hash(rc::ReactionComplex,h::UInt) = Base.hash(rc.speciesids,Base.hash(rc.speciesstoichs,h))
Base.size(rc::ReactionComplex) = size(rc.speciesids)
Base.length(rc::ReactionComplex) = length(rc.speciesids)
Base.getindex(rc::ReactionComplex, i...) = 
        ReactionComplexElement(getindex(rc.speciesids, i...), getindex(rc.speciesstoichs, i...))
Base.setindex!(rc::ReactionComplex, t::ReactionComplexElement, i...) = 
    (setindex!(rc.speciesids, t.speciesid, i...); setindex!(rc.speciesstoichs, t.speciesstoich, i...); rc) 
Base.isless(a::ReactionComplexElement, b::ReactionComplexElement) = isless(a.speciesid, b.speciesid)
Base.Sort.defalg(::ReactionComplex{T}) where {T <: Integer} = Base.DEFAULT_UNSTABLE

"""
    reactioncomplexmap(rn::ReactionSystem; smap=speciesmap(rn))

Find each [`ReactionComplex`](@ref) within the specified system, constructing a
mapping from the complex to vectors that indicate which reactions it appears in
as substrates and products.

Notes:
- Each [`ReactionComplex`](@ref) is mapped to a vector of pairs, with each pair
  having the form `reactionidx => ± 1`, where `-1` indicates the complex appears
  as a substrate and `+1` as a product in the reaction with integer label
  `reactionidx`.
"""
function reactioncomplexmap(rn::ReactionSystem; smap=speciesmap(rn))
    isempty(get_systems(rn)) || error("reactioncomplexmap does not currently support subsystems.")
    rxs = reactions(rn)
    numreactions(rn) > 0 || error("There must be at least one reaction to find reaction complexes.")
    complextorxsmap = OrderedDict{ReactionComplex{eltype(rxs[1].substoich)},Vector{Pair{Int,Int}}}()
    for (i,rx) in enumerate(rxs)
        reactantids = isempty(rx.substrates) ? Vector{Int}() : [smap[sub] for sub in rx.substrates]
        subrc = sort!(ReactionComplex(reactantids, copy(rx.substoich)))
        if haskey(complextorxsmap, subrc)
            push!(complextorxsmap[subrc], i => -1)
        else
            complextorxsmap[subrc] = [i => -1]
        end

        productids = isempty(rx.products) ? Vector{Int}() : [smap[prod] for prod in rx.products]
        prodrc = sort!(ReactionComplex(productids, copy(rx.prodstoich)))
        if haskey(complextorxsmap, prodrc)
            push!(complextorxsmap[prodrc], i => 1)
        else
            complextorxsmap[prodrc] = [i => 1]
        end
    end
    complextorxsmap
end


function reactioncomplexes(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem, smap, complextorxsmap)
    complexes = collect(keys(complextorxsmap))
    Is=Int[];  Js=Int[];  Vs=Int[];
	for (i,c) in enumerate(complexes)
        for (j,σ) in complextorxsmap[c]
			push!(Is, i)
			push!(Js, j)
			push!(Vs, σ)
        end
    end
	B = sparse(Is,Js,Vs,length(complexes),numreactions(rn))
    complexes,B
end
function reactioncomplexes(::Type{Matrix{Int}}, rn::ReactionSystem, smap, complextorxsmap)
    complexes = collect(keys(complextorxsmap))
    B = zeros(Int, length(complexes), numreactions(rn));
    for (i,c) in enumerate(complexes)
        for (j,σ) in complextorxsmap[c]
            B[i,j] = σ
        end
    end
    complexes,B
end

@doc raw"""
    reactioncomplexes(network::ReactionSystem; sparse=false, smap=speciesmap(rn), 
                      complextorxsmap=reactioncomplexmap(rn; smap=smap))

Calculate the reaction complexes and complex incidence matrix for the given
[`ReactionSystem`](@ref). 

Notes:
- returns a pair of a vector of [`ReactionComplex`](@ref)s and the complex
  incidence matrix.
- An empty [`ReactionComplex`](@ref) denotes the null (∅) state (from reactions
  like ∅ -> A or A -> ∅).
- The complex incidence matrix, ``B``, is number of complexes by number of
  reactions with
```math
B_{i j} = \begin{cases}
-1, &\text{if the i'th complex is the substrate of the j'th reaction},\\
1, &\text{if the i'th complex is the product of the j'th reaction},\\
0, &\text{otherwise.}
\end{cases}
```
- Set sparse=true for a sparse matrix representation of the incidence matrix
"""
function reactioncomplexes(rn::ReactionSystem; sparse=false, smap=speciesmap(rn), 
                           complextorxsmap=reactioncomplexmap(rn; smap=smap))
    isempty(get_systems(rn)) || error("reactioncomplexes does not currently support subsystems.")
	sparse ? reactioncomplexes(SparseMatrixCSC{Int,Int}, rn, smap, complextorxsmap) :
             reactioncomplexes(Matrix{Int}, rn, smap, complextorxsmap)
end


function complexstoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem, rcs)
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (i,rc) in enumerate(rcs)
        for rcel in rc
			push!(Is, rcel.speciesid)
			push!(Js, i)
			push!(Vs, rcel.speciesstoich)
        end
    end
    Z = sparse(Is,Js,Vs, numspecies(rn), length(rcs))
end
function complexstoichmat(::Type{Matrix{Int}}, rn::ReactionSystem, rcs)
    Z=zeros(Int, numspecies(rn), length(rcs))
    for (i,rc) in enumerate(rcs)
        for rcel in rc
            Z[rcel.speciesid,i] = rcel.speciesstoich
        end
    end
    Z
end

"""
    complexstoichmat(network::ReactionSystem; sparse=false, rcs=keys(reactioncomplexmap(rn)))

Given a [`ReactionSystem`](@ref) and vector of reaction complexes, return a
matrix with positive entries of size number of species by number of complexes,
where the non-zero positive entries in the kth column denote stoichiometric
coefficients of the species participating in the kth reaction complex.

Notes:
- `rcs` correspond to an iterable of the `ReactionComplexes`, i.e.
  `rcs=keys(reactioncomplexmap(rn))` or `reactioncomplexes(rn)[1]`.
- Set sparse=true for a sparse matrix representation
"""
function complexstoichmat(rn::ReactionSystem; sparse=false, rcs=keys(reactioncomplexmap(rn)))
    isempty(get_systems(rn)) || error("complexstoichmat does not currently support subsystems.")
    sparse ? complexstoichmat(SparseMatrixCSC{Int,Int}, rn, rcs) : 
             complexstoichmat(Matrix{Int}, rn, rcs)
end


function complexoutgoingmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem, B)    
    n = size(B,2)
	rows = rowvals(B)
	vals = nonzeros(B)
    Is = Int[]; Js = Int[]; Vs = Int[]
    sizehint!(Is, div(length(vals),2))
    sizehint!(Js, div(length(vals),2))
    sizehint!(Vs, div(length(vals),2))
	for j = 1:n
	   for i in nzrange(B, j)
	      if vals[i] != one(eltype(vals)) 
            push!(Is, rows[i])
            push!(Js, j) 
            push!(Vs, vals[i])
          end
	   end
	end
    sparse(Is,Js,Vs,size(B,1),size(B,2))
end
function complexoutgoingmat(::Type{Matrix{Int}}, rn::ReactionSystem, B)
    Δ = copy(B)
    for (I,b) in pairs(Δ)
        (b == 1) && (Δ[I] = 0)
    end
    Δ
end

@doc raw"""
    complexoutgoingmat(network; sparse=false, B=reactioncomplexes(rn)[2])

Given a [`ReactionSystem`](@ref) and complex incidence matrix, ``B``, return a
matrix of size num of complexes by num of reactions that identifies substrate
complexes.

Notes:
- The complex outgoing matrix, ``\Delta``, is defined by 
```math
\Delta_{i j} = \begin{cases}
    = 0,    &\text{if } B_{i j} = 1, \\
    = B_{i j}, &\text{otherwise.}
\end{cases}
```
- Set sparse=true for a sparse matrix representation
"""
function complexoutgoingmat(rn::ReactionSystem; sparse=false, B=reactioncomplexes(rn,sparse=sparse)[2])
    isempty(get_systems(rn)) || error("complexoutgoingmat does not currently support subsystems.")
    sparse ? complexoutgoingmat(SparseMatrixCSC{Int,Int}, rn, B) : 
             complexoutgoingmat(Matrix{Int}, rn, B)
end


"""
    incidencematgraph(incidencemat)   

Given an incidence matrix of a reaction-network, construct a directed simple
graph where nodes correspond to reaction complexes and directed edges to
reactions converting between two complexes.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end β ν
rcs,incidencemat = reactioncomplexes(sir)
incidencegraph   = incidencematgraph(incidencemat)
```
"""
function incidencematgraph(incidencemat::Matrix{Int})
    @assert all(∈([-1,0,1]) ,incidencemat)
    n = size(incidencemat,1)  # no. of nodes/complexes
    graph = LG.DiGraph(n)
    for col in eachcol(incidencemat)
        src = 0; dst = 0;
        for i in eachindex(col)
                (col[i] == -1) && (src = i)
                (col[i] == 1) && (dst = i)
                (src != 0) && (dst != 0) && break
        end
        LG.add_edge!(graph, src, dst)
    end
    return graph
end
function incidencematgraph(incidencemat::SparseMatrixCSC{Int,Int})
    @assert all(∈([-1,0,1]) ,incidencemat)
    m,n = size(incidencemat)  
    graph = LG.DiGraph(m)
    rows = rowvals(incidencemat)
    vals = nonzeros(incidencemat)
    for j = 1:n
        inds=nzrange(incidencemat, j)
        row = rows[inds];
        val = vals[inds];
        if val[1] == -1
            LG.add_edge!(graph, row[1], row[2])
        else
            LG.add_edge!(graph, row[2], row[1])
        end
    end
    return graph
end


"""
    linkageclasses(incidencegraph)

Given the incidence graph of a reaction network, return a vector of the
connected components of the graph (i.e. sub-groups of reaction complexes that
are connected in the incidence graph).

For example, continuing the example from [`incidencematgraph`](@ref)
```julia
julia> linkageclasses(incidencegraph)
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3, 4]
```
"""
function linkageclasses(incidencegraph)
    LG.connected_components(incidencegraph)
end


@doc raw"""
    deficiency(netstoich_mat, incidence_graph, linkage_classes)

Calculate the deficiency of a reaction network. 

Here the deficiency, ``\delta``, of a network with ``n`` reaction complexes, 
``\ell`` linkage classes and a rank ``s`` stoichiometric matrix is

```math
\delta = n - \ell - s
```

For example, 
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end β ν
rcs,incidencemat = reactioncomplexes(sir)
incidence_graph  = incidencematgraph(incidencemat)
linkage_classes   = linkageclasses(incidence_graph)
netstoich_mat    = netstoichmat(sir)
δ = deficiency(netstoich_mat, incidence_graph, linkage_classes)
```
"""
function deficiency(ns, ig, lc)
    LG.nv(ig) - length(lc) - AA.rank(AA.matrix(AA.ZZ,ns))
end

function subnetworkmapping(linkageclass, allrxs, complextorxmap, p)
    rxinds  = sort!(collect(Set(rxidx for rcidx in linkageclass for rxidx in complextorxmap[rcidx])))
    rxs     = allrxs[rxinds]
    specset = Set(substrate for rx in rxs for substrate in rx.substrates)
    for rx in rxs
        for product in rx.products
            push!(specset, product)
        end
    end
    specs = collect(specset)
    newps = Vector{eltype(p)}()
    for rx in rxs
        Symbolics.get_variables!(newps, rx.rate, p)
    end
    rxs, specs, newps   # reactions and species involved in reactions of subnetwork
end
	
"""
    subnetworks(network, linkage_classes ; rxs = reactions(network),
                  complextorxmap = collect(values(reactioncomplexmap(network))),
                  p = parameters(network))

Find subnetworks corresponding to the each linkage class of reaction network

For example, continuing the example from [`deficiency`](@ref)
```julia
   subnets = subnetworks(sir, linkage_classes)
```
"""
function subnetworks(rs::ReactionSystem, lcs::AbstractVector;
                  rxs = reactions(rs),
                  complextorxmap = [map(first,rcmap) for rcmap in values(reactioncomplexmap(rs))],
                  p = parameters(rs))
    isempty(get_systems(rs)) || error("subnetworks does not currently support subsystems.")
    t = get_iv(rs)
    subreac = Vector{ReactionSystem}()
    for i in 1:length(lcs)
        reacs,specs,newps = subnetworkmapping(lcs[i], rxs, complextorxmap, p)
        newname = Symbol(nameof(rs), "_", i)
        push!(subreac, ReactionSystem(reacs, t, specs, newps; name=newname))
    end
    subreac
end


"""
    linkagedeficiencies(subnetworks::AbstractVector, linkage_classes::AbstractVector)

Calculates the deficiency of each sub-reaction network defined by a collection
of linkage_classes.

For example, continuing the example from [`deficiency`](@ref)
```julia
subnets = subnetworks(sir, linkage_classes)
linkage_deficiencies = linkagedeficiency(subnets, linkage_classes)
```
"""
function linkagedeficiencies(subnets, lcs)
    δ = zeros(Int,length(lcs))
    for (i,subnet) in enumerate(subnets)
        ns_sub = netstoichmat(subnet)
        δ[i] = length(lcs[i]) - 1 - AA.rank(AA.matrix(AA.ZZ, ns_sub))
    end
    δ
end
	
						
"""
    isreversible(incidencegraph)

Given an incidence graph of the reaction network, returns if the network is reversible or not.
For example, continuing the example from [`linkagedeficiencies`](@ref)
```julia
isreversible(incidence_graph)
```
"""
function isreversible(ig::LG.SimpleDiGraph)
    LG.reverse(ig) == ig
end

"""
    isweaklyreversible(subnetworks)

Given the subnetworks corresponding to the each linkage class of reaction network,
determines if the reaction network is weakly reversible or not.
For example, continuing the example from [`isreversible`](@ref)
```julia
isweaklyreversible(subnets)
```
"""
function isweaklyreversible(subnets::Vector{ReactionSystem})
    igs = [incidencematgraph(reactioncomplexes(subrs)[2]) for subrs in subnets]
    all([LG.is_strongly_connected(ig) for ig in igs])
end
								
								
################################################################################################
######################## conservation laws ###############################

""" 
    conservationlaws(netstoichmat::AbstractMatrix)::Matrix

Given the net stoichiometry matrix of a reaction system, computes a matrix of
conservation laws, each represented as a row in the output. 
"""
function conservationlaws(nsm::AbstractMatrix)

    # We basically have to compute the left null space of the matrix
    # over the integers. We do this using AbstractAlgebra's integer (ZZ) interface.
    N   = AA.nullspace(AA.matrix(AA.ZZ, nsm'))[2]
    
    # to save allocations we manually take the adjoint when converting back
    # to a Julia integer matrix from the Nemo matrix. 
    ret = [convert(Int,N[i,j]) for j=1:size(N,2), i=1:size(N,1)]  

    # If all coefficients for a conservation law are negative
    # we might as well flip them to become positive
    for retrow in eachrow(ret)
        all(r -> r <= 0, retrow) && (retrow .*= -1)
    end
    
    ret
end

"""
    conservedquantities(state, cons_laws)

Compute conserved quantities for a system with the given conservation laws.
"""
conservedquantities(state, cons_laws) = cons_laws * state

######################## reaction network operators #######################

"""
    ==(rx1::Reaction, rx2::Reaction)

Tests whether two [`Reaction`](@ref)s are identical.

Notes:
- Ignores the order in which stoichiometry components are listed.
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be
    considered different than `(A+1)^2`.
"""
function (==)(rx1::Reaction, rx2::Reaction)
    isequal(rx1.rate, rx2.rate) || return false
    issetequal(zip(rx1.substrates,rx1.substoich), zip(rx2.substrates,rx2.substoich)) || return false
    issetequal(zip(rx1.products,rx1.prodstoich), zip(rx2.products,rx2.prodstoich)) || return false
    issetequal(rx1.netstoich, rx2.netstoich) || return false
    rx1.only_use_rate == rx2.only_use_rate
end

function hash(rx::Reaction, h::UInt) 
    h = Base.hash(rx.rate, h)
    h = Base.hash(rx.substrates, h)
    h = Base.hash(rx.products, h)
    h = Base.hash(rx.prodstoich, h)
    h = Base.hash(rx.substoich, h)
    h = Base.hash(rx.netstoich, h)
    Base.hash(rx.only_use_rate, h)
end

"""
    isequal_ignore_names(rn1::ReactionSystem, rn2::ReactionSystem)

Tests whether the underlying species, parameters and reactions are the same in
the two [`ReactionSystem`](@ref)s. Ignores the names of the systems in testing
equality.

Notes:
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be
    considered different than `(A+1)^2`.
- Does not include `defaults` in determining equality.
"""
function isequal_ignore_names(rn1::ReactionSystem, rn2::ReactionSystem)
    isequal(get_iv(rn1), get_iv(rn2)) || return false
    issetequal(get_states(rn1), get_states(rn2)) || return false
    issetequal(get_ps(rn1), get_ps(rn2)) || return false
    issetequal(MT.get_observed(rn1), MT.get_observed(rn2)) || return false
    issetequal(get_eqs(rn1), get_eqs(rn2)) || return false
    (get_constraints(rn1) == get_constraints(rn2)) || return false

    # subsystems
    (length(get_systems(rn1)) == length(get_systems(rn2))) || return false
    issetequal(get_systems(rn1), get_systems(rn2)) || return false

    true
end

"""
    ==(rn1::ReactionSystem, rn2::ReactionSystem)

Tests whether the underlying species, parameters and reactions are the same in
the two [`ReactionSystem`](@ref)s. Requires the systems to have the same names
too.

Notes:
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be
    considered different than `(A+1)^2`.
- Does not include `defaults` in determining equality.
"""
function (==)(rn1::ReactionSystem, rn2::ReactionSystem)
    (nameof(rn1) == nameof(rn2)) || return false
    isequal_ignore_names(rn1,rn2)
end


######################## functions to extend a network ####################

"""
    make_empty_network(; iv=DEFAULT_IV, name=gensym(:ReactionSystem))

Construct an empty [`ReactionSystem`](@ref). `iv` is the independent variable,
usually time, and `name` is the name to give the `ReactionSystem`.
"""
function make_empty_network(; iv=DEFAULT_IV, name=gensym(:ReactionSystem))
    ReactionSystem(Reaction[], iv, [], []; name=name)
end

"""
    addspecies!(network::ReactionSystem, s::Symbolic; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the species corresponding to the variable
`s` to the network (if it is not already defined). Returns the integer id of the
species within the system.

Notes:
- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addspecies!(network::ReactionSystem, s::Symbolic; disablechecks=false)

    # we don't check subsystems since we will add it to the top-level system...
    curidx = disablechecks ? nothing : findfirst(S -> isequal(S, s), get_states(network))
    if curidx === nothing
        push!(get_states(network), s)
        return length(get_states(network))
    else
        return curidx
    end
end


"""
    addspecies!(network::ReactionSystem, s::Num; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the species corresponding to the
variable `s` to the network (if it is not already defined). Returns the
integer id of the species within the system.

- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addspecies!(network::ReactionSystem, s::Num; disablechecks=false)
    addspecies!(network, value(s), disablechecks=disablechecks)
end

"""
    reorder_states!(rn, neworder)

Given a [`ReactionSystem`](@ref) and a vector `neworder`, orders the states of `rn` accordingly to `neworder`.

Notes:
- Currently only supports `ReactionSystem`s without constraints or subsystems.
"""
function reorder_states!(rn, neworder) 
   (get_constraints(rn) === nothing) && isempty(get_systems(rn)) || 
           error("Reordering of states is only supported for systems without constraints or subsystems.")
    permute!(get_states(rn), neworder)
end

"""
    addparam!(network::ReactionSystem, p::Symbolic; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the parameter corresponding to the
variable `p` to the network (if it is not already defined). Returns the integer
id of the parameter within the system.

- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addparam!(network::ReactionSystem, p::Symbolic; disablechecks=false)
    # we don't check subsystems since we will add it to the top-level system...
    if istree(p) && !(operation(p) isa Sym)
        error("If the passed in parameter is an expression, it must correspond to an underlying Variable.")
    end
    curidx = disablechecks ? nothing : findfirst(S -> isequal(S, p), get_ps(network))
    if curidx === nothing
        push!(get_ps(network), p)
        return length(get_ps(network))
    else
        return curidx
    end
end

"""
    addparam!(network::ReactionSystem, p::Num; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the parameter corresponding to the
variable `p` to the network (if it is not already defined). Returns the
integer id of the parameter within the system.

- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addparam!(network::ReactionSystem, p::Num; disablechecks=false)
    addparam!(network, value(p); disablechecks=disablechecks)
end

"""
    addreaction!(network::ReactionSystem, rx::Reaction)

Add the passed in reaction to the [`ReactionSystem`](@ref). Returns the
integer id of `rx` in the list of `Reaction`s within `network`.

Notes:
- Any new species or parameters used in `rx` should be separately added to
    `network` using [`addspecies!`](@ref) and [`addparam!`](@ref).
"""
function addreaction!(network::ReactionSystem, rx::Reaction)
    push!(get_eqs(network), rx)
    length(get_eqs(network))
end

"""
    merge!(network1::ReactionSystem, network2::ReactionSystem)

Merge `network2` into `network1`.

Notes:
- Duplicate reactions between the two networks are not filtered out.
- [`Reaction`](@ref)s are not deepcopied to minimize allocations, so both
  networks will share underlying data arrays.
- Subsystems are not deepcopied between the two networks and will hence be
  shared.
- Returns `network1`.
"""
function Base.merge!(network1::ReactionSystem, network2::ReactionSystem)
    ((get_constraints(network1) === nothing) && (get_constraints(network2) === nothing)) ||
        error("merge! does not currently support ReactionSystems with constraints, consider ModelingToolkit.extend instead.")
    isequal(get_iv(network1), get_iv(network2)) || 
        error("Reaction networks must have the same independent variable to be mergable.")
    append!(get_eqs(network1), get_eqs(network2))
    union!(get_states(network1), get_states(network2))
    union!(get_ps(network1), get_ps(network2))
    append!(get_observed(network1), get_observed(network2))    
    append!(get_systems(network1), get_systems(network2))
    merge!(get_defaults(network1), get_defaults(network2))
    network1
end

bracket(s::Num) = bracket(MT.value(s))

function bracket(s)     
    if MT.istree(s)
        op = MT.operation(s)
        args = MT.arguments(s)
        (length(args) == 1) || error("Unable to indentify a unique scalar independent variable.")
        ssym = Symbol("⟦", Symbol(op), "⟧")
        return (@variables $ssym(t))[1]
    else
        ssym = Symbol("⟦", Symbol(s), "⟧")
        return (@variables $ssym)[1] 
    end
end

function addconcentrations(rn, vol, specs=species(rn))
    obs = get_observed(rn)
    t   = get_iv(rn)
    for s in specs        
        push!(obs, Equation(bracket(s), s/vol))
    end
    nothing
end

###############################   units   #####################################

"""
    validate(rx::Reaction; info::String = "")     

Check that all substrates and products within the given [`Reaction`](@ref) have
the same units, and that the units of the reaction's rate expression are
internally consistent (i.e. if the rate involves sums, each term in the sum has
the same units).

"""
function validate(rx::Reaction; info::String = "")     
    validated = MT._validate([rx.rate], [string(rx, ": rate")], info = info)
    
    subunits = isempty(rx.substrates) ? nothing : get_unit(rx.substrates[1])
    for i in 2:length(rx.substrates)
        if get_unit(rx.substrates[i]) != subunits
            validated = false
            @warn(string("In ", rx, " the substrates have differing units."))
        end
    end

    produnits = isempty(rx.products) ? nothing : get_unit(rx.products[1])
    for i in 2:length(rx.products)
        if get_unit(rx.products[i]) != produnits
            validated = false
            @warn(string("In ", rx, " the products have differing units."))
        end
    end

    if (subunits !== nothing) && (produnits !== nothing) && (subunits != produnits)
        validated = false
        @warn(string("in ", rx, " the substrate units are not consistent with the product units."))
    end

    validated
end

"""
    validate(rs::ReactionSystem, info::String="")

Check that all species in the [`ReactionSystem`](@ref) have the same units, and
that the rate laws of all reactions reduce to units of (species units) / (time
units).

Notes:
- Does not check subsystems too. 
"""
function validate(rs::ReactionSystem, info::String="")
    specs = get_states(rs)

    # if there are no species we don't check units on the system
    isempty(specs) && return true   

    specunits = get_unit(specs[1])
    validated = true
    for spec in specs
        if get_unit(spec) != specunits
            validated = false 
            @warn(string("Species are expected to have units of ", specunits, " however, species ", spec, " has units ", get_unit(spec), "."))
        end
    end
    timeunits = get_unit(get_iv(rs))
    rateunits = specunits / timeunits

    for rx in get_eqs(rs)
        rxunits = get_unit(rx.rate)
        for (i,sub) in enumerate(rx.substrates)
            rxunits *= get_unit(sub) ^ rx.substoich[i]
        end

        if rxunits != rateunits
            validated = false
            @warn(string("Reaction rate laws are expected to have units of ", rateunits, " however, ", rx, " has units of ", rxunits, "."))
        end
    end

    validated
end
