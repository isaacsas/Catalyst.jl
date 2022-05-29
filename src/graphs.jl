#######################################################################
# Contains code from  Catlab.jl:
# https://raw.githubusercontent.com/AlgebraicJulia/Catlab.jl/master/src/graphics/Graphviz.jl
#
# That license for that code is:
#
# The MIT License
#
# Copyright (c) 2017-2020: Evan Patterson.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#######################################################################
""" AST and pretty printer for Graphviz's DOT language.

References:

- DOT grammar: http://www.graphviz.org/doc/info/lang.html
- DOT language guide: http://www.graphviz.org/pdf/dotguide.pdf
"""

# AST
#####

abstract type Expression end
abstract type Statement <: Expression end

""" AST type for Graphviz's "HTML-like" node labels.

For now, the HTML content is just a string.
"""
struct Html
  content::String
end
Base.print(io::IO, html::Html) = print(io, html.content)

const AttributeValue = Union{String,Html}
const Attributes = OrderedDict{Symbol,AttributeValue}

as_attributes(attrs::Attributes) = attrs
as_attributes(d::OrderedDict) = Attributes(Symbol(k) => d[k] for k in keys(d))
as_attributes(d::AbstractDict) =
  Attributes(Symbol(k) => d[k] for k in sort!(collect(keys(d))))

@with_kw_noshow struct Graph <: Expression
  name::String
  directed::Bool
  prog::String = "dot"
  stmts::Vector{Statement} = Statement[]
  graph_attrs::Attributes = Attributes()
  node_attrs::Attributes = Attributes()
  edge_attrs::Attributes = Attributes()
end

Graph(name::String, stmts::Vector{Statement}; kw...) =
  Graph(; name=name, directed=false, stmts=stmts, kw...)
Graph(name::String, stmts::Vararg{Statement}; kw...) =
  Graph(; name=name, directed=false, stmts=collect(stmts), kw...)
Digraph(name::String, stmts::Vector{Statement}; kw...) =
  Graph(; name=name, directed=true, stmts=stmts, kw...)
Digraph(name::String, stmts::Vararg{Statement}; kw...) =
  Graph(; name=name, directed=true, stmts=collect(stmts), kw...)

@with_kw_noshow struct Subgraph <: Statement
  name::String = "" # Subgraphs can be anonymous
  stmts::Vector{Statement} = Statement[]
  graph_attrs::Attributes = Attributes()
  node_attrs::Attributes = Attributes()
  edge_attrs::Attributes = Attributes()
end

Subgraph(stmts::Vector{Statement}; kw...) = Subgraph(; stmts=stmts, kw...)
Subgraph(stmts::Vararg{Statement}; kw...) = Subgraph(; stmts=collect(stmts), kw...)
Subgraph(name::String, stmts::Vector{Statement}; kw...) =
  Subgraph(; name=name, stmts=stmts, kw...)
Subgraph(name::String, stmts::Vararg{Statement}; kw...) =
  Subgraph(; name=name, stmts=collect(stmts), kw...)

struct Node <: Statement
  name::String
  attrs::Attributes
end
Node(name::String, attrs::AbstractDict) = Node(name, as_attributes(attrs))
Node(name::String; attrs...) = Node(name, attrs)

struct NodeID <: Expression
  name::String
  port::String
  anchor::String
  NodeID(name::String, port::String="", anchor::String="") = new(name, port, anchor)
end

struct Edge <: Statement
  path::Vector{NodeID}
  attrs::Attributes
end
Edge(path::Vector{NodeID}, attrs::AbstractDict) = Edge(path, as_attributes(attrs))
Edge(path::Vector{NodeID}; attrs...) = Edge(path, attrs)
Edge(path::Vararg{NodeID}; attrs...) = Edge(collect(path), attrs)
Edge(path::Vector{String}, attrs::AbstractDict) = Edge(map(NodeID, path), attrs)
Edge(path::Vector{String}; attrs...) = Edge(map(NodeID, path), attrs)
Edge(path::Vararg{String}; attrs...) = Edge(map(NodeID, collect(path)), attrs)

# Bindings
##########

""" Run a Graphviz program.

Invokes Graphviz through its command-line interface. If the `Graphviz_jll`
package is installed and loaded, it is used; otherwise, Graphviz must be
installed on the local system.

For bindings to the Graphviz C API, see the the package
[GraphViz.jl](https://github.com/Keno/GraphViz.jl). At the time of this writing,
GraphViz.jl is unmaintained.
"""
function run_graphviz(io::IO, graph::Graph; prog::Union{String,Nothing}=nothing,
  format::String="json0")
  if isnothing(prog)
    prog = graph.prog
  end
  @assert prog in ("dot", "neato", "fdp", "sfdp", "twopi", "circo")
#   if USE_GV_JLL[]
    fun = getfield(Graphviz_jll, Symbol(prog))
    fun() do exe
      open(`$exe -T$format`, io, write=true) do gv
        pprint(gv, graph)
      end
    end
#   else
#     open(`$prog -T$format`, io, write=true) do gv
#       pprint(gv, graph)
#     end
#   end
end
function run_graphviz(graph::Graph; kw...)
  io = IOBuffer()
  run_graphviz(io, graph; kw...)
  seekstart(io)
end

function Base.show(io::IO, ::MIME"image/svg+xml", graph::Graph)
  run_graphviz(io, graph, format="svg")
end

# Pretty-print
##############

""" Pretty-print the Graphviz expression.
"""
pprint(expr::Expression) = pprint(stdout, expr)
pprint(io::IO, expr::Expression) = pprint(io, expr, 0)

function pprint(io::IO, graph::Graph, n::Int)
  indent(io, n)
  print(io, graph.directed ? "digraph " : "graph ")
  print(io, graph.name)
  println(io, " {")
  pprint_attrs(io, graph.graph_attrs, n + 2; pre="graph", post=";\n")
  pprint_attrs(io, graph.node_attrs, n + 2; pre="node", post=";\n")
  pprint_attrs(io, graph.edge_attrs, n + 2; pre="edge", post=";\n")
  for stmt in graph.stmts
    pprint(io, stmt, n + 2, directed=graph.directed)
    println(io)
  end
  indent(io, n)
  println(io, "}")
end

function pprint(io::IO, subgraph::Subgraph, n::Int; directed::Bool=false)
  indent(io, n)
  if isempty(subgraph.name)
    println(io, "{")
        else
    print(io, "subgraph ")
    print(io, subgraph.name)
    println(io, " {")
  end
  pprint_attrs(io, subgraph.graph_attrs, n + 2; pre="graph", post=";\n")
  pprint_attrs(io, subgraph.node_attrs, n + 2; pre="node", post=";\n")
  pprint_attrs(io, subgraph.edge_attrs, n + 2; pre="edge", post=";\n")
  for stmt in subgraph.stmts
    pprint(io, stmt, n + 2, directed=directed)
    println(io)
  end
  indent(io, n)
  print(io, "}")
end

function pprint(io::IO, node::Node, n::Int; directed::Bool=false)
  indent(io, n)
  print(io, node.name)
  pprint_attrs(io, node.attrs)
  print(io, ";")
end

function pprint(io::IO, node::NodeID, n::Int)
  print(io, node.name)
  if !isempty(node.port)
    print(io, ":")
    print(io, node.port)
  end
  if !isempty(node.anchor)
    print(io, ":")
    print(io, node.anchor)
end
end

function pprint(io::IO, edge::Edge, n::Int; directed::Bool=false)
  indent(io, n)
  for (i, node) in enumerate(edge.path)
    if i > 1
      print(io, directed ? " -> " : " -- ")
    end
    pprint(io, node, n)
  end
  pprint_attrs(io, edge.attrs)
  print(io, ";")
end

function pprint_attrs(io::IO, attrs::Attributes, n::Int=0;
                      pre::String="", post::String="")
  if !isempty(attrs)
    indent(io, n)
        print(io, pre)
    print(io, " [")
    for (i, (key, value)) in enumerate(attrs)
      if (i > 1) print(io, ",") end
      print(io, key)
      print(io, "=")
      print(io, value isa Html ? "<" : "\"")
      print(io, value)
      print(io, value isa Html ? ">" : "\"")
    end
    print(io, "]")
    print(io, post)
  end
end

indent(io::IO, n::Int) = print(io, " "^n)


#######################################################################
# following is adapted from Petri.jl
# https://github.com/mehalter/Petri.jl
#######################################################################

const graph_attrs = Attributes(:rankdir => "LR")
const node_attrs  = Attributes(:shape => "plain", :style => "filled", :color => "white")
const edge_attrs  = Attributes(:splines => "splines")

function edgify(δ, i, reverse::Bool)
    attr = Attributes()
    return map(δ) do p
        val = String(p[1].f.name)
      weight = "$(p[2])"
      attr = Attributes(:label => weight, :labelfontsize => "6")
      return Edge(reverse ? ["rx_$i", "$val"] :
                            ["$val", "rx_$i"], attr)
    end
end

# make distinguished edge based on rate constant
function edgifyrates(rxs, specs)
    es = Edge[]
    for (i, rx) in enumerate(rxs)
        deps = rx.rate isa Number ? Any[] : get_variables(rx.rate, specs)
        for dep in deps
            val = String(dep.f.name)
            attr = Attributes(:color => "#d91111", :style => "dashed")
            e = Edge(["$val", "rx_$i"], attr)
            push!(es, e)
        end
    end
    es
end

# create edges from one reaction complex to another reaction complex
function edgifycomplex(δ,attr)
    return map(δ) do p
        return Edge([p[1], p[2]] , attr)
    end
end

# modify vector of string of complexes into graphviz compatible strings
function modifystrcomp(strcomp::Vector{String})
  for i in 1:length(strcomp)
    #   cannot allow (t) as per graphviz
      if occursin("(t)",strcomp[i])
          strcomp[i] = replace(strcomp[i], "(t)" => "")
      end
      #   prettify NUll with ∅
      if occursin("0",strcomp[i])
          strcomp[i] = replace(strcomp[i], "0" => "∅")
      end
  end
  strcomp = "<".*strcomp.*">"
end

"""
    complexgraph(rn::ReactionSystem; complexdata=reactioncomplexes(rn))

Creates a Graphviz graph of the [`ReactionComplex`](@ref)s in `rn`. Reactions
correspond to arrows and reaction complexes to blue circles.

Notes:
- Black arrows from complexes to complexes indicate reactions whose rate is a
  parameter or a `Number`. i.e. `k, A --> B`.
- Red dashed arrows from complexes to complexes indicate reactions whose rate
  depends on species. i.e. `k*C, A --> B` for `C` a species.
- Requires the Graphviz jll to be installed, or Graphviz to be installed and
  commandline accessible.
"""
function complexgraph(rn::ReactionSystem; complexdata=reactioncomplexes(rn))
    rxs   = reactions(rn);
    specs = species(rn);
    complexes,B = complexdata;
    fun = rcel -> specs[rcel.speciesid]*rcel.speciesstoich;
    compfun(rc) = rc == Catalyst.ReactionComplex{Int64}[] ? 0 : sum(fun, rc);

    strcomp = [string(compfun(rc)) for rc in complexes];

    newstrcomp = modifystrcomp(strcomp)
    compnodes = [Node(str, Attributes(:shape => "circle",:color => "#6C9AC3")) for str in newstrcomp]

    edges = []
    for (i,r) in enumerate(rxs)
        subcomp = newstrcomp[argmin(@view B[:,i])]
        prodcomp = newstrcomp[argmax(@view B[:,i])]
        deps = get_variables(r.rate, specs)
        if deps != Any[]
            attr = Attributes(:color => "#d91111", :style => "dashed")
            push!(edges, edgifycomplex(zip([subcomp],[prodcomp]),attr))
        else
            attr = Attributes()
            push!(edges,edgifycomplex(zip([subcomp],[prodcomp]),attr))
        end
    end
    stmts2 = Vector{Statement}()
    append!(stmts2, compnodes)
    append!(stmts2, collect(Iterators.flatten(edges)))
    g = Digraph("G", stmts2; graph_attrs=graph_attrs, node_attrs=node_attrs,edge_attrs=edge_attrs)
    return g
end

"""
    Graph(rn::ReactionSystem)

Converts a [`ReactionSystem`](@ref) into a Graphviz graph.
Reactions correspond to small green circles, and species to blue circles.

Notes:
- Black arrows from species to reactions indicate reactants, and are labelled
  with their input stoichiometry.
- Black arrows from reactions to species indicate products, and are labelled
  with their output stoichiometry.
- Red arrows from species to reactions indicate that species is used within the
  rate expression. For example, in the reaction `k*A, B --> C`, there would be a
  red arrow from `A` to the reaction node. In `k*A, A+B --> C`, there would be
  red and black arrows from `A` to the reaction node.
- Requires the Graphviz jll to be installed, or Graphviz to be installed and
  commandline accessible.
"""
function Graph(rn::ReactionSystem)
    rxs   = reactions(rn)
    specs = species(rn)
    statenodes = [Node(string(s.f.name), Attributes(:shape => "circle", :color => "#6C9AC3")) for s in specs]
    transnodes = [Node(string("rx_$i"), Attributes(:shape => "point", :color => "#E28F41", :width => ".1")) for (i, r) in enumerate(rxs)]

    stmts = vcat(statenodes, transnodes)
    edges = map(enumerate(rxs)) do (i, r)
      vcat(edgify(zip(r.substrates, r.substoich), i, false),
           edgify(zip(r.products, r.prodstoich), i, true))
    end
    es = edgifyrates(rxs, specs)
    (!isempty(es)) && push!(edges, es)

    stmts2 = Vector{Statement}()
    append!(stmts2, stmts)
    append!(stmts2, collect(Iterators.flatten(edges)))
    g = Digraph("G", stmts2; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
    return g
end


"""
    savegraph(g::Graph, fname, fmt="png")

Given a `Graph` generated by [`Graph`](@ref), save the graph to the file with
name `fname` and extension `fmt`.

Notes:
- `fmt="png"` is the default output format.
- Requires the Graphviz jll to be installed, or Graphviz to be installed and commandline accessible.
"""
function savegraph(g::Graph, fname, fmt="png")
    open(fname, "w") do io
        run_graphviz(io, g, format=fmt)
    end
    nothing
end
