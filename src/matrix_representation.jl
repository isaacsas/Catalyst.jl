struct MatrixNetwork{T <: AbstractMatrix, W <: AbstractVector{Int}} 
    ns::T
    ss::T
    ridxs::W
    rn::DiffEqBase.AbstractReactionNetwork
end

function MatrixNetwork{T}(rn::DiffEqBase.AbstractReactionNetwork) where T <: SparseMatrixCSC
    ss = DiffEqBiological.substratestoich_sparsemat(rn)
    ns = DiffEqBiological.netstoich_sparsemat(rn, substoich=ss)    
    rateidxs = map(i -> paramsmap(rn)[rateexpr(rn,i).args[end]], 1:numreactions(rn))
    MatrixNetwork(ns, ss, rateidxs, rn)
end

struct MatrixNetworkODEs{U,V,W <: MatrixNetwork{U,V}}
    mn::W
    syms
end


function (mnodes::MatrixNetworkODEs{U,V,W})(du, u, p, t) where {U<:SparseMatrixCSC,V,W <: MatrixNetwork{U,V}}
    mn = mnodes.mn
    @unpack ns, ss, ridxs = mn

    fill!(du, zero(eltype(du)))
    numrx = size(ns)[2]
    srows = rowvals(ss)
    svals = nonzeros(ss)
    nrows = rowvals(ns)
    nvals = nonzeros(ns)
    @inbounds for j = 1:numrx

        # calculate rate law from substrate stoich
        ratelaw = p[ridxs[j]]
        coef = one(eltype(svals))
        @inbounds for ir in nzrange(ss, j)
            stoich = svals[ir]
            ratelaw *= isone(stoich) ? u[srows[ir]] : u[srows[ir]]^stoich
            # if isone(stoich)
            #     ratelaw *= u[srows[ir]] 
            # else
            #     ratelaw *= u[srows[ir]]^stoich
            #     #coef *= factorial(stoich)
            # end            
        end
        #!isone(coef) && (ratelaw /= coef)

        # add in to du using net stoich
        @inbounds for ir in nzrange(ns, j)
            du[nrows[ir]] += nvals[ir]*ratelaw
        end
    end

    nothing

end

# function oderhs!(du, u, p, t, mn::MatrixNetwork{T,W}) where 
# end