struct MatrixNetwork{T <: AbstractMatrix, W <: AbstractVector{Int}} 
    ns::T
    ss::T
    ridxs::W
    rn::DiffEqBase.AbstractReactionNetwork
end

function MatrixNetwork{T}(rn::DiffEqBase.AbstractReactionNetwork) where T <: SparseMatrixCSC
    ss = DiffEqBiological.substratestoich_sparsemat(rn)
    ns = DiffEqBiological.netstoich_sparsemat(rn, substoich=ss)    
    rateidxs = map(i -> paramsmap(rs)[rateexpr(rs,i)], 1:numreactions(rs))
    MatrixNetwork(ns, ss, rateidxs, rn)
end

function oderhs!(du, u, p, t, mn::MatrixNetwork{T,W}) where {T <: SparseMatrixCSC,W}
    @unpack ns, ss, ridxs, rn = mn

    fill!(du, zero(eltype(du)))
    numrx = numreactions(rn)
    srows = rowvals(ss)
    svals = nonzeros(ss)
    nrows = rowvals(ns)
    nvals = nonzeros(ns)
    for j = 1:numrx

        # calculate rate law from substrate stoich
        ratelaw = p[ridxs[j]]
        coef = 1
        for ir in nzrange(ss, j)
            i = srows[ir]
            stoich = svals[ir]
            if stoich == 1
                ratelaw *= u[i] 
            else
                ratelaw *= u[i]^stoich
                coef *= factorial(stoich)
            end
        end

        # add in to du using net stoich

end