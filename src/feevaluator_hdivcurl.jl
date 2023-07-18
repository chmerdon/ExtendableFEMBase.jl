# IDENTITY Hdivcurl
function update_basis!(FEBE::SingleFEEvaluator{<:Real,<:Real,<:Integer,<:Identity,<:AbstractHdivcurlFiniteElement})
    L2GM = _update_piola!(FEBE)
    L2GAinv = _update_trafo!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    refbasisvals = FEBE.refbasisvals
    cvals = FEBE.cvals
    fill!(cvals,0)
    for i = 1 : size(cvals,3), dof_i = 1 : size(cvals,2), k = 1 : size(cvals,1)
        # todo: L2GAinv * refbasisvals * L2GM
        #for l = 1 : size(L2GAinv,2)
        #    cvals[k,dof_i,i] += L2GAinv[k,l] * refbasisvals[i][subset[dof_i],l]
        #end    
        cvals[k,dof_i,i] *= coefficients[k,dof_i]
    end
    return nothing
end