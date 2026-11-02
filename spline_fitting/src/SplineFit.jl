module SplineFit
using LinearAlgebra, MixedModels, StatsModels
export build_model, fit_reml, restore_observations, predict_response
include("SharedVariance.jl")

"""One common observation rotation, with coefficient rotations restricted to each component.

For concatenated Z=Q*[R;0], let R=U*D*V'. Use T=blockdiag(U',I)*Q'.
For each component choose its own orthogonal V_j, and retain Zr_j=T*Z_j*V_j.
For one component Zr is diagonal; for several components it is generally dense.
Source: MixedModels 5.9.0 src/remat.jl and src/linearmixedmodel.jl.
"""
function rotate_inputs(yf,Xf,Z_blocks)
    all_Z=hcat(Z_blocks...)
    q=size(all_Z,2)
    qrZ=qr(all_Z)
    R=Matrix(qrZ.R)
    complete=svd(R)
    minimum(complete.S)>1e-10*maximum(complete.S) || throw(ArgumentError("Combined penalized columns must be linearly independent"))
    rotated=qrZ.Q'*hcat(yf,Xf)
    rotated[1:q,:]=complete.U'*rotated[1:q,:]
    V_blocks=Matrix{Float64}[];Zr=Matrix{Float64}[];singular_values=Vector{Float64}[]
    offset=0
    for Z in Z_blocks
        r=size(Z,2);cols=(offset+1):(offset+r)
        parts=svd(R[:,cols])
        top=complete.U'*R[:,cols]*parts.V
        push!(Zr,vcat(top,zeros(length(yf)-q,r)))
        push!(V_blocks,parts.V);push!(singular_values,parts.S)
        offset+=r
    end
    return (;yr=rotated[:,1],Xr=rotated[:,2:end],Zr,qrZ,U=complete.U,V_blocks,
             singular_values,rotated_rows=q)
end

"""Native dense ReMat for each smooth, with one tied diagonal variance per component."""
function build_model_terms(yr,Xr,Zr,names)
    n,p=size(Xr)
    fixed_names=Symbol.("fixed_".*string.(1:p))
    entries=Pair{Symbol,Any}[:response=>yr]
    append!(entries,[key=>Xr[:,j] for (j,key) in enumerate(fixed_names)])
    fixed=foldl(+,term.(fixed_names);init=term(0))
    random_terms=Any[]
    for (j,Z) in enumerate(Zr)
        group=Symbol("smooth_"*lpad(string(j),length(string(length(Zr))),'0'))
        push!(entries,group=>fill("spline",n))
        columns=Symbol.("z"*string(j)*"_".*string.(1:size(Z,2)))
        append!(entries,[key=>Z[:,k] for (k,key) in enumerate(columns)])
        lhs=foldl(+,term.(columns);init=term(0))
        push!(random_terms,lhs|term(group))
    end
    tbl=NamedTuple{Tuple(first.(entries))}(Tuple(last.(entries)))
    formula=FormulaTerm(term(:response),fixed+sum(random_terms))
    form=MixedModels.schematize(formula,tbl,Dict{Symbol,Any}())
    _,design=StatsModels.modelcols(form,tbl)
    reterms=[SharedVarianceTerm(t) for t in design if t isa MixedModels.ReMat]
    feterm=MixedModels.FeTerm(Xr,String.(fixed_names))
    group_names=[MixedModels.fname(t) for t in reterms]
    return (;feterm,reterms,form,group_names)
end

"""Accept a variable-length list of spline designs, one shared variance per entry."""
function build_model(y::AbstractVector,X::AbstractMatrix,Z_blocks::AbstractVector{<:AbstractMatrix};names=["x"*string(j) for j in eachindex(Z_blocks)])
    yf,Xf=Vector{Float64}(y),Matrix{Float64}(X)
    Zs=Matrix{Float64}.(Z_blocks)
    n,p=size(Xf)
    !isempty(Zs) && length(names)==length(Zs) && length(unique(names))==length(names) || throw(ArgumentError("Supply distinct names and at least one component"))
    length(yf)==n && all(size(Z,1)==n for Z in Zs) || throw(DimensionMismatch("All inputs must share observation rows"))
    p>0 && all(size(Z,2)>0 for Z in Zs) && n>p+sum(size(Z,2) for Z in Zs) || throw(ArgumentError("Require fixed and penalized columns and n > p + total penalized columns"))
    all(isfinite,yf) && all(isfinite,Xf) && all(all(isfinite,Z) for Z in Zs) || throw(ArgumentError("Model inputs must be finite"))
    rank(Xf)==p || throw(ArgumentError("Fixed-effect columns must be linearly independent"))
    rotated=rotate_inputs(yf,Xf,Zs)
    terms=build_model_terms(rotated.yr,rotated.Xr,rotated.Zr,names)
    model=MixedModels.LinearMixedModel(rotated.yr,terms.feterm,terms.reterms,terms.form,Float64[],nothing,false,1000,false)
    length(model.θ)==length(Zs) || error("Expected one variance parameter per component")
    # The constructor sorts terms by size. Preserve the user-facing component order.
    group_names=MixedModels.fname.(model.reterms)
    model_order=[only(findall(==(key),group_names)) for key in terms.group_names]
    for j in eachindex(Zs)
        isapprox(Matrix(model.reterms[model_order[j]]),rotated.Zr[j];rtol=1e-12,atol=1e-12) || error("Stored random design changed")
    end
    return (;model,names=String.(names),y=yf,X=Xf,yr=rotated.yr,Xr=rotated.Xr,Z_blocks=Zs,Zr_blocks=rotated.Zr,
             qrZ=rotated.qrZ,U=rotated.U,V_blocks=rotated.V_blocks,
             singular_values=rotated.singular_values,rotated_rows=rotated.rotated_rows,model_order)
end
build_model(y::AbstractVector,X::AbstractMatrix,Z::AbstractMatrix;kwargs...) = build_model(y,X,[Z];kwargs...)

"""Undo the common observation rotation once, regardless of component count."""
function restore_observations(prepared,values::AbstractVecOrMat)
    out=copy(values);q=prepared.rotated_rows
    if out isa AbstractVector
        out[1:q]=prepared.U*out[1:q]
    else
        out[1:q,:]=prepared.U*out[1:q,:]
    end
    return prepared.qrZ.Q*out
end

"""Joint REML fitting through MixedModels; each component has its own r² and lambda."""
function fit_reml(y::AbstractVector,X::AbstractMatrix,Z_blocks::AbstractVector{<:AbstractMatrix};kwargs...)
    prepared=build_model(y,X,Z_blocks;kwargs...);m=prepared.model
    # Multi-parameter reference checks need tighter objective stopping tolerances.
    m.optsum.ftol_rel=1e-14
    m.optsum.ftol_abs=1e-12
    m.optsum.xtol_abs=fill(1e-10,length(m.θ))
    MixedModels.fit!(m;REML=true,progress=false)
    status=m.optsum.returnvalue
    status in (:SUCCESS,:STOPVAL_REACHED,:FTOL_REACHED,:XTOL_REACHED) || error("MixedModels did not converge: $status")
    theta=[abs(m.reterms[i].λ[1,1]) for i in prepared.model_order]
    sigma2=MixedModels.varest(m)
    r2=sigma2.*theta.^2;r=sqrt.(r2)
    lambda=[iszero(t) ? Inf : inv(t^2) for t in theta]
    b=MixedModels.coef(m) # Combined fixed coefficients: b=[β; α₁; ...; αₘ].
    random=MixedModels.ranef(m)
    tau=[prepared.V_blocks[j]*vec(random[prepared.model_order[j]]) for j in eachindex(Z_blocks)]
    fitted=restore_observations(prepared,MixedModels.fitted(m))
    reconstructed=predict_response((;b,tau),prepared.X,prepared.Z_blocks)
    isapprox(fitted,reconstructed;rtol=1e-10,atol=1e-10) || error("Fitted-response reconstruction failed")
    isfinite(sigma2) && sigma2>0 && all(isfinite,fitted) || error("Invalid fitted model")
    return (;prepared,lambda,sigma2,r2,r,theta,b,tau,fitted,status,
             objective=MixedModels.objective(m),evaluations=m.optsum.feval)
end
fit_reml(y::AbstractVector,X::AbstractMatrix,Z::AbstractMatrix;kwargs...) = fit_reml(y,X,[Z];kwargs...)

function predict_response(fit,X,Z_blocks)
    length(Z_blocks)==length(fit.tau) || throw(DimensionMismatch("Component counts disagree"))
    return X*fit.b+sum(Z*coef for (Z,coef) in zip(Z_blocks,fit.tau))
end
end
