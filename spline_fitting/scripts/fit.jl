using DelimitedFiles, LinearAlgebra, MixedModels, StatsModels
include(joinpath(@__DIR__,"..","..","spline_basis","src","BSpline.jl"))
using .BSplineBasis
include(joinpath(@__DIR__,"..","..","model_preparation","src","ModelInputs.jl"))
using .ModelInputs
include(joinpath(@__DIR__,"..","src","SplineFit.jl"))
using .SplineFit

function main(args)
    length(args)==2 || error("Expected shared input directory and new Julia output directory")
    inputs,out=abspath.(args)
    isdir(out)&&!isempty(readdir(out))&&error("Julia output directory is not empty")
    cases=readdlm(joinpath(inputs,"cases.tsv"),'\t',String;skipstart=1)
    mkpath(out)
    rows=NamedTuple[]
    for i in axes(cases,1)
        name,k,mode=cases[i,1],parse(Int,cases[i,2]),cases[i,3]
        src=joinpath(inputs,name); dst=joinpath(out,name); mkpath(dst)
        dat=readdlm(joinpath(src,"data.csv"),',',Float64;skipstart=1)
        x,y,truth=dat[:,1],dat[:,2],dat[:,3]
        # Independent Julia construction: no R-generated matrices or estimates are read.
        knots=mode=="default" ? default_knots(x,k) : vec(readdlm(joinpath(src,"knots.csv"),',',Float64))
        basis=build_spline(x,knots)
        prepared=prepare_inputs(basis.B,basis.S;expected_rank=k-2)
        X_fixed=hcat(ones(length(x)),prepared.X0)
        fit=fit_reml(y,X_fixed,prepared.Z)
        open(joinpath(dst,"fitted.csv"),"w") do io
            println(io,"x,y,truth,fitted")
            writedlm(io,hcat(x,y,truth,fit.fitted),',')
        end
        # Keep original spline coordinates for interpretation and later prediction.
        for (label,value) in ["knots"=>knots,"H"=>basis.H,"C"=>prepared.C,
                              "fixed_coefficients"=>fit.beta,"random_coefficients"=>fit.r,
                              "spline_coefficients"=>prepared.C*vcat(fit.r,fit.beta[2:end])]
            writedlm(joinpath(dst,label*".csv"),value,',')
        end
        push!(rows,(;case=name,lambda=fit.lambda,sigma2=fit.sigma2,tau2=fit.tau2,
          theta=fit.theta,random_directions=size(prepared.Z,2),
          variance_parameters=length(fit.prepared.model.θ),converged=true,method="REML",
          status=fit.status,objective=fit.objective,evaluations=fit.evaluations))
        println("Julia REML: ",name," | lambda = ",fit.lambda)
    end
    open(joinpath(out,"fit_summary.tsv"),"w") do io
        println(io,join(string.(keys(first(rows))),'\t'))
        for row in rows
            println(io,join(string.(values(row)),'\t'))
        end
    end
    open(joinpath(out,"versions.tsv"),"w") do io
        println(io,"key\tvalue")
        println(io,"Julia\t",VERSION)
        println(io,"MixedModels\t",pkgversion(MixedModels))
        println(io,"StatsModels\t",pkgversion(StatsModels))
        println(io,"method\tREML")
        println(io,"random_covariance\ttau_squared * identity")
        println(io,"residual_covariance\tsigma_squared * identity")
    end
end
main(ARGS)
