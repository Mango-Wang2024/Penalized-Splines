using DelimitedFiles, LinearAlgebra, MixedModels, StatsModels
import BSplineKit
include(joinpath(@__DIR__,"../../spline_basis/src/BSpline.jl"));using .BSplineBasis
include(joinpath(@__DIR__,"../../model_preparation/src/ModelInputs.jl"));using .ModelInputs
include(joinpath(@__DIR__,"../src/SplineFit.jl"));using .SplineFit
include(joinpath(@__DIR__,"../../data/SharedData.jl"));using .SharedData
function write_columns(path,labels,values)
    open(path,"w") do io;println(io,join(labels,','));writedlm(io,values,',');end
end
function main(args)
    length(args)==2||error("Expected shared inputs and a new output directory")
    inputs,out=abspath.(args)
    isdir(out)&&!isempty(readdir(out))&&error("Output directory is not empty")
    mkpath(out);rows=NamedTuple[];estimates=NamedTuple[]
    for name in case_names(inputs)
        case=load_case(inputs,name)
        bases=[build_spline(t.x,t.mode=="default" ? default_knots(t.x,t.k) : t.knots) for t in case.terms]
        p=prepare_additive(bases;names=case.names)
        fit=fit_reml(case.y,p.X_fixed,p.Z_blocks;names=case.names)
        dst=joinpath(out,name);mkpath(dst)
        x=hcat((t.x for t in case.terms)...)
        write_columns(joinpath(dst,"fitted.csv"),vcat(case.names,["y","truth","fitted"]),hcat(x,case.y,case.truth,fit.fitted))
        curves=hcat((p.X_fixed[:,p.fixed_ranges[j]]*fit.b[p.fixed_ranges[j]]+p.Z_blocks[j]*fit.tau[j] for j in eachindex(bases))...)
        write_columns(joinpath(dst,"components.csv"),case.names,curves)
        grid=evaluation_grid(case)
        new=prediction_inputs(bases,p,[evaluate_basis(b.knots,g) for (b,g) in zip(bases,grid)])
        predicted=predict_response(fit,new.X_fixed,new.Z_blocks)
        write_columns(joinpath(dst,"grid.csv"),vcat(case.names,["fitted"]),hcat(grid...,predicted))
        for (label,M) in ["fixed_coefficients"=>fit.b,"Xf"=>p.X_fixed,"Xr"=>fit.prepared.Xr,
                          "Zf"=>hcat(p.Z_blocks...),"Zr"=>hcat(fit.prepared.Zr_blocks...),
                          "relative_loading"=>Matrix(Diagonal(vcat((fill(fit.theta[j],size(p.Z_blocks[j],2)) for j in eachindex(bases))...)))]
            writedlm(joinpath(dst,label*".csv"),M,',')
        end
        for j in eachindex(bases)
            path=joinpath(dst,case.names[j]);mkpath(path)
            for (label,M) in ["knots"=>bases[j].knots,"H"=>bases[j].H,"C"=>p.components[j].C,
                              "random_coefficients"=>fit.tau[j],
                              "spline_coefficients"=>p.components[j].C*vcat(fit.tau[j],fit.b[p.fixed_ranges[j]]),
                              "Zf"=>p.Z_blocks[j],"Zr"=>fit.prepared.Zr_blocks[j]]
                writedlm(joinpath(path,label*".csv"),M,',')
            end
            push!(estimates,(;case=name,component=case.names[j],lambda=fit.lambda[j],r2=fit.r2[j],r=fit.r[j],theta=fit.theta[j],random_directions=size(p.Z_blocks[j],2)))
        end
        push!(rows,(;case=name,components=length(bases),sigma2=fit.sigma2,fixed_columns=size(p.X_fixed,2),
          random_directions=sum(size(Z,2) for Z in p.Z_blocks),variance_parameters=length(fit.prepared.model.θ),
          converged=true,method="REML",status=fit.status,objective=fit.objective,evaluations=fit.evaluations))
        println("Julia joint REML: ",name," | lambda = ",fit.lambda)
    end
    for (file,records) in (("fit_summary.tsv",rows),("component_estimates.tsv",estimates))
        open(joinpath(out,file),"w") do io
            println(io,join(string.(keys(first(records))),'\t'))
            for r in records;println(io,join(string.(values(r)),'\t'));end
        end
    end
    open(joinpath(out,"versions.tsv"),"w") do io
        println(io,"key\tvalue\nJulia\t",VERSION,"\nBSplineKit\t",pkgversion(BSplineKit),"\nMixedModels\t",pkgversion(MixedModels),"\nStatsModels\t",pkgversion(StatsModels),"\nmethod\tREML")
    end
end
main(ARGS)
