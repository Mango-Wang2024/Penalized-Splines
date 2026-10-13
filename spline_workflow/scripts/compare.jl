using DelimitedFiles, LinearAlgebra
include(joinpath(@__DIR__,"..","..","spline_basis","src","NaturalCubic.jl"))
using .NaturalCubic
include(joinpath(@__DIR__,"..","..","spline_conversion","src","SplineConversion.jl"))
using .SplineConversion

function main(args)
    root=dirname(@__DIR__)
    inputs=length(args)>=1 ? abspath(args[1]) : joinpath(root,"..","data","inputs")
    reference=length(args)>=2 ? abspath(args[2]) : joinpath(root,"outputs","shared","r_reference")
    out=length(args)>=3 ? abspath(args[3]) : joinpath(root,"outputs","shared","julia")
    isdir(out)&&!isempty(readdir(out))&&error("Julia output directory is not empty; choose a new run directory")
    cases=readdlm(joinpath(inputs,"cases.tsv"),'\t',String;skipstart=1)
    meta=readdlm(joinpath(reference,"metadata.tsv"),'\t',String;skipstart=1)
    md=Dict(meta[i,1]=>meta[i,2] for i in axes(meta,1))
    threshold=1e-8
    records=NamedTuple[]
    function check(case,label,a,b)
        size(a)==size(b)||error("$case: $label dimensions disagree")
        err=relative_error(a,b)
        push!(records,(;case,check=label,scaled_error=err,max_abs_error=maximum(abs,a-b;init=0.0),passed=err<threshold))
    end
    mkpath(out)
    for i in axes(cases,1)
        name=cases[i,1];k=parse(Int,cases[i,2]);mode=cases[i,3]
        p=joinpath(inputs,name);rpath=joinpath(reference,name)
        data=readdlm(joinpath(p,"data.csv"),',',Float64;skipstart=1)
        x,y=data[:,1],data[:,2]
        # Construction reads the shared raw inputs only, never R basis outputs.
        knots=mode=="default" ? default_knots(x,k) : vec(readdlm(joinpath(p,"knots.csv"),',',Float64))
        basis=build_spline(x,knots)
        converted=convert_spline(basis.B,basis.S;expected_rank=k-2)
        X_fixed=hcat(ones(length(x)),converted.X0)
        load(name)=Matrix{Float64}(readdlm(joinpath(rpath,name*".csv"),',',Float64))
        Xr,X0r,Zr,Cr=load("X_fixed"),load("X0"),load("Z"),load("C")
        target=Matrix(Diagonal(vcat(ones(k-2),zeros(1))))
        check(name,"fixed space: complete Julia vs complete R",space_projector(X_fixed),space_projector(Xr))
        check(name,"random covariance: complete Julia vs complete R",converted.Z*converted.Z',Zr*Zr')
        check(name,"curve preserved by Julia conversion",basis.B*converted.C,hcat(converted.Z,converted.X0))
        check(name,"penalty preserved by Julia conversion",converted.C'*basis.S*converted.C,target)
        check(name,"all basis directions retained",space_projector(basis.B),space_projector(hcat(converted.Z,converted.X0)))
        # Compare all coefficient directions after mapping into R coordinates.
        mapping=Cr\converted.C
        check(name,"curve in common coefficient coordinates",hcat(converted.Z,converted.X0),hcat(Zr,X0r)*mapping)
        check(name,"penalty in common coefficient coordinates",mapping'*target*mapping,target)
        rank(converted.C)==k-1||error("$name: conversion lost coefficient directions")
        dst=joinpath(out,name);mkpath(dst)
        for (label,M) in ["y"=>reshape(y,:,1),"B"=>basis.B,"S"=>basis.S,
                         "X_fixed"=>X_fixed,"X0"=>converted.X0,"Z"=>converted.Z,
                         "C"=>converted.C,"knots"=>reshape(knots,:,1),"H"=>basis.H,
                         "penalty_scale"=>reshape([basis.penalty_scale],1,1)]
            writedlm(joinpath(dst,label*".csv"),M,',')
        end
        open(joinpath(dst,"model_settings.tsv"),"w") do io
            println(io,"key\tvalue")
            println(io,"spline_type\tcr")
            println(io,"observations\t",length(x))
            println(io,"penalty_rank\t",converted.penalty_rank)
            println(io,"random_covariance\tsigma_u_squared * identity")
            println(io,"residual_covariance\tsigma_error_squared * identity")
            println(io,"fitting_status\tnot fitted; lambda not estimated")
        end
        println("Combined workflow: ",name)
    end
    open(joinpath(out,"checks.tsv"),"w") do io
        println(io,"case\tcheck\tscaled_error\tmax_abs_error\tpassed")
        for r in records
            println(io,r.case,'\t',r.check,'\t',r.scaled_error,'\t',r.max_abs_error,'\t',r.passed)
        end
    end
    passed=count(r->r.passed,records)
    worst=maximum(r.scaled_error for r in records)
    result=["# Combined construction and conversion","",
      "R constructs and converts the spline. Julia independently constructs and converts it from the same raw inputs.","",
      "**Result: $passed/$(length(records)) checks passed.** Largest scaled error: $worst. Threshold: $threshold.","",
      "Versions: $(md["R_version"]), mgcv $(md["mgcv_version"]), Julia $(VERSION).","",
      "| Case | Checks passed | Largest scaled error |","| --- | --- | --- |"]
    for name in cases[:,1]
        rs=filter(r->r.case==name,records)
        push!(result,"| $name | $(count(r->r.passed,rs))/$(length(rs)) | $(maximum(r.scaled_error for r in rs)) |")
    end
    append!(result,["","Checks compare the fixed space and random covariance, preserve curve and penalty, and map all coefficient directions into common coordinates. Raw eigenvectors are not compared.","",
      "Scaled errors use norm(a-b)/max(1, norm(a), norm(b)).","",
      "Each Julia case exports y, X_fixed and Z for Stage 2, with the matching B, S, coefficient map C and training settings. No coefficients or smoothing strength have been estimated.","",
      "[Detailed checks](julia/checks.tsv)"])
    write(joinpath(dirname(out),"RESULTS.md"),join(result,"\n")*"\n")
    println("$passed/$(length(records)) combined-workflow checks passed; largest scaled error = $worst")
    all(r->r.passed,records)||error("Combined workflow comparisons failed; see checks.tsv")
end
main(ARGS)
