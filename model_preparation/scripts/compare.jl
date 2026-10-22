using DelimitedFiles, LinearAlgebra
include(joinpath(@__DIR__,"..","..","spline_basis","src","BSpline.jl"))
using .BSplineBasis
include(joinpath(@__DIR__,"..","src","ModelInputs.jl"))
using .ModelInputs

function main(args)
    root=dirname(@__DIR__)
    inputs=length(args)>=1 ? abspath(args[1]) : joinpath(root,"..","data","inputs")
    reference=length(args)>=2 ? abspath(args[2]) : joinpath(root,"outputs","shared","r_reference")
    out=length(args)>=3 ? abspath(args[3]) : joinpath(root,"outputs","shared","julia")
    isdir(out)&&!isempty(readdir(out))&&error("Julia output directory is not empty; choose a new run directory")
    cases=readdlm(joinpath(inputs,"cases.tsv"),'\t',String;skipstart=1)
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
        converted=prepare_inputs(basis.B,basis.S;expected_rank=k-2)
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
            println(io,"spline_type\tbs")
            println(io,"degree\t3")
            println(io,"penalty_derivative_order\t2")
            println(io,"observations\t",length(x))
            println(io,"penalty_rank\t",converted.penalty_rank)
            println(io,"random_covariance\tr_squared * identity")
            println(io,"residual_covariance\tsigma_squared * identity")
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
    println("$passed/$(length(records)) combined-workflow checks passed; largest scaled error = $worst")
    all(r->r.passed,records)||error("Combined workflow comparisons failed; see checks.tsv")
end
main(ARGS)
