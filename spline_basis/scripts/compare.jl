using DelimitedFiles, LinearAlgebra, Statistics
import BSplineKit
include(joinpath(@__DIR__,"..","src","BSpline.jl"))
using .BSplineBasis

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
        size(a)==size(b)||error("$case: $label dimensions disagree: $(size(a)) vs $(size(b))")
        err=norm(a-b)/max(1.0,norm(a),norm(b))
        push!(records,(;case,check=label,scaled_error=err,max_abs_error=maximum(abs,a-b;init=0.0),passed=err<threshold))
    end
    mkpath(out)
    for i in axes(cases,1)
        name=cases[i,1];k=parse(Int,cases[i,2]);mode=cases[i,3]
        p=joinpath(inputs,name);rpath=joinpath(reference,name)
        data=readdlm(joinpath(p,"data.csv"),',',Float64;skipstart=1)
        x=data[:,1]
        # These knots come from Julia's settings/default rule, never R's output.
        knots=mode=="default" ? default_knots(x,k) : vec(readdlm(joinpath(p,"knots.csv"),',',Float64))
        c=build_spline(x,knots)
        xe=vec(readdlm(joinpath(p,"evaluation_x.csv"),',',Float64))
        Be=evaluate_basis(c.knots,xe)
        values=["knots"=>reshape(c.knots,:,1),"D2"=>c.D2,"B_raw"=>c.B_raw,"S_raw"=>c.S_raw,
          "constraint"=>reshape(c.constraint,:,1),"penalty_scale"=>reshape([c.penalty_scale],1,1),
          "S_scaled"=>c.S_scaled,"H"=>c.H,"B"=>c.B,"S"=>c.S,
          "B_eval_raw"=>Be,"B_eval"=>Be*c.H]
        dst=joinpath(out,name);mkpath(dst)
        for (label,M) in values
            R=Matrix{Float64}(readdlm(joinpath(rpath,label*".csv"),',',Float64))
            check(name,label,M,R)
            writedlm(joinpath(dst,label*".csv"),M,',')
        end
        check(name,"mean-zero prepared basis",vec(mean(c.B;dims=1)),zeros(k-1))
        println("Julia comparison: ",name)
    end
    open(joinpath(out,"checks.tsv"),"w") do io
        println(io,"case\tcheck\tscaled_error\tmax_abs_error\tpassed")
        for r in records
            println(io,r.case,'\t',r.check,'\t',r.scaled_error,'\t',r.max_abs_error,'\t',r.passed)
        end
    end
    passed=count(r->r.passed,records)
    worst=maximum(r.scaled_error for r in records)
    println("$passed/$(length(records)) checks passed; largest scaled error = $worst")
    all(r->r.passed,records)||error("Matrix comparisons failed; see checks.tsv")
    open(joinpath(out,"versions.tsv"),"w") do io
        println(io,"key\tvalue")
        println(io,"Julia\t",VERSION)
        println(io,"BSplineKit\t",pkgversion(BSplineKit))
    end
end
main(ARGS)
