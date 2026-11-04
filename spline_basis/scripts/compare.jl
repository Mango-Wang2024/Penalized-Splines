using DelimitedFiles, LinearAlgebra, Statistics
import BSplineKit
include(joinpath(@__DIR__,"../../spline_basis/src/BSpline.jl"));using .BSplineBasis
include(joinpath(@__DIR__,"../../data/SharedData.jl"));using .SharedData
function main(args)
    length(args)==3 || error("Expected shared inputs, R reference and a new output directory")
    inputs,reference,out=abspath.(args)
    isdir(out)&&!isempty(readdir(out))&&error("Output directory is not empty")
    mkpath(out);records=NamedTuple[]
    function check(case,component,label,a,b)
        size(a)==size(b)||error("$case/$component: dimensions disagree for $label")
        err=norm(a-b)/max(1.0,norm(a),norm(b))
        push!(records,(;case,component,check=label,scaled_error=err,max_abs_error=maximum(abs,a-b;init=0.0),passed=err<1e-8))
    end
    for name in case_names(inputs)
        case=load_case(inputs,name)
        for term in case.terms
            knots=term.mode=="default" ? default_knots(term.x,term.k) : term.knots
            c=build_spline(term.x,knots);Be=evaluate_basis(c.knots,term.evaluation_x)
            values=["knots"=>reshape(c.knots,:,1),"D2"=>c.D2,"B_raw"=>c.B_raw,"S_raw"=>c.S_raw,
              "constraint"=>reshape(c.constraint,:,1),"penalty_scale"=>reshape([c.penalty_scale],1,1),
              "S_scaled"=>c.S_scaled,"H"=>c.H,"B"=>c.B,"S"=>c.S,"B_eval_raw"=>Be,"B_eval"=>Be*c.H]
            dst=joinpath(out,name,term.name);mkpath(dst)
            for (label,M) in values
                R=readdlm(joinpath(reference,name,term.name,label*".csv"),',',Float64)
                check(name,term.name,label,M,R);writedlm(joinpath(dst,label*".csv"),M,',')
            end
            check(name,term.name,"mean-zero prepared basis",vec(mean(c.B;dims=1)),zeros(term.k-1))
            println("Julia basis: ",name," ",term.name)
        end
    end
    open(joinpath(out,"checks.tsv"),"w") do io
        println(io,"case\tcomponent\tcheck\tscaled_error\tmax_abs_error\tpassed")
        for r in records;println(io,join(string.(values(r)),'\t'));end
    end
    open(joinpath(out,"versions.tsv"),"w") do io
        println(io,"key\tvalue\nJulia\t",VERSION,"\nBSplineKit\t",pkgversion(BSplineKit))
    end
    all(r.passed for r in records)||error("Basis comparisons failed")
    println(length(records),"/",length(records)," basis checks passed")
end
main(ARGS)
