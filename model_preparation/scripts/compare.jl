using DelimitedFiles, LinearAlgebra
include(joinpath(@__DIR__,"../../spline_basis/src/BSpline.jl"));using .BSplineBasis
include(joinpath(@__DIR__,"../src/ModelInputs.jl"));using .ModelInputs
include(joinpath(@__DIR__,"../../data/SharedData.jl"));using .SharedData
function main(args)
    length(args)==3 || error("Expected shared inputs, R reference and a new output directory")
    inputs,reference,out=abspath.(args)
    isdir(out)&&!isempty(readdir(out))&&error("Output directory is not empty")
    mkpath(out);records=NamedTuple[]
    function check(case,component,label,a,b)
        size(a)==size(b)||error("$case/$component: dimensions disagree")
        err=relative_error(a,b)
        push!(records,(;case,component,check=label,scaled_error=err,max_abs_error=maximum(abs,a-b;init=0.0),passed=err<1e-8))
    end
    for name in case_names(inputs)
        case=load_case(inputs,name)
        bases=[build_spline(t.x,t.mode=="default" ? default_knots(t.x,t.k) : t.knots) for t in case.terms]
        X=ones(length(case.y),1)
        p=prepare_additive(bases,X)
        Xr=readdlm(joinpath(reference,name,"X_fixed.csv"),',',Float64)
        check(name,"all","combined fixed space",space_projector(p.X_fixed),space_projector(Xr))
        dst=joinpath(out,name);mkpath(dst)
        writedlm(joinpath(dst,"X_fixed.csv"),p.X_fixed,',');writedlm(joinpath(dst,"y.csv"),case.y,',')
        cp(joinpath(case.path,"terms.tsv"),joinpath(dst,"terms.tsv"))
        for j in eachindex(bases)
            term=case.terms[j];basis=bases[j];c=p.components[j]
            load(label)=readdlm(joinpath(reference,name,term.name,label*".csv"),',',Float64)
            Zr,X0r,Cr=load("Z"),load("X0"),load("C")
            target=Matrix(Diagonal(vcat(ones(c.penalty_rank),zeros(c.null_rank))))
            check(name,term.name,"random covariance",c.Z*c.Z',Zr*Zr')
            check(name,term.name,"curve preserved",basis.B*c.C,hcat(c.Z,c.X0))
            check(name,term.name,"normalized penalty",c.C'*basis.S*c.C,target)
            check(name,term.name,"all spline directions retained",space_projector(basis.B),space_projector(hcat(c.Z,c.X0)))
            mapping=Cr\c.C
            check(name,term.name,"curve in common coordinates",hcat(c.Z,c.X0),hcat(Zr,X0r)*mapping)
            check(name,term.name,"penalty in common coordinates",mapping'*target*mapping,target)
            path=joinpath(dst,term.name);mkpath(path)
            for (label,M) in ["B_raw"=>basis.B_raw,"B"=>basis.B,"S"=>basis.S,"X0"=>c.X0,"Z"=>c.Z,
                              "C"=>c.C,"knots"=>reshape(basis.knots,:,1),"H"=>basis.H,"penalty_scale"=>reshape([basis.penalty_scale],1,1)]
                writedlm(joinpath(path,label*".csv"),M,',')
            end
        end
        open(joinpath(dst,"model_settings.tsv"),"w") do io
            println(io,"key\tvalue\ncomponents\t",length(bases),"\nobservations\t",length(case.y),
                    "\nrandom_covariance\tone r_j_squared * identity per component\nresidual_covariance\tsigma_squared * identity")
        end
        println("Combined preparation: ",name)
    end
    open(joinpath(out,"checks.tsv"),"w") do io
        println(io,"case\tcomponent\tcheck\tscaled_error\tmax_abs_error\tpassed")
        for r in records;println(io,join(string.(values(r)),'\t'));end
    end
    all(r.passed for r in records)||error("Combined preparation comparisons failed")
    println(length(records),"/",length(records)," combined preparation checks passed")
end
main(ARGS)
