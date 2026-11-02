using Random, DelimitedFiles
include(joinpath(@__DIR__,"..","spline_basis","src","BSpline.jl"))
using .BSplineBasis

# One entry per predictor; loaders and fitting functions do not assume a fixed count.
const COMPONENTS=[
    (name="x1",truth=x->sin.(1.7 .* x)+0.15 .* x),
    (name="x2",truth=x->0.5 .* cos.(1.2 .* x))]
function main(args)
    out=isempty(args) ? joinpath(@__DIR__,"inputs") : abspath(args[1])
    isdir(out)&&!isempty(readdir(out))&&error("Data already exist; supply a new directory")
    n=200
    x1=sort(vcat(-2.,4 .* rand(MersenneTwister(42),n-2).-2,2.))
    # Sorting only whole observation rows preserves the pairing of independent inputs.
    x2=4 .* rand(MersenneTwister(43),n).-2
    noise=0.35 .* randn(MersenneTwister(44),n)
    repeated=sort(repeat(collect(range(-2,2;length=40));inner=5))
    uneven=extend_knots([-2.004,-1.7,-1.1,-0.55,-0.1,0.6,1.3,2.004])
    cases=[
      (name="uniform_k8",x=[x1,x2],k=[8,8],mode=["explicit","explicit"]),
      (name="default_k8",x=[x1,x2],k=[8,8],mode=["default","default"]),
      (name="uneven_k10",x=[x1,x2],k=[10,6],mode=["uneven","explicit"]),
      (name="duplicates_k6",x=[repeated,x2],k=[6,8],mode=["explicit","explicit"]),
      (name="small_k5",x=[x1,x2],k=[5,5],mode=["explicit","explicit"])]
    mkpath(out)
    open(joinpath(out,"cases.tsv"),"w") do io
        println(io,"case")
        for case in cases
            println(io,case.name)
            dst=joinpath(out,case.name);mkpath(dst)
            truths=[COMPONENTS[j].truth(case.x[j]) for j in eachindex(COMPONENTS)]
            truth=sum(truths);y=truth+noise
            labels=vcat(getproperty.(COMPONENTS,:name),["y"],"truth_".*getproperty.(COMPONENTS,:name),["truth"])
            open(joinpath(dst,"data.csv"),"w") do f
                println(f,join(labels,','));writedlm(f,hcat(case.x...,y,truths...,truth),',')
            end
            axes=[collect(range(minimum(v),maximum(v);length=31)) for v in case.x]
            grid=[repeat(axes[1];outer=31),repeat(axes[2];inner=31)]
            grid_truth=sum(COMPONENTS[j].truth(grid[j]) for j in eachindex(COMPONENTS))
            open(joinpath(dst,"grid.csv"),"w") do f
                println(f,join(vcat(getproperty.(COMPONENTS,:name),["truth"]),','))
                writedlm(f,hcat(grid...,grid_truth),',')
            end
            open(joinpath(dst,"terms.tsv"),"w") do f
                println(f,"predictor\tk\tknot_mode")
                for j in eachindex(COMPONENTS)
                    key=COMPONENTS[j].name
                    mode=case.mode[j]=="default" ? "default" : "explicit"
                    println(f,key,'\t',case.k[j],'\t',mode)
                    knots=case.mode[j]=="uneven" ? uneven : default_knots(case.x[j],case.k[j])
                    writedlm(joinpath(dst,"knots_"*key*".csv"),knots,',')
                    writedlm(joinpath(dst,"evaluation_"*key*".csv"),collect(range(-2.7,2.7;length=61)),',')
                end
            end
        end
    end
    write(joinpath(out,"ORIGIN.txt"),"Generated once in Julia $(VERSION). Predictor seeds: 42 and 43; noise seed: 44.\nEach row contains x1,x2 and one response y.\ntruth=sin(1.7*x1)+0.15*x1+0.5*cos(1.2*x2); noise SD=0.35.\nBoth languages read the same saved values. Components are additive, without interactions.\n")
    println("Shared additive data saved: ",out)
end
main(ARGS)
