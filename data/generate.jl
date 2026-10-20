using Random, DelimitedFiles
include(joinpath(@__DIR__,"..","spline_basis","src","BSpline.jl"))
using .BSplineBasis

function main(args)
    out=isempty(args) ? joinpath(@__DIR__,"inputs") : abspath(args[1])
    isdir(out)&&!isempty(readdir(out))&&error("Data already exist; reuse them or choose a new directory")
    rng=MersenneTwister(42)
    x=sort(vcat(-2.,4 .* rand(rng,198).-2,2.))
    truth=sin.(1.7 .* x)+0.15 .* x
    y=truth+0.35 .* randn(rng,length(x))
    repeated=sort(repeat(collect(range(-2,2;length=40));inner=5))
    repeated_truth=sin.(1.7 .* repeated)+0.15 .* repeated
    repeated_y=repeated_truth+0.35 .* randn(rng,length(repeated))
    cases=[
      (name="uniform_k8",x=x,y=y,truth=truth,t=default_knots(x,8),mode="explicit"),
      (name="default_k8",x=x,y=y,truth=truth,t=default_knots(x,8),mode="default"),
      (name="uneven_k10",x=x,y=y,truth=truth,t=extend_knots([-2.004,-1.7,-1.1,-0.55,-0.1,0.6,1.3,2.004]),mode="explicit"),
      (name="duplicates_k6",x=repeated,y=repeated_y,truth=repeated_truth,t=extend_knots([-2.004,-0.8,0.6,2.004]),mode="explicit"),
      (name="small_k5",x=x,y=y,truth=truth,t=default_knots(x,5),mode="explicit")]
    mkpath(out)
    open(joinpath(out,"cases.tsv"),"w") do io
        println(io,"case\tk\tknot_mode")
        for c in cases
            dst=joinpath(out,c.name);mkpath(dst)
            open(joinpath(dst,"data.csv"),"w") do f
                println(f,"x,y,truth")
                writedlm(f,hcat(c.x,c.y,c.truth),',')
            end
            writedlm(joinpath(dst,"knots.csv"),c.t,',')
            writedlm(joinpath(dst,"evaluation_x.csv"),collect(range(-2.7,2.7;length=61)),',')
            println(io,c.name,'\t',length(c.t)-4,'\t',c.mode)
        end
    end
    write(joinpath(out,"ORIGIN.txt"),"Generated once in Julia $(VERSION), seed 42.\nBoth checks use uniform_k8/data.csv and uniform_k8/knots.csv for the main example.\nFor cubic B-splines k is the basis size; each full knot vector has k+4 entries. The other knot cases reuse the same x/y values. duplicates_k6 is a separate stress-test dataset.\n")
    println("Shared Julia data saved: ",out)
end
main(ARGS)
