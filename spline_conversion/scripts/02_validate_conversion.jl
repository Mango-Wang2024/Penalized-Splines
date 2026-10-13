using DelimitedFiles, LinearAlgebra, Random
include(joinpath(@__DIR__,"..","src","SplineConversion.jl"))
using .SplineConversion

function main(args)
    root = dirname(@__DIR__)
    ref = isempty(args) ? joinpath(root,"outputs","shared","r_reference") : abspath(args[1])
    isfile(joinpath(ref,"metadata.tsv")) || error("Missing R reference. Run scripts/01_inspect_gamm_inputs.R first.")
    out = length(args)>=2 ? abspath(args[2]) : joinpath(root,"outputs","shared","julia_conversion")
    isdir(out) && !isempty(readdir(out)) && error("Output directory is not empty. Supply a new directory.")
    metadata = readdlm(joinpath(ref,"metadata.tsv"),'\t',String;skipstart=1)
    md = Dict(metadata[i,1]=>metadata[i,2] for i in axes(metadata,1))
    expected_rank = parse(Int,md["penalty_rank"])
    load(name) = Matrix{Float64}(readdlm(joinpath(ref,name*".csv"),',',Float64))
    B,S = load("B"),load("S")
    X,X0r,Zr,Cr = load("X_ordinary"),load("X0_r"),load("Z_r"),load("C_r")
    Xactual,Zactual = load("X_fixed_lme"),load("Z_lme")
    c = convert_spline(B,S;expected_rank)
    p = size(B,2)
    target_penalty = Matrix(Diagonal(vcat(ones(expected_rank),zeros(p-expected_rank))))
    checks = Pair{String,Float64}[]
    push!(checks,"curve matrix identity"=>relative_error(B*c.C,hcat(c.Z,c.X0)))
    push!(checks,"penalty identity"=>relative_error(c.C'*S*c.C,target_penalty))
    push!(checks,"all basis directions retained"=>relative_error(space_projector(B),space_projector(hcat(c.Z,c.X0))))
    push!(checks,"fixed-effect space vs mgcv conversion"=>relative_error(space_projector(hcat(X,c.X0)),space_projector(hcat(X,X0r))))
    push!(checks,"random covariance vs mgcv conversion"=>relative_error(c.Z*c.Z',Zr*Zr'))
    push!(checks,"fixed-effect space vs actual lme inputs"=>relative_error(space_projector(hcat(X,c.X0)),space_projector(Xactual)))
    push!(checks,"random covariance vs actual lme inputs"=>relative_error(c.Z*c.Z',Zactual*Zactual'))
    # Map coefficients between representations rather than compare eigenvectors.
    rng = MersenneTwister(42)
    coeff = randn(rng,p)
    z = c.C*coeff
    r_coeff = Cr\z
    push!(checks,"curve after mapping coefficients to R coordinates"=>relative_error(hcat(c.Z,c.X0)*coeff,hcat(Zr,X0r)*r_coeff))
    push!(checks,"penalty after mapping coefficients to R coordinates"=>relative_error([sum(abs2,coeff[1:expected_rank])],[sum(abs2,r_coeff[1:expected_rank])]))
    rank(c.C)==p || error("Coefficient transformation lost directions")
    threshold = 1e-8
    mkpath(out)
    open(joinpath(out,"validation.tsv"),"w") do io
        println(io,"check\trelative_error\tpassed")
        for (name,err) in checks
            println(io,name,'\t',err,'\t',err<threshold)
        end
    end
    for (name,M) in ["X0_julia"=>c.X0,"Z_julia"=>c.Z,"C_julia"=>c.C]
        writedlm(joinpath(out,name*".csv"),M,',')
    end
    lines = ["STAGE 1: CONVERSION ONLY", "",
      "R reference: "*ref,"mgcv version: "*md["mgcv_version"],
      "Julia version: "*string(VERSION),
      "Basis dimensions: "*string(size(B)),"Penalty rank: "*string(c.penalty_rank),
      "Unpenalized spline columns: "*string(c.null_rank),"",
      "Raw eigenvector signs are not compared. Fixed spaces and Z*Z' are compared.",
      "Relative errors use max(1, norm(left), norm(right)) as the denominator.",""]
    append!(lines,[string(err<threshold ? "PASS" : "FAIL"," | ",name," | error=",err) for (name,err) in checks])
    append!(lines,["","R and Julia convert the exact same R-generated basis and penalty.",
      "This check validates conversion only. Independent construction and the combined workflow have separate checks.",
      "Julia fitting, lambda estimation and final fitted predictions remain later stages."])
    write(joinpath(out,"RESULTS_EXPLAINED.txt"),join(lines,"\n")*"\n")
    println(join(lines,"\n"))
    all(last(c)<threshold for c in checks) || error("Conversion checks failed; see validation.tsv")
end
main(ARGS)
