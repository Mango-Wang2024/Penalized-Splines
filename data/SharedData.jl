module SharedData
using DelimitedFiles
export case_names, load_case, evaluation_grid
case_names(root) = vec(readdlm(joinpath(root,"cases.tsv"),'\t',String;skipstart=1)[:,1])
function load_case(root,name)
    path=joinpath(root,name)
    values,header=readdlm(joinpath(path,"data.csv"),',',Float64;header=true)
    columns=Dict(string(key)=>values[:,i] for (i,key) in enumerate(vec(header)))
    spec=readdlm(joinpath(path,"terms.tsv"),'\t',String;skipstart=1)
    names=vec(spec[:,1])
    !isempty(names) && length(unique(names))==length(names) || error("Use distinct predictor names")
    all(occursin(r"^[A-Za-z][A-Za-z0-9_]*$",key) for key in names) || error("Invalid predictor name")
    terms=[(;name=names[i],k=parse(Int,spec[i,2]),mode=spec[i,3],x=columns[names[i]],
             knots=vec(readdlm(joinpath(path,"knots_"*names[i]*".csv"),',',Float64)),
             evaluation_x=vec(readdlm(joinpath(path,"evaluation_"*names[i]*".csv"),',',Float64)),
             truth=columns["truth_"*names[i]]) for i in axes(spec,1)]
    return (;name,path,terms,names,y=columns["y"],truth=columns["truth"])
end
"""Read the supplied prediction grid in predictor-list order."""
function evaluation_grid(case)
    values,header=readdlm(joinpath(case.path,"grid.csv"),',',Float64;header=true)
    keys=vec(string.(header))
    return [values[:,only(findall(==(name),keys))] for name in case.names]
end
end
