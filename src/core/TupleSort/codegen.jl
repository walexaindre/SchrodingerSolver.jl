using JSON
using Base.Iterators

const network_dir = pwd()*"/src/core/TupleSort/Network"
const output_dir = pwd()*"/src/core/TupleSort/"
#######################
#
#  This code generator is based on the work at bertdobbelaere.github.io
#
#####################

#######
#   normalize(x::T) where {T} = x<=zero(T) ? one(T)-T(2)*x : T(2)*x
#   min_max(x, y, lt = isless, by = identity) = lt(by(x), by(y)) ? (x, y) : (y, x)
#######



function generator!(io,ordercollection)
    cnum = sum(length, ordercollection)
    dnum = length(ordercollection)
    n = maximum(map(maximum,flatten(ordercollection))) + 1

    ##Info:
    info ="""
    number of elements: $n
                 depth: $dnum
                 steps: $cnum
    """
    println(info)

    varlist = ["i$idx" for idx in 1:n]
    indexvarx = ["x[$idx]" for idx in 1:n]
    indexvarxlist = join(indexvarx, ", ")
    function_name = " = min_max(:x,:y,lt,by)"
    output_params = ":o1,:o2"
    spacing = "    "
    vspacing = "\n\n"

    commentaryheader = "#Compares: $cnum Depth $dnum \n"
    finalfunctionheader = "@inline function swapsort(:args;lt=isless,by=identity)"
    auxiliaryfunction1 = "@inline swapsort(x::NTuple{$n,T};lt=isless,by=identity) where {T} = swapsort($indexvarxlist,lt=lt,by=by)\n"
    auxiliaryfunction2 = "@inline swapsort(::Type{Val{$n}},x::Vec; lt=isless, by=identity) where {Vec<:AbstractVector} = swapsort($indexvarxlist,lt=lt,by=by)\n"
    arglist = join(varlist, ", ")

    write(io, commentaryheader)
    write(io, replace(finalfunctionheader, ":args" => arglist))
    for seq in ordercollection
        write(io, "\n")
        for (var1, var2) in seq
            v1 = var1 + 1
            v2 = var2 + 1
            strv1 = varlist[v1]
            strv2 = varlist[v2]

            outputv = replace(output_params, ":o1" => strv1, ":o2" => strv2)
            funcall = replace(function_name, ":x" => strv1, ":y" => strv2)
            write(io, spacing)
            write(io, outputv)
            write(io, funcall)
            write(io, "\n")
        end
    end
    write(io, spacing)
    write(io, "return $arglist\n")
    write(io, "end")
    write(io, vspacing)
    write(io, auxiliaryfunction1)
    write(io, vspacing)
    write(io, auxiliaryfunction2)
    write(io, vspacing)
    io
end

function assemblycode()
    io=IOBuffer()

    for dir in readdir(network_dir)
        jsonfile = open(io->read(io,String),network_dir*"/"*dir)
        re = JSON.parse(jsonfile)
        #return re
        generator!(io,re["nw"])
    end
    io
end

re = assemblycode()
open(output_dir*"swapsort.jl","w+") do io
write(io,take!(re))
end