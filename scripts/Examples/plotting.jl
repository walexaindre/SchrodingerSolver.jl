#Parameters to put figures in the right place and size

julia> save("a.png",figure1,px_per_unit=3)


function solvertime(Stats,ounit)
    p=_solvertime(Stats,ounit,fontsize=50)

    p.axis.xlabel = "τ"
    p.axis.ygridvisible = false
    p.axis.xgridvisible = false
    p.axis.yminorgridvisible = false
    p.axis.ylabel = "Time"
    p.axis.yticks = WilkinsonTicks(6,k_min=5)
    p.axis.xticks = WilkinsonTicks(6,k_min=5)
    p.axis.xlabelsize = p.axis.ylabelsize = 24
    p.axis.xticklabelsize=p.axis.yticklabelsize=24
    p.axis.ytickformat = xs -> [string(transform(v,ounit)) for v in xs]
    p    
end
