@recipe(SparsityPattern,I,J,V) do scene
    Attributes(
        colormap=theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function convert_arguments(::Type{<:SparsityPattern}, I, J, V)
    return I,J,V
end

function Makie.plot!(p::SparsityPattern)

    updated_values = p.V[]

    unique_values = unique(updated_values)
    subs = unique_values .=> 1:length(unique_values)
    replace!(updated_values,subs...)

    rect = lift(p.I,p.J) do I,J
        ext = extrema(J)
        [Rect(i-1,ext[2]-j,1,1) for (i,j) in zip(I,J)]
        
    end

    poly!(p,rect,color=updated_values,colormap=p.colormap,inspectable=p.inspectable,visible=p.visible)
end

export sparsitypattern


@recipe(Stencil,offsets) do scene
    Attributes(
        colormap=theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function convert_arguments(::Type{<:Stencil}, offsets)
    return offsets
end

function Makie.plot!(p::Stencil)
    color = zeros(Int,0)

    rect = lift(p.offsets) do offsets
        dims = ndims(offsets)
        stencil_center = ntuple(Returns(0.0),dims)
        offbydim = offsets.offsets
        rectangle_type = typeof(Rect(stencil_center...,ntuple(Returns(0.8),dims)...))
        out::Vector{rectangle_type} = []

        for (dim,off) in zip(1:dims,offbydim)

            for ele in off.offsets
                newidx = Base.setindex(stencil_center,ele,dim)
                push!(out,Rect(newidx...,ntuple(Returns(0.80),dims)...))
                push!(color,abs(ele))
            end
        end
        out
    end
    println(rect)
    poly!(p,rect,color=color,colormap=p.colormap,inspectable=p.inspectable,visible=p.visible)
end

export stencil



@recipe(ApproximateSparsityPattern,I,J,V) do scene
    Attributes(
        colorscale=identity,
        colormap=theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function convert_arguments(::Type{<:ApproximateSparsityPattern}, I, J, V)
    return I,J,V
end

function Makie.plot!(p::ApproximateSparsityPattern)

    updated_values = copy(p.V[])
    updated_values = map(abs,updated_values)
    updated_values = map(x->round(x,sigdigits=2),updated_values)


    unique_values = unique(updated_values)
    #subs = unique_values .=> 1:length(unique_values)
    #replace!(updated_values,subs...)

    rect = lift(p.I,p.J) do I,J
        ext = extrema(J)
        [Rect(i-1,ext[2]-j,1,1) for (i,j) in zip(I,J)]
        
    end
    poly!(p,rect, colorscale=p.colorscale , colorrange=(minimum(unique_values),maximum(unique_values)) ,color=updated_values,colormap=p.colormap,inspectable=p.inspectable,visible=p.visible)
end

export approximatesparsitypattern


@recipe(_SystemEnergy,Stats) do scene
    Attributes(
        comparewithstartup = true,
        color = theme(scene,:color),
        colormap = theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function toscientific(x)
    result = @sprintf("%.1E",x)
    splitres = split(result,"E")

    rich("$(splitres[1])×10",superscript("$(splitres[2])"))
end

function Makie.plot!(p::_SystemEnergy)
    Stats = p[1]
    τ = Stats[].τ


    comparewithzero = p.comparewithstartup
    points = Observable(Float64[])    
    timerange = Observable(range(0,step=τ,length=4))

    function update_plot(Stats)
        empty!(points[])
        maxindex = Stats[].store_index-1
        timerange[] = range(Stats[].τ,step=Stats[].τ,length=maxindex)

        if comparewithzero[]
            append!(points[], cbrt.(calculate_diff_system_energy(Stats[])))
        else
            append!(points,Stats[].system_energy[1:maxindex])
        end
        
    end

    #Makie.Observables.onany(update_plot,Stats)

    update_plot(Stats)
    
    Makie.scatter!(p,timerange,points,  color = :blue, colormap = p.colormap, inspectable = p.inspectable, visible = p.visible)
    p
end

function systemenergy(Stats)
    p = _systemenergy(Stats)
    p.axis.yscale = Makie.pseudolog10
    p.axis.xlabel = "τ"
    #p.axis.ylabel = 
    p.axis.ytickformat = xs -> [toscientific(v^3) for v in xs]
    ylims!(p.axis,(cbrt(1e-17),cbrt(1e-12)))
    p
end


export systemenergy


@recipe(SystemND,Memory,Grid) do scene
    Attributes(
        colorscale=identity,
        colormap=theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function Makie.plot!(Sys::SystemND)
    Grid = Sys[2][]
    Memory = Sys[1][]

    current_state = Memory.current_state

    points = abs2.(current_state) 

    points = sum(points,dims=2)
    
    points = points |> Array

    z = reshape(points,size(Grid))
    x = get_range(Grid,1)
    y = get_range(Grid,2)

    surface!(Sys,x,y, z;colormap=Sys.colormap,inspectable=Sys.inspectable,visible=Sys.visible)
end

export systemnd


@recipe(_ExecutionTime,Stats) do scene
    Attributes(
        color = theme(scene,:color),
        colormap = theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function Makie.plot!(Sys::_ExecutionTime)
    Stats = Sys[1]
    τ = Stats[].τ

    points = Observable(Float64[])    
    timerange = Observable(range(τ,step=τ,length=4))
    avg = Observable(0.0)

    function update_plot(Stats)
        empty!(points[])
        maxindex = Stats[].store_index-1
        timerange[] = range(Stats[].τ,step=Stats[].τ,length=maxindex)
        append!(points[],Stats[].step_time[1:maxindex])
        avg[] = sum(points[])/maxindex
    end

    update_plot(Stats)
    
    Makie.scatter!(Sys,timerange,points,  color = RGBf(0.0039,0.239216,0.5333) , colormap = Sys.colormap, inspectable = Sys.inspectable, visible = Sys.visible)
    Makie.ablines!(Sys,avg[],0,color= RGBf(0.698,0.168,0.0745),linestyle=:dash,linewidth=2,fontsize=Sys.fontsize)
    Sys

end

function executiontime(Stats,ounit)
    p=_executiontime(Stats,ounit,fontsize=50)

    p.axis.xlabel = rich("t",subscript("n"))
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


export executiontime


@recipe(_SolverTime,Stats,ounit) do scene
    Attributes(
        fontsize = theme(scene,:fontsize),
        color = theme(scene,:color),
        colormap = theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function transform(x,unit)
    r = uconvert(unit,x*u"s")
    r = round(typeof(r),r,digits=0)
    r = Quantity{Int64}[r][1]
end

function Makie.plot!(Sys::_SolverTime)
    Stats = Sys[1]
    τ = Stats[].τ
    outputunit = Sys[2][]

    points = Observable(Float64[])    
    timerange = Observable(range(τ,step=τ,length=4))
    avg = Observable(0.0)
    function update_plot(Stats)
        empty!(points[])
        maxindex = Stats[].store_index-1
        timerange[] = range(Stats[].τ,step=Stats[].τ,length=maxindex)
        append!(points[],Stats[].solver_time[1:maxindex])
        avg[] = sum(points[])/maxindex
    end

    update_plot(Stats)
    
    Makie.scatter!(Sys,timerange,points,  color = :blue, colormap = Sys.colormap, inspectable = Sys.inspectable, visible = Sys.visible)
    Makie.ablines!(Sys,avg[],0,color=:red,linestyle=:dash,linewidth=2,fontsize=Sys.fontsize)
    Sys

end


function solvertime(Stats,ounit)
    p=_solvertime(Stats,ounit,fontsize=50)

    p.axis.xlabel = rich("t",subscript("n"))
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


export solvertime

@recipe(SolverIterations,Stats) do scene
    Attributes(
        color = theme(scene,:color),
        colormap = theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function Makie.plot!(Sys::SolverIterations)
    Stats = Sys[1]
    τ = Stats[].τ

    points = Observable(Float64[])    
    timerange = Observable(range(τ,step=τ,length=4))

    function update_plot(Stats)
        empty!(points[])
        maxindex = Stats[].store_index-1
        timerange[] = range(Stats[].τ,step=Stats[].τ,length=maxindex)

        append!(points[],Stats[].solver_iterations[1:maxindex])
    end

    update_plot(Stats)
    
    Makie.scatter!(Sys,timerange,points,  color = :blue, colormap = Sys.colormap, inspectable = Sys.inspectable, visible = Sys.visible)
    Sys

end

export solveriterations