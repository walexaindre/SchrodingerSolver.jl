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

    z = reshape(points,size(Grid))
    x = get_range(Grid,1)
    y = get_range(Grid,2)

    surface!(Sys,x,y, z;colormap=Sys.colormap,inspectable=Sys.inspectable,visible=Sys.visible)
end

export systemnd