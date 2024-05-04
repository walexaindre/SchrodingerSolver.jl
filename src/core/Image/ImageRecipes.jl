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