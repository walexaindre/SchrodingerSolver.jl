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
    rect = lift(p.I,p.J) do I,J
        ext = extrema(J)
        [Rect(i-1,ext[2]-j,1,1) for (i,j) in zip(I,J)]
    end

    poly!(p,rect,color=p.V,colormap=p.colormap,inspectable=p.inspectable,visible=p.visible)
end

export sparsitypattern