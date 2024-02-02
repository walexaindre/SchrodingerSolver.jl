#1D

@inline function linear_indexing(index::Int, metadata::MetaMesh1D)
    index
end

@inline function linear_size(metadata::MetaMesh1D)
    metadata.M
end

@inline function get_stencil(index::Int,start_depth::Int,end_depth::Int,metadata::MetaMesh1D)
    depth =  end_depth-start_depth+1

    out = zeros(Int64,1+depth*2)

    out[1] = index

    for dist in start_depth:end_depth
        out[(2*(dist-1)+2):(2*(dist)+1)] = [step_x(index,dist,metadata)...]
    end
    out    
end

@inline function get_linear_stencil(index::Int,start_depth::Int,end_depth::Int,metadata::MetaMesh1D)
    depth =  end_depth-start_depth+1

    out = zeros(Int64,1+depth*2)

    out[1] = index

    x = linear_indexing(index,metadata)
    for dist in 1:depth

        lx,rx =step_x(x,start_depth+dist-1,metadata)
        olx = linear_indexing(lx,metadata)
        orx = linear_indexing(rx,metadata)

        out[(2*dist):(2*(dist)+1)] = [olx,orx]
    end
    out    
end


#Base methods




