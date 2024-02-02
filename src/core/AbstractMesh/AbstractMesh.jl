include("AbstractMesh2.jl")
include("AbstractMeshTypes.jl")
include("AbstractMeshBase.jl")

export linear_indexing, step, step_x, step_y, step_z, linear_size

export MetaMesh2D, AbstractMesh







@inline function step(index::Int, lim::Int)
    r1 = (index - 1 < 1) ? lim : index - 1
    r2 = (index + 1 > lim) ? 1 : index + 1
    (r1, r2)
end

@inline function step(index::Int, lim::Int, step_size::Int)
    offset = index - 1
    r1 = mod(offset - step_size, lim) + 1
    r2 = mod(offset + step_size, lim) + 1
    (r1, r2)
end


@inline function step_x(index::Int, metadata::MetaMesh1D)
    step(index, metadata.M)
end

@inline function step_x(index::Int, metadata::MetaMesh2D)
    step(index, metadata.M)
end

@inline function step_x(index::Int, metadata::MetaMesh3D)
    step(index, metadata.M)
end

@inline function step_y(index::Int, metadata::MetaMesh2D)
    step(index, metadata.N)
end

@inline function step_y(index::Int, metadata::MetaMesh3D)
    step(index, metadata.N)
end

@inline function step_z(index::Int, metadata::MetaMesh3D)
    step(index, metadata.L)
end

@inline function step_x(index::Int, step_size::Int, metadata::MetaMesh)
    step(index, metadata.M, step_size)
end

@inline function step_y(index::Int, step_size::Int, metadata::MetaMesh2D)
    step(index, metadata.N, step_size)
end

@inline function step_y(index::Int, step_size::Int, metadata::MetaMesh3D)
    step(index, metadata.N, step_size)
end

@inline function step_z(index::Int, step_size::Int, metadata::MetaMesh3D)
    step(index, metadata.L, step_size)
end


include("AbstractMesh1D.jl")
include("AbstractMesh2D.jl")
include("AbstractMesh3D.jl")


@inline function get_stencil(index::Int,start_depth::Int,end_depth::Int,metadata::MetaMesh2D)
    depth =  end_depth-start_depth+1

    out = zeros(Int64,1+depth*4)

    out[1] = index

    for dist in start_depth:end_depth
        @show step_x(index,dist,metadata)...,step_y(index,dist,metadata)...
        out[(4*(dist-1)+2):(4*(dist)+1)] = [step_x(index,dist,metadata)...,step_y(index,dist,metadata)...]
    end
    out    
end

@inline function get_stencil(index::Int,start_depth::Int,end_depth::Int,metadata::MetaMesh3D)
    depth =  end_depth-start_depth+1

    out = zeros(Int64,1+depth*6)

    out[1] = index

    for dist in start_depth:end_depth
        out[(6*(dist-1)+2):(6*(dist)+1)] = [step_x(index,dist,metadata)...,step_y(index,dist,metadata)...,step_z(index,dist,metadata)...]
    end
    out    
end




@inline function get_linear_stencil(index::Int,start_depth::Int,end_depth::Int,metadata::MetaMesh2D)
    depth =  end_depth-start_depth+1

    out = zeros(Int64,1+depth*4)

    out[1] = index
    
    x,y = linear_indexing(index,metadata)

    for dist in 1:depth
        step_mov = start_depth+dist-1

        lx,rx =step_x(x,step_mov,metadata)
        olx = linear_indexing(lx,y,metadata)
        orx = linear_indexing(rx,y,metadata)

        ly,ry =step_y(y,step_mov,metadata)
        oly = linear_indexing(x,ly,metadata)
        ory = linear_indexing(x,ry,metadata)

        out[(4*dist-2):(4*dist+1)] = [olx,orx,oly,ory]
    end
    out  
end

@inline function get_linear_stencil(index::Int,start_depth::Int,end_depth::Int,metadata::MetaMesh3D)
    depth =  end_depth-start_depth+1
    outdim = 1+depth*6
    #@show outdim,depth,end_depth,start_depth
    out = zeros(Int64,outdim)

    out[1] = index

    x,y,z = linear_indexing(index,metadata)

    for dist in 1:depth
        step_mov = start_depth+dist-1

        lx,rx =step_x(x,step_mov,metadata)
        olx = linear_indexing(lx,y,z,metadata)
        orx = linear_indexing(rx,y,z,metadata)

        ly,ry =step_y(y,step_mov,metadata)
        oly = linear_indexing(x,ly,z,metadata)
        ory = linear_indexing(x,ry,z,metadata)

        lz,rz =step_z(z,step_mov,metadata)
        olz = linear_indexing(x,y,lz,metadata)
        orz = linear_indexing(x,y,rz,metadata)


        out[(6*dist-4):(6*dist+1)] = [olx,orx,oly,ory,olz,orz]
    end
    out 
end


#Base methods






#################################################################