const USE_GPU = true
using GLMakie
using ParallelStencil
using ParallelStencil.FiniteDifferences3D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end

@parallel function diffusion3D_step!(T2, T, Ci, lam, dt, dx, dy, dz)
    @inn(T2) = @inn(T) +
               dt * (lam * @inn(Ci) *
                     (@d2_xi(T) / dx^2 + @d2_yi(T) / dy^2 + @d2_zi(T) / dz^2))
    return
end

function diffusion3D()
    # Physics
    lam = 1.0                                        # Thermal conductivity
    cp_min = 1.0                                        # Minimal heat capacity
    lx, ly, lz = 10.0, 10.0, 10.0                           # Length of domain in dimensions x, y and z.

    # Numerics
    nx, ny, nz = 256, 256, 256                              # Number of gridpoints dimensions x, y and z.
    nt = 10000                                        # Number of time steps
    dx = lx / (nx - 1)                                  # Space step in x-dimension
    dy = ly / (ny - 1)                                  # Space step in y-dimension
    dz = lz / (nz - 1)                                  # Space step in z-dimension

    # Array initializations
    T = @zeros(nx, ny, nz)
    T2 = @zeros(nx, ny, nz)
    Ci = @zeros(nx, ny, nz)

    # Initial conditions (heat capacity and temperature with two Gaussian anomalies each)
    Ci .= 1.0 ./ (cp_min .+ Data.Array([5 * exp(-(((ix - 1) * dx - lx / 1.5))^2 -
                                                (((iy - 1) * dy - ly / 2))^2 -
                                                (((iz - 1) * dz - lz / 1.5))^2) +
                                        5 * exp(-(((ix - 1) * dx - lx / 3.0))^2 -
                                                (((iy - 1) * dy - ly / 2))^2 -
                                                (((iz - 1) * dz - lz / 1.5))^2)
                                        for ix in 1:size(T, 1), iy in 1:size(T, 2), iz in 1:size(T, 3)]))
    T .= Data.Array([100 * exp(-(((ix - 1) * dx - lx / 2) / 2)^2 -
                               (((iy - 1) * dy - ly / 2) / 2)^2 -
                               (((iz - 1) * dz - lz / 3.0) / 2)^2) +
                     50 * exp(-(((ix - 1) * dx - lx / 2) / 2)^2 -
                              (((iy - 1) * dy - ly / 2) / 2)^2 -
                              (((iz - 1) * dz - lz / 1.5) / 2)^2)
                     for ix in 1:size(T, 1), iy in 1:size(T, 2), iz in 1:size(T, 3)])
    T2 .= T                                                 # Assign also T2 to get correct boundary conditions.

    # Time loop
    dt = min(dx^2, dy^2, dz^2) * cp_min / lam / 8.1                 # Time step for the 3D Heat diffusion
    #for it = 1:nt
    #    
    #end
    println("Time step: ", dt)
    # Plotting
    fig = Figure(; size = (800, 800))
    ax = Axis3(fig[1, 1])
    Makie.record(fig, "file.mp4", 1:nt; framerate = 20) do _
        @parallel diffusion3D_step!(T2, T, Ci, lam, dt, dx, dy, dz)
        T, T2 = T2, T
        Tx = Array(T)
        empty!(ax)
        volume!(ax, Tx; colormap = :thermal)
    end
end

diffusion3D()

