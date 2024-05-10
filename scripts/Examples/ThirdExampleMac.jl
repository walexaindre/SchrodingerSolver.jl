include("PDE_ex2D_3.jl")
using Dates
using Unitful
using JSON
using GLMakie
using Base.Iterators

stateidx = 1

function pde_parameters(ncomp, τv, hv, sord, tord)
    params = Dict()
    params["time"] = "$(Dates.now(UTC))"
    params["components"] = ncomp
    params["backend"] = string(backend)
    params["τ"] = τv
    params["h"] = hv
    params["TF"] = TF
    params["space_ord"] = sord
    params["time_ord"] = tord
    return params
end

const log_freq = 10
const backend = CUDABackendF64
const τ = 0.05
const h = 0.16
const TF = 25.0

const hvalues = [h, h / 2, h / 4, h / 8]
const τvalues = [τ, τ / 2, τ / 4, τ / 8]

const spaceorder = [:ord2, :ord4, :ord6]

const timeorder = [:tord2_1_1, :tord4_5_1, :tord6_7_1]

const plot_times = [0.0, 5.0, 6.75, 7.7, 10.25, 15.5, 20.45]
const norm_times = [1.0, 10.25, 20.0]
const plot_iter = map(x -> Int(floor(x)), plot_times ./ τ)
const norm_iter = map(x -> Int(floor(x)), norm_times ./ τ)
const outputfolder = pwd() * "/scripts/Examples/Results/P3"

function output_name(ind, dim)
    global stateidx
    return outputfolder * "$(dim)output_$(stateidx)_" * string(ind) * "."
end

function test(τv, hv, spord, tord)
    global stateidx

    ind = 1
    function nextind()
        ind += 1
        return ind
    end

    function grabstate(Mem,Mesh,idx,outfile)
        g = Figure(; size = (800, 600), fontsize = 25)
        ax3 = Axis3(g[1, 1])
        ax3.zlabel = ""
        ax3.yticks = WilkinsonTicks(6; k_min = 5)
        ax3.xticks = WilkinsonTicks(6; k_min = 5)
        systemnd!(ax3, Mem, Mesh, idx)
        save(outfile * "c$idx.png", g; px_per_unit = 3)
    end

    PDE = SchrodingerPDEPolynomic((Ω, Ω), (C1, C2, C3, C4), F, N, TF)
    Mesh = PeriodicGrid(backend, PDE, τv, (hv, hv))
    Params = DefaultSolver(backend, 2, (spord, spord), tord)
    Method, Mem = PaulMethod3(backend, PDE, Mesh, Params)
    Stats = initialize_stats(Vector{Float64}, PDE, Mesh, Mem,
                             IterativeLinearSolver(backend), log_freq)
    Problem = SchrodingerProblem(PDE, Params, Method, Mem, Stats)

    tsteps = estimate_timesteps(PDE, Mesh)

    for index in 0:(tsteps - 1)
        step!(Problem)

        if index in plot_iter
            println("Plotting at time: ", index * τv)

            grabstate(Mem,Mesh,1,output_name(nextind(),3))
            grabstate(Mem,Mesh,2,output_name(nextind(),3))
            grabstate(Mem,Mesh,3,output_name(nextind(),3))
            grabstate(Mem,Mesh,4,output_name(nextind(),3))
        end

        if index in norm_iter
            data = dump_backend_memory(Mem)
            tosavejson = Dict()
            tosavejson["time"] = index * τv
            tosavejson["comp1"] = data[:, 1]
            tosavejson["comp2"] = data[:, 2]
            tosavejson["comp3"] = data[:,3]
            tosavejson["comp4"] = data[:,4]
            open(output_name(nextind(),3) * "norm.json", "w") do io
                JSON.print(io, tosavejson)
            end
        end
    end

    px = solveriterations(Stats)
    save(output_name(nextind(),3) * "solviter.png", px; px_per_unit = 3)
    px = solvertime(Stats, u"ms")
    save(output_name(nextind(),3) * "solvtime.png", px; px_per_unit = 3)
    px = executiontime(Stats, u"ms")
    save(output_name(nextind(),3) * "exectime.png", px; px_per_unit = 3)
    px = systemenergy(Stats)
    save(output_name(nextind(),3) * "sysenergy.png", px; px_per_unit = 3)
    px = componentpower(Stats, 1)
    save(output_name(nextind(),3) * "syspowerc1.png", px; px_per_unit = 3)
    px = componentpower(Stats, 2)
    save(output_name(nextind(),3) * "syspowerc2.png", px; px_per_unit = 3)
    px = componentpower(Stats, 3)
    save(output_name(nextind(),3) * "syspowerc3.png", px; px_per_unit = 3)
    px = componentpower(Stats, 4)
    save(output_name(nextind(),3) * "syspowerc4.png", px; px_per_unit = 3)

    store(Stats, output_name(nextind(),3) * "stats.json")

    towrite = pde_parameters(2, τv, hv, spord, tord)
    (f -> open(f, "w") do io
         JSON.print(io, towrite)
     end)(output_name(nextind(),3) * "meta.json")
    stateidx += 1
end

#test(τ,h,spaceorder[1],timeorder[1])
#=
for spord in spaceorder
    for tord in timeorder
        for τv in τvalues
            for hv in hvalues
                try
                    test(τv, hv, spord, tord)
                catch e
                    println("Error: ", e)
                end
            end
        end
    end
end

=#

for (spord,tord,τv,hv) in zip(spaceorder,timeorder,τvalues,hvalues)
    test(τv, hv, spord, tord)

    try
    catch e
        println("Error: ", e)
    end
    break
end