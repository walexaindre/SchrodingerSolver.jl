include("PDE_ex2D_1.jl")
using Dates
using Unitful
using JSON
using CUDA
using GLMakie
device!(1)

stateidx = 1

function pde_parameters(ncomp,τv,hv,sord,tord)
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

const log_freq = 20
const backend = CUDABackendF64
const τ = 0.05
const h = 0.16
const TF = 25.0

const hvalues = [h,h/2,h/4,h/8,h/16,h/32]
const τvalues = [τ,τ/2,τ/4,τ/8,τ/16,τ/32]

const spaceorder = [:ord2,:ord4,:ord6,:ord8]

const timeorder = [:tord2_1_1,:tord4_5_1,:tord6_7_1,:tord8_15_1]

const plot_times = [0.0,2.5,5.0,6.75,7.7,10.25,15.5,20.45,24.5]
const norm_times = [1.0,10.25,20.0]
const plot_iter = map(x->floor(x)|>Int,plot_times ./ τ)
const norm_iter = map(x->floor(x)|>Int,norm_times ./ τ)
const outputfolder=pwd()*"/scripts/Examples/Results/"

function output_name(ind)
    global stateidx
    return outputfolder*"$(stateidx)output_"*string(ind)*"."
end

function test(τv,hv,spord,tord)
    global stateidx

    ind = 1
    function nextind()
        ind+=1
        return ind
    end

    PDE = SchrodingerPDEPolynomic((Ω,Ω),(C1,C2),F,N,TF)
    Mesh = PeriodicGrid(backend,PDE,τv,(hv,hv))
    Params = DefaultSolver(backend,2,(spord,spord),tord)
    Method,Mem = PaulMethod3(backend,PDE,Mesh,Params)
    Stats = initialize_stats(Vector{Float64},PDE,Mesh,Mem,IterativeLinearSolver(backend),log_freq)
    Problem = SchrodingerProblem(PDE,Params,Method,Mem,Stats)
    
    tsteps = estimate_timesteps(PDE, Mesh)

    for index in 0:tsteps-1
        step!(Problem)

        if index in plot_iter
            println("Plotting at time: ",index*τv)
            g = Figure(size = (800, 600);fontsize = 25)
            ax3 = Axis3(g[1, 1])
            ax3.zlabel = ""
            ax3.yticks = WilkinsonTicks(6,k_min=5)
            ax3.xticks = WilkinsonTicks(6,k_min=5)
            systemnd!(ax3,Mem,Mesh,1)
            save(output_name(nextind())*"cfirst.png",g,px_per_unit=3)
            gz = Figure(size = (800, 600);fontsize = 25)
            zax3 = Axis3(gz[1, 1])
            zax3.zlabel = ""
            zax3.yticks = WilkinsonTicks(6,k_min=5)
            zax3.xticks = WilkinsonTicks(6,k_min=5)
            systemnd!(zax3,Mem,Mesh,2)
            save(output_name(nextind())*"clast.png",gz,px_per_unit=3)
        end

        if index in norm_iter
            data = dump_backend_memory(Mem)
            tosavejson = Dict()
            tosavejson["time"] = index*τv
            tosavejson["comp1"] = data[:,1]
            tosavejson["comp2"] = data[:,2]
            open(output_name(nextind())*"norm.json","w") do io
                JSON.print(io,tosavejson)
            end
        end
    end




    px = solveriterations(Stats)
    save(output_name(nextind())*"solviter.png",px,px_per_unit=3)
    px = solvertime(Stats,u"ms")
    save(output_name(nextind())*"solvtime.png",px,px_per_unit=3)
    px = executiontime(Stats,u"ms")
    save(output_name(nextind())*"exectime.png",px,px_per_unit=3)
    px = systemenergy(Stats)
    save(output_name(nextind())*"sysenergy.png",px,px_per_unit=3)
    px = componentpower(Stats,1)
    save(output_name(nextind())*"syspowerc1.png",px,px_per_unit=3)
    px = componentpower(Stats,2)
    save(output_name(nextind())*"syspowerc2.png",px,px_per_unit=3)

    store(Stats,output_name(nextind())*"stats.json")




    towrite=pde_parameters(2,τv,hv,spord,tord)
    output_name(nextind())*"meta.json"|>f->open(f,"w") do io
        JSON.print(io,towrite)
    end
    stateidx+=1
end

#test(τ,h,spaceorder[1],timeorder[1])


for τv in τvalues
    for hv in hvalues
        for spord in spaceorder
            for tord in timeorder
                try
                    test(τv,hv,spord,tord)
                catch e
                    println("Error: ",e)
                end
            end
        end
    end
end
