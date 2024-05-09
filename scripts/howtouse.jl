using SchrodingerSolver
using GLMakie
using CUDA
device!(1)

const α₁ = 1.0
const α₂ = 1.0
const σ₁₁ = 1.0
const σ₂₂ = 1.0
const σ₁₂ = 2.0/3.0
const σ₂₁ = σ₁₂

const T = 40.0

const Ω = (-8.0,8.0) #Domain

function f₁(x)
    @. @views σ₁₁*x[:,1] + σ₁₂*x[:,2]
end

function f₂(x)
    @. @views σ₂₁*x[:,1] + σ₂₂*x[:,2]
end

function F(x)
    @. @views 0.5*(σ₁₁*x[:,1]^2 + σ₂₂*x[:,2]^2) + σ₁₂*x[:,1]*x[:,2]
end

function ψ₁(x)
    @. @views 1/sqrt(pi)*exp(-(x[:,1]-2.0)^2-(x[:,2]-2.0)^2)
end

function ψ₂(x)
    @. @views 1/sqrt(pi)*exp(-(x[:,1]+2.0)^2-(x[:,2]+2.0)^2)
end

#x: 
#-0.5r*σ₁₁ - 0.5x*σ₁₁ - y*σ₁₂
#y: 
#-0.5r*σ₂₂ - x*σ₁₂ - 0.5y*σ₂₂
function N(prev,next,idx)
    if idx == 1
        return @. @views -0.5*σ₁₁*next-0.5*σ₁₁*prev[:,1]-σ₁₂*prev[:,2]
    else
        return @. @views -0.5*σ₂₂*next-σ₁₂*prev[:,1]-0.5*σ₂₂*prev[:,2]
    end
end

C1 = SchrodingerPDEComponent(α₁,f₁,ψ₁)
C2 = SchrodingerPDEComponent(α₂,f₂,ψ₂)

PDE = SchrodingerPDEPolynomic((Ω,Ω),(C1,C2),F,N,T)

Mesh = PeriodicGrid(CPUBackendF64,PDE,0.01,(1000,1000))

Params = DefaultSolver(CPUBackendF64,2,(:ord4,:ord4))



Method,Mem = PaulMethod1(CPUBackendF64,PDE,Mesh,Params)
Method,Mem = PaulMethod2(CPUBackendF64,PDE,Mesh,Params)
Method,Mem = PaulMethod3(CUDABackendF64,PDE,Mesh,Params)

Stats = initialize_stats(Vector{Float64},PDE,Mesh,Mem,IterativeLinearSolver(CPUBackendF64),1)

Problem = SchrodingerProblem(PDE,Params,Method,Mem,Stats)

step!(Problem)

for i in 1:900
    step!(Problem)
end

save("executiontimeGPU.png",fig,px_per_unit=3)

#sA = Mem.current_state_abs2

#szA = size(sA,1)

#plot(1:szA,real(sA[:,1]))

f = Figure(size = (800, 600))
ax = Axis(f[1, 1],yscale = log10)
plot!(ax,3:300,abs.(Stats.system_power[1][2].-Stats.system_power[1][3:300]),color = :red)
plot!(ax,3:300,abs.(Stats.system_power[2][2].-Stats.system_power[2][3:300]),color = :blue)
plot!(ax,3:300,abs.(Stats.system_energy[2].-Stats.system_energy[3:300]),color = :purple)
plot!(ax,1:10000,abs2.(Mem.current_state[:,2]),color = :green)
f


g = Figure(size = (800, 600))
ax3 = Axis3(g[1, 1])
systemnd!(ax3,Mem,Mesh)


surface!(ax3,Mesh[:,1],Mesh[:,2],abs2.(Mem.current_state[:,1]))
g


Makie.record(g,"file.mp4",1:80,framerate = 5) do _
    empty!(ax3)
    CUDA.@time for i in 1:5
        step!(Problem)
    end
    systemnd!(ax3,Mem,Mesh)
end


### Tuning figure parameters

f = solvertime(Stats,u"ms")
save("solvertime.png",f)

f.figure[:size] = (800, 600)












BI,BJ,BV = findnz(opB)
brs = sparse(BI|>Array,BJ|>Array,BV|>Array)
lop =  drop(brs,Method.Mesh|>PeriodicAbstractMesh,0.1)
brf = sparse(preBI,preBJ,preBV)