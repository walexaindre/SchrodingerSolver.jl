using SchrodingerSolver
using GLMakie
const α = 1.5
const β = -0.5
const γ = α

#Always define functions in a vectorial way
function f1(x)
    @views  α*x[:,1]+β*x[:,2]
end

function f2(x)
    @views β*x[:,1]+γ*x[:,2]
end

function psi1(x)
    @views 1/sqrt(pi)*exp.(-x[:,1].^2-x[:,2].^2)
end

function psi2(x)
    1im*psi1(x)
end

C1 = SchrodingerPDEComponent(0.5,f1,psi1)
C2 = SchrodingerPDEComponent(0.5,f2,psi1)

function F(x)
    @views 0.5*α*x[:,1].^2+0.5*γ*x[:,2].^2+β*x[:,1].*x[:,2]
end

function N(prev,next,idx)
    if idx == 1
        return -0.5*α*next-0.5*α*prev[:,1]-β*prev[:,2]
    else
        return -0.5*γ*next-β*prev[:,1]-0.5*γ*prev[:,2]
    end
end

PDE = SchrodingerPDEPolynomic(((-8.0,8.0),(-8.0,8.0)),(C1,C2),F,N,40.0)

Mesh = PeriodicGrid(CPUBackendF64,PDE,0.1,(100,100))

Params = DefaultSolver(CPUBackendF64,2)

Method,Mem = PaulMethod1(CPUBackendF64,PDE,Mesh,Params)

evaluate_ψ(PDE,Mesh,Mem)

Stats = initialize_stats(Int64,Float64,Vector{Float64},ncomponents(PDE),4,1:3000)

update_stats!(Stats,0.0,PDE,Mesh,Mem,Method.linear_solve_params)

Problem = SchrodingerProblem(PDE,DefaultSolver(CPUBackendF64,2),Method,Mem,Stats)

step!(Problem)

for i in 1:300
    step!(Problem)
end


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
surface!(ax3,Mesh[:,1],Mesh[:,2],abs2.(Mem.current_state[:,1]))
g