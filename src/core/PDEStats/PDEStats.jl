function initialize_stats(::Type{IntType},::Type{FloatType},::Type{VectorType},ncomponents,freq,seq) where {IntType<:Integer,FloatType<:AbstractFloat,VectorType<:AbstractArray{FloatType}}

    sys_energy = vzeros(VectorType,length(seq))
    sys_power = ntuple(Returns(ComponentPower(similar(sys_energy))),ncomponents)
    sys_time = similar(sys_energy)
    
    first_v = IntType(first(seq))
    last_v = IntType(last(seq))
    step = IntType(freq)
    first_v=first_v-step

    RuntimeStats(sys_energy,sys_power,sys_time,first_v,last_v,step,first_v)
end

@inline function system_power(Grid, Memory)
    measure(Grid) * vec(sum(abs2.(Memory.current_state), dims = 1))
end

@inline function system_energy(PDE,Mem,ItStop)
    preA = Mem.preA
    A = Mem.opA
    D = Mem.opD
    σ_collection =  get_σ(PDE)

    for (comp,σ) in enumerate(σ_collection)
        @views comp = components[:, idx]
        b = D * comp

        gmres!(Memory.solver_memory,
            A,
            b,
            restart = true,
            N = preA,
            atol = ItStop.atol,
            rtol = ItStop.rtol)

        energy -= σ * dot(comp, Memory.solver_memory.x)
    end

end

function update_power!(stats::Stats,power,idx::IntType) where {IntType<:Integer,FloatType<:AbstractFloat,Stats<:RuntimeStats{IntType,FloatType}}
    for comp in 1:length(stats.system_power)
        stats.system_power[comp][idx] = power[comp]
    end
end

function update_system_energy!(stats::Stats,energy::FloatType,idx::IntType) where {IntType<:Integer,FloatType<:AbstractFloat,Stats<:RuntimeStats{IntType,FloatType}}
    stats.system_energy[idx] = energy
end

function update_stats!(stats::Stats,time,PDE,Grid,Mem,ItStop) where {Stats<:RuntimeStats}
    stats.current_iter+=stats.step

    idx  = div(stats.current_iter-stats.first,stats.step)+1

    stats.step_time[idx]=time
    update_power!(stats,system_power(Grid,Mem),idx)
    update_system_energy!(stats,system_energy(PDE,Mem,ItStop),idx)
end