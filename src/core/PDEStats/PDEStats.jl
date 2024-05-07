Base.length(R::RuntimeStats{IntType,FloatType}) where {IntType,FloatType} = R.store_index-1
Base.size(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = size(C.power)
Base.length(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = length(C.power)

function Base.getindex(C::ComponentPower{FloatType,
                                         ArrayType},
                       index) where {FloatType,ArrayType}
    C.power[index]
end

function Base.setindex!(C::ComponentPower{FloatType,ArrayType},
                        value::FloatType,
                        index) where {FloatType,ArrayType}
    C.power[index] = value
end

function initialize_stats(::Type{FloatType}, ::Type{VectorType},
                          ncomponents::IntType, log_freq::IntType,
                          timesteps::IntType,
                          log_solverinfo::Bool = true) where {IntType<:Integer,
                                                              FloatType<:AbstractFloat,
                                                              VectorType<:AbstractArray{FloatType}}
    seq = div(timesteps, log_freq) + 1
    seq_solver = log_solverinfo ? seq : 0
    sys_energy = vundef(VectorType, seq)
    sys_power = ntuple(Returns(ComponentPower(similar(sys_energy))), ncomponents)
    sys_time = vundef(VectorType, seq)
    solver_time = vzeros(VectorType, seq_solver)
    solver_iterations = vzeros(VectorType, seq_solver)

    RuntimeStats(sys_energy, sys_power, sys_time, solver_time, solver_iterations,
                 log_freq, false, zero(IntType), one(IntType))
end

@inline function system_power(Grid, Memory)
    measure(Grid) * vec(sum(abs2.(Memory.current_state); dims = 1))
end

@inline function system_energy(PDE, Grid, Mem, ItStop)
    preA = Mem.preA
    A = Mem.opA
    D = Mem.opD
    components = Mem.current_state
    state_abs2 = Mem.current_state_abs2
    stage1 = Mem.stage1

    @. state_abs2 = abs2(components)

    F = get_field(PDE)

    energy = zero(eltype(components))

    σ_collection = get_σ(PDE)

    for (idx, σ) in enumerate(σ_collection)
        @views comp = components[:, idx]
        b = D * comp

        gmres!(Mem.solver_memory,
               A,
               b;
               restart = true,
               N = preA,
               atol = get_atol(ItStop),
               rtol = get_rtol(ItStop),
               itmax= get_max_iterations(ItStop))

        energy -= σ * dot(comp, Mem.solver_memory.x)
    end
    stage1 .= F(state_abs2)
    energy += sum(stage1; dims = 1)[1]
    real(energy) * measure(Grid)
end

@inline function advance_iteration(R::RuntimeStats{IntType,FloatType}) where {IntType<:Integer,
                                                                              FloatType<:AbstractFloat}
    R.current_iteration += 1
    if mod(R.current_iteration, R.log_frequency) == 0
        R.log_data = true
    end
end

function update_power!(stats::Stats, power,
                       idx::IntType) where {IntType<:Integer,
                                            FloatType<:AbstractFloat,
                                            Stats<:RuntimeStats{IntType,FloatType}}
    for comp in 1:length(stats.system_power)
        stats.system_power[comp][idx] = power[comp]
    end
end

function update_system_energy!(stats::Stats, energy::FloatType,
                               idx::IntType) where {IntType<:Integer,
                                                    FloatType<:AbstractFloat,
                                                    Stats<:RuntimeStats{IntType,
                                                                        FloatType}}
    stats.system_energy[idx] = energy
end

function update_solver_info!(stats::Stats, time,
                             iterations) where {Stats<:RuntimeStats}
    if stats.log_data && length(stats.solver_time) > 0
        idx = stats.store_index
        stats.solver_time[idx] += time
        stats.solver_iterations[idx] += iterations
    end
end

function update_stats!(stats::Stats, time, PDE, Grid, Mem,
                       ItStop) where {Stats<:RuntimeStats}
    if stats.log_data
        idx = stats.store_index
        stats.step_time[idx] = time
        update_power!(stats, system_power(Grid, Mem), idx)
        update_system_energy!(stats, system_energy(PDE, Grid, Mem, ItStop), idx)
        stats.store_index += 1
        stats.log_data = false
    end
    advance_iteration(stats)
end

function startup_stats(stats::Stats, PDE, Grid, Mem,
    ItStop)
    last_index = lastindex(stats.step_time)

    evaluate_ψ(PDE,Grid,Mem)

    update_power!(stats, system_power(Grid, Mem), last_index)
    update_system_energy!(stats, system_energy(PDE, Grid, Mem, ItStop), last_index)
    advance_iteration(stats)
end

function initialize_stats(::Type{VectorType},PDE,Grid,Mem,log_freq::IntType, log_solverinfo::Bool = true)
    tsteps = estimate_timesteps(PDE,Grid)
    initialize_stats(eltype(Grid),VectorType, ncomponents(PDE), log_freq, tsteps, log_solverinfo)
    startup_stats(stats, PDE, Grid, Mem, ItStop)
end

export initialize_stats, update_stats!, update_power!, update_system_energy!,
       system_power, system_energy