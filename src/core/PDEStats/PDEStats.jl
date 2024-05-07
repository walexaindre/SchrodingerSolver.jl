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

function initialize_stats(::Type{IntType}, ::Type{FloatType}, ::Type{VectorType},
                          ncomponents, freq,
                          seq) where {IntType<:Integer,FloatType<:AbstractFloat,
                                      VectorType<:AbstractArray{FloatType}}
    sys_energy = vzeros(VectorType, length(seq))
    sys_power = ntuple(Returns(ComponentPower(similar(sys_energy))), ncomponents)
    sys_time = similar(sys_energy)

    first_v = IntType(first(seq))
    last_v = IntType(last(seq))
    step = IntType(freq)
    first_v = first_v - step

    RuntimeStats(sys_energy, sys_power, sys_time, first_v, last_v, step, first_v)
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
               atol = ItStop.atol,
               rtol = ItStop.rtol)

        energy -= σ * dot(comp, Mem.solver_memory.x)
    end
    stage1 .= F(state_abs2)
    energy+=sum(stage1,dims=1)[1]
    real(energy)*measure(Grid)
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

function update_stats!(stats::Stats, time, PDE, Grid, Mem,
                       ItStop) where {Stats<:RuntimeStats}
    stats.current_iter += stats.step
    idx = div(stats.current_iter - stats.first, stats.step)

    stats.step_time[idx] = time
    update_power!(stats, system_power(Grid, Mem), idx)
    update_system_energy!(stats, system_energy(PDE, Grid,Mem, ItStop), idx)
end

export initialize_stats, update_stats!, update_power!, update_system_energy!,
       system_power, system_energy