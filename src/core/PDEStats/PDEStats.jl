Base.length(R::RuntimeStats{IntType,FloatType}) where {IntType,FloatType} = R.store_index -
                                                                            1
Base.size(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = size(C.power)
Base.length(C::ComponentPower{FloatType,ArrayType}) where {FloatType,ArrayType} = length(C.power)

@inline islocked(R::RuntimeStats) = R.locked

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

function initialize_stats(::Type{VectorType},
                          ncomponents::IntType, log_freq::IntType, τ::FloatType,
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
                 log_freq, τ, false, false, zero(IntType), one(IntType))
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
               itmax = get_max_iterations(ItStop))

        energy -= σ * dot(comp, Mem.solver_memory.x)
    end
    stage1 .= F(state_abs2)
    vecenergy = sum(stage1; dims = 1)|>Vector
    energy += vecenergy[1]
    real(energy) * measure(Grid)
end

@inline function advance_iteration(R::RuntimeStats{IntType,FloatType}) where {IntType<:Integer,
                                                                              FloatType<:AbstractFloat}
    if !islocked(R)
        R.current_iteration += 1
        if mod(R.current_iteration, R.log_frequency) == 0
            R.log_data = true
        end
    else
        error("Stats are locked")
    end
end

function update_power!(stats::Stats, power,
                       idx::IntType) where {IntType<:Integer,
                                            FloatType<:AbstractFloat,
                                            Stats<:RuntimeStats{IntType,FloatType}}
    power = power|>Vector

    if !islocked(stats)
        for comp in 1:length(stats.system_power)
            stats.system_power[comp][idx] = power[comp]
        end
    else
        error("Stats are locked")
    end
end

function update_system_energy!(stats::Stats, energy::FloatType,
                               idx::IntType) where {IntType<:Integer,
                                                    FloatType<:AbstractFloat,
                                                    Stats<:RuntimeStats{IntType,
                                                                        FloatType}}
    if !islocked(stats)
        stats.system_energy[idx] = energy
    else
        error("Stats are locked")
    end
end

function update_solver_info!(stats::Stats, time,
                             iterations) where {Stats<:RuntimeStats}
    if !islocked(stats)
        if stats.log_data && length(stats.solver_time) > 0
            idx = stats.store_index
            stats.solver_time[idx] += time
            stats.solver_iterations[idx] += iterations
        end
    else
        error("Stats are locked")
    end
end

function update_stats!(stats::Stats, time, PDE, Grid, Mem,
                       ItStop) where {Stats<:RuntimeStats}
    if !islocked(stats)
        if stats.log_data
            idx = stats.store_index
            stats.step_time[idx] = time
            update_power!(stats, system_power(Grid, Mem), idx)
            update_system_energy!(stats, system_energy(PDE, Grid, Mem, ItStop), idx)
            stats.store_index += 1
            stats.log_data = false
        end
        advance_iteration(stats)
    else
        error("Stats are locked")
    end
end

function startup_stats(stats::Stats, PDE, Grid, Mem,
                       ItStop) where {Stats<:RuntimeStats}
    last_index = lastindex(stats.step_time)

    evaluate_ψ(PDE, Grid, Mem)

    update_power!(stats, system_power(Grid, Mem), last_index)
    update_system_energy!(stats, system_energy(PDE, Grid, Mem, ItStop), last_index)
    advance_iteration(stats)
end

function initialize_stats(::Type{VectorType}, PDE, Grid::PerGrid, Mem, ItStop,
                          log_freq::IntType,
                          log_solverinfo::Bool = true) where {IntType,FloatType,
                                                              PerGrid<:PeriodicGrid{IntType,
                                                                                    FloatType},
                                                              VectorType<:AbstractArray}
    tsteps = estimate_timesteps(PDE, Grid)
    τ = Grid.τ
    stats = initialize_stats(VectorType, ncomponents(PDE), log_freq,τ,
                             tsteps, log_solverinfo)
    startup_stats(stats, PDE, Grid, Mem, ItStop)
    stats
end

function serialize(array, len)
    output = vzeros(typeof(array), len + 1)
    output[2:(len + 1)] .= array[1:len]
    output[1] = array[end]
    output
end

function deserialize(::Type{VectorType}, array) where {VectorType<:AbstractVector}
    sz = length(array)
    output = vundef(VectorType, sz)
    output[1:(sz - 1)] .= array[2:end]
    output[end] = array[1]
    output
end

function store(stats::Stats, path) where {Stats<:RuntimeStats}
    max_elem = stats.store_index - 1

    serialized_energy = serialize(stats.system_energy, max_elem)
    serialized_power = map(x -> serialize(x.power, max_elem), stats.system_power)
    serialized_time = serialize(stats.step_time, max_elem)

    solver_time = serialize(stats.solver_time, max_elem)
    solver_iterations = serialize(stats.solver_iterations, max_elem)

    logfreq = stats.log_frequency
    current_iter = stats.current_iteration
    store_index = stats.store_index
    τ = stats.τ


    parameters = Dict()
    parameters["log_frequency"] = logfreq
    parameters["current_iteration"] = current_iter
    parameters["τ"] = τ
    parameters["store_index"] = store_index
    parameters["date"] = Dates.now(UTC)
    parameters["struct"] = "RuntimeStats"
    parameters["float_type"] = string(eltype(solver_time))
    parameters["int_type"] = string(eltype(max_elem))
    tojson = Dict()

    tojson["system_energy"] = serialized_energy
    tojson["system_power"] = serialized_power
    tojson["step_time"] = serialized_time
    tojson["solver_time"] = solver_time
    tojson["solver_iterations"] = solver_iterations
    tojson["metadata"] = parameters

    open(path, "w") do io
        JSON.print(io, tojson, 4)
    end
end

function deserialize(::Type{VectorType}, ::Type{IntType},
                     path::String) where {VectorType<:AbstractVector,
                                          IntType<:Integer}
    tojson = JSON.parsefile(path)
    system_energy = deserialize(VectorType, tojson["system_energy"])
    system_power = Tuple(map(x -> ComponentPower(deserialize(VectorType, x)),
                             tojson["system_power"]))
    step_time = deserialize(VectorType, tojson["step_time"])
    solver_time = deserialize(VectorType, tojson["solver_time"])
    solver_iterations = deserialize(VectorType, tojson["solver_iterations"])

    parameters = tojson["parameters"]

    logfreq = IntType(parameters["log_frequency"])
    current_iter = IntType(parameters["current_iteration"])
    store_index = IntType(parameters["store_index"])
    τ = FloatType(parameters["τ"])

    RuntimeStats(system_energy, system_power, step_time, solver_time,
                 solver_iterations, logfreq, τ, true, false, current_iter, store_index)
end

function calculate_diff_system_energy(stats::Stats) where {IntType,FloatType,
                                                           ArrayType,
                                                           Stats<:RuntimeStats{IntType,
                                                                               FloatType,
                                                                               ArrayType}}
    max_elem = stats.store_index - 1

    startup_energy = stats.system_energy[end]

    output = vundef(ArrayType, max_elem)

    @inbounds @simd for i in 1:max_elem
        output[i] = abs(stats.system_energy[i] - startup_energy)
    end
    output
end

export initialize_stats, update_stats!, update_power!, update_system_energy!,
       system_power, system_energy, store, deserialize, serialize