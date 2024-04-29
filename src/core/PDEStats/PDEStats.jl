function initialize_stats(::Type{IntType},::Type{FloatType},::Type{ArrayType},freq,time_steps)
    seq = 1:freq:time_steps

    sys_energy = vzeros(ArrayType{FloatType,1},length(seq))
    sys_power = similar(sys_energy)
    sys_time = similar(sys_energy)
    
    first_v = IntType(first(seq))
    last_v = IntType(last(seq))
    step = IntType(freq)
    first_v=first_v-step

    RuntimeStats(sys_energy,sys_power,sys_time,first_v,last_v,step,first_v)
end

@inline function system_power(Grid, Memory)
    get_measure(Grid) * vec(sum(abs2.(Memory.current_state), dims = 1))
end

@inline function system_energy(Opt,Memory)

end

function update_stats(stats::Stats,memory) where {Stats<:RuntimeStats}


    
    stats.current_iter+=stats.step
end