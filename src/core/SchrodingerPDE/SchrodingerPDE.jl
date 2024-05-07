"Number of dimensions for the Schrodinger PDE"
@inline Base.ndims(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = N

"Number of dimensions for the Schrodinger PDE"
@inline Base.ndims(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = N

"Number of components for the Schrodinger PDE"
#@inline Base.length(SPDE::SchrodingerPDEPolynomic{N,M,Tv,Comp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = M

"Number of components for the Schrodinger PDE"
#@inline Base.length(SPDE::SchrodingerPDENonPolynomic{N,M,Tv,Comp,Potential}) where {N,Tv,MComp,Potential} = M

"Type of the elements for the Schrodinger PDE"
@inline Base.eltype(::Type{SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}}) where {N,Tv,MComp,Potential,Optimized} = Tv

"Type of the elements for the Schrodinger PDE"
@inline Base.eltype(::Type{SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}}) where {N,Tv,MComp,Potential} = Tv

"Boundary start and end point at index dimension"
@inline get_boundary(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.boundaries[index]

"Boundary start and end point at index dimension"
@inline get_boundary(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.boundaries[index]

"Obtain the component of the Schrodinger PDE at index"
@inline get_component(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index]

"Obtain the component of the Schrodinger PDE at index"
@inline get_component(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index]

"Dispersion coefficient for the Schrodinger PDE component at index"
@inline get_σ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].σ

@inline get_σ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = ntuple(idx -> SPDE.components[idx].σ,
                                                                                                                                         N)

"Dispersion coefficient for the Schrodinger PDE component at index"
@inline get_σ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].σ

@inline get_σ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = ntuple(idx -> SPDE.components[idx].σ,
                                                                                                            N)

"Function for the Schrodinger PDE component"
@inline get_f(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].f

"Function for the Schrodinger PDE component"
@inline get_f(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].f

"Initial condition for the Schrodinger PDE component"
@inline get_ψ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}, index::Int) where {N,Tv,MComp,Potential,Optimized} = SPDE.components[index].ψ

"Initial condition for the Schrodinger PDE component"
@inline get_ψ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}, index::Int) where {N,Tv,MComp,Potential} = SPDE.components[index].ψ

"Potential for the Schrodinger PDE"
@inline get_field(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.F

"Potential for the Schrodinger PDE"
@inline get_field(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = SPDE.F

@inline get_optimized(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.N

@inline get_optimized(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = nothing

"Finish time for the Schrodinger PDE"
@inline get_time_boundary(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized}) where {N,Tv,MComp,Potential,Optimized} = SPDE.T

"Finish time for the Schrodinger PDE"
@inline get_time_boundary(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential}) where {N,Tv,MComp,Potential} = SPDE.T

@inline function evaluate_ψ(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized},P::PGrid,Memory) where {N,Tv,MComp,Potential,Optimized,PGrid<:PeriodicGrid}
    cstate = Memory.current_state
    points = collect_points(P)|>typeof(cstate)
    for i in 1:ncomponents(SPDE)
        ψ_ = view(cstate,:,i)
        func = get_ψ(SPDE,i)
        ψ_ .= func(points)
    end
    nothing
end

@inline function evaluate_ψ(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential},P::PGrid,Memory) where {N,Tv,MComp,Potential,PGrid<:PeriodicGrid}
    cstate = Memory.current_state
    points = collect_points(P)|>typeof(cstate)
    for i in 1:ncomponents(SPDE)
        ψ_ = view(cstate,:,i)
        func = get_ψ(SPDE,i)
        ψ_ .= func(points)
    end
    nothing
end

@inline ncomponents(SPDE::PDEeq) where {PDEeq<:SchrodingerPDE} = length(SPDE.components)

@inline function estimate_timesteps(SPDE::SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized},P::PGrid) where {N,Tv,MComp,Potential,Optimized,PGrid<:PeriodicGrid} 
    τ = P.τ
    T = get_time_boundary(SPDE)
    rank = 0:τ:T
    return length(rank)-1
end

@inline function estimate_timesteps(SPDE::SchrodingerPDENonPolynomic{N,Tv,MComp,Potential},P::PGrid) where {N,Tv,MComp,Potential,PGrid<:PeriodicGrid} 
    τ = P.τ
    T = get_time_boundary(SPDE)
    rank = 0:τ:T
    return length(rank)-1
end

export get_boundary,get_component,get_σ,get_f,get_ψ,get_field,get_time_boundary,ncomponents,get_optimized,evaluate_ψ