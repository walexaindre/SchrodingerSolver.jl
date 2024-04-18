abstract type SchrodingerPDE end

"PDE component"
struct SchrodingerPDEComponent{Tv,Fn,Hn}
    σ::Tv #Dispersion coefficient
    f::Fn 
    ψ::Hn #Initial condition
end


"Generic handler for non polynomic potentials"
struct SchrodingerPDENonPolynomic{N,Tv,MComp<:SchrodingerPDEComponent{Tv},Potential} <: SchrodingerPDE
    boundaries::NTuple{N,NTuple{2,Tv}}
    components::MComp
    F::Potential
    T::Tv
end

"Optimized structure for polynomical potentials"
struct SchrodingerPDEPolynomic{N,Tv,MComp,Potential,Optimized} <: SchrodingerPDE
    boundaries::NTuple{N,NTuple{2,Tv}}
    components::MComp
    F::Potential
    N::Optimized
    T::Tv
end

export SchrodingerPDE, SchrodingerPDEComponent, SchrodingerPDEPolynomic, SchrodingerPDENonPolynomic


