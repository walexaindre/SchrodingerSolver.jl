abstract type SchrodingerPDE end

struct SchrodingerPDEComponent{Tv,Fn,Hn}
    σ::Tv #Dispersion coefficient
    f::Fn 
    ψ::Hn #Initial condition
end



struct SchrodingerPDEPolynomic{N,M,Tv,Comp<:SchrodingerPDEComponent{Tv}}
    components::NTuple{M,Comp}
end

export SchrodingerPDE, SchrodingerPDEComponent, SchrodingerPDEPolynomic