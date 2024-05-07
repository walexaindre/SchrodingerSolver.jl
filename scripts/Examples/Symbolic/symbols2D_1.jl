using Symbolics

@variables x,y,r,σ₁₁,σ₂₂,σ₁₂,σ₂₁


function calculate_component(nvar)
    F= 0.5*σ₁₁*x^2 + 0.5*σ₂₂*y^2 + σ₁₂*x*y

    F=-F
    Fnew = substitute(F,Dict([nvar=>r]))
    Div = (Fnew-F)/(r-nvar) 
    R = simplify(Div,expand=true,threaded=true)
    simplify(R)
end

println("x: ")
display(calculate_component(x))
println("y: ")
display(calculate_component(y))