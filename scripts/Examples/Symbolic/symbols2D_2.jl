using Symbolics

@variables x,y,z,r,σ₁₁,σ₂₂,σ₃₃,σ₁₂,σ₂₁,σ₁₃,σ₃₁,σ₂₃,σ₃₂

function calculate_component(nvar)
    F= 0.5*σ₁₁*x^2 + 0.5*σ₂₂*y^2 +0.5*σ₃₃*z^2 + σ₃₁*x*z + σ₃₂*y*z + σ₂₁*x*y

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
println("z: ")
display(calculate_component(z))
