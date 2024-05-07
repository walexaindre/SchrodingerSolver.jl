using Symbolics

@variables x,y,w,z,r,σ₁₁,σ₂₂,σ₃₃,σ₄₄,σ₁₂,σ₂₁,σ₁₃,σ₃₁,σ₁₄,σ₄₁,σ₂₃,σ₃₂,σ₂₄,σ₄₂,σ₃₄,σ₄₃

function calculate_component(nvar)
    F= 0.5*σ₁₁*x^2 + 0.5*σ₂₂*y^2 +0.5*σ₃₃*z^2 + 0.5*σ₄₄*w^2 + σ₃₁*x*z + σ₃₂*y*z + σ₂₁*x*y + σ₄₁*x*w + σ₄₂*y*w + σ₄₃*z*w

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
println("w: ")
display(calculate_component(w))