using Symbolics

@variables x,y,z,w,r,α,e,k,β,γ

function FieldF(A)
    0.5 * (A[:, 1] .^ 2 + A[:, 2] .^ 2 + A[:, 3] .^ 2) +
    (2.0 / 3.0) * (A[:, 1] .* A[:, 2] + A[:, 2] .* A[:, 3]) + A[:, 1] .* A[:, 3]
end





function calculate_component(nvar)
    F= 0.5*α*x^2+ 0.5*γ*y^2+β*x*y
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
