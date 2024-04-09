using SchrodingerSolver

rs = get_coefficients(TimeComposition,:tord6_9_1)

get_available(TimeComposition)


for i in 1:7
    if i>4
        println(8-i)
    else
        println(i)
    end
end

for i in 1:7
    println(mod1(8-i,4))
end

