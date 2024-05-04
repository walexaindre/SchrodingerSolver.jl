A=(-2,-1,0,1,2)
ta = A[findall(x->x>=0,A)]
@nexpr 
count(A)
for a in ta
    if a in A && -a in A
        println("Gotcha $a ")
    end
end

reduce(min,A)

