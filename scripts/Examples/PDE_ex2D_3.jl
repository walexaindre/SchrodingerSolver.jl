using SchrodingerSolver

const α₁ = 1.0
const α₂ = 1.0
const α₃ = 1.0
const σ₁₁ = 1.0
const σ₂₂ = 2.0
const σ₃₃ = 1.0
const σ₄₄ = 2.0
const σ₁₂ = 2.0/3.0
const σ₂₁ = σ₁₂
const σ₁₃ = 2.0
const σ₃₁ = σ₁₃
const σ₁₄ = 3.0/2.0
const σ₄₁ = σ₁₄
const σ₂₃ = 3.0/2.0
const σ₃₂ = σ₂₃
const σ₂₄ = 1.0
const σ₄₂ = σ₂₄
const σ₃₄ = 2.0/3.0
const σ₄₃ = σ₃₄

const T = 40.0

const Ω = (-8.0,8.0) #Domain

function f₁(x)
    @. @views σ₁₁*x[:,1] + σ₁₂*x[:,2] + σ₁₃*x[:,3] + σ₁₄*x[:,4]
end

function f₂(x)
    @. @views σ₂₁*x[:,1] + σ₂₂*x[:,2] + σ₂₃*x[:,3] + σ₂₄*x[:,4]
end

function f₃(x)
    @. @views σ₃₁*x[:,1] + σ₃₂*x[:,2] + σ₃₃*x[:,3] + σ₃₄*x[:,4]
end

function f₄(x)
    @. @views σ₄₁*x[:,1] + σ₄₂*x[:,2] + σ₄₃*x[:,3] + σ₄₄*x[:,4]
end

#F 0.5*σ₁₁*x^2 + 0.5*σ₂₂*y^2 +0.5*σ₃₃*z^2 + 0.5*σ₄₄*w^2 + σ₃₁*x*z + σ₃₂*y*z + σ₂₁*x*y + σ₄₁*x*w + σ₄₂*y*w + σ₄₃*z*w
function F(x)
    @. @views 0.5*(σ₁₁*x[:,1]^2 + σ₂₂*x[:,2]^2 + σ₃₃*x[:,3]^2 + σ₄₄*x[:,4]^2) + σ₁₂*x[:,1]*x[:,2] + σ₁₃*x[:,1]*x[:,3] + σ₁₄*x[:,1]*x[:,4] + σ₂₃*x[:,2]*x[:,3] + σ₂₄*x[:,2]*x[:,4] + σ₃₄*x[:,3]*x[:,4]
end

function ψ₁(x)
    @. @views 1/sqrt(pi)*exp(-(x[:,1]-2.0)^2-(x[:,2]-2.0)^2)
end

function ψ₂(x)
    @. @views 1/sqrt(pi)*exp(-(x[:,1]+2.0)^2-(x[:,2]+2.0)^2)
end

function ψ₃(x)
    @. @views 1/sqrt(pi)*exp(-(x[:,1]-2.0)^2-(x[:,2]+2.0)^2)
end

function ψ₄(x)
    @. @views 1/sqrt(pi)*exp(-(x[:,1]+2.0)^2-(x[:,2]-2.0)^2)
end

C1 = SchrodingerPDEComponent(α₁,f₁,ψ₁)
C2 = SchrodingerPDEComponent(α₂,f₂,ψ₂)
C3 = SchrodingerPDEComponent(α₃,f₃,ψ₃)
C4 = SchrodingerPDEComponent(α₃,f₄,ψ₄)

#x:
#-0.5r*σ₁₁ - w*σ₄₁ - 0.5x*σ₁₁ - y*σ₂₁ - z*σ₃₁
#y:
#-0.5r*σ₂₂ - w*σ₄₂ - x*σ₂₁ - 0.5y*σ₂₂ - z*σ₃₂
#z:
#-0.5r*σ₃₃ - w*σ₄₃ - x*σ₃₁ - y*σ₃₂ - 0.5z*σ₃₃
#w:
#-0.5r*σ₄₄ - 0.5w*σ₄₄ - x*σ₄₁ - y*σ₄₂ - z*σ₄₃
function N(prev,next,idx)
    if idx == 1
        return @. @views -0.5*σ₁₁*next-0.5*σ₁₁*prev[:,1]-σ₁₂*prev[:,2]-σ₁₃*prev[:,3]-σ₁₄*prev[:,4]
    elseif idx == 2
        return @. @views -0.5*σ₂₂*next-σ₂₁*prev[:,1]-0.5*σ₂₂*prev[:,2]-σ₂₃*prev[:,3]-σ₂₄*prev[:,4]
    elseif idx == 3
        return @. @views -0.5*σ₃₃*next-σ₃₁*prev[:,1]-σ₃₂*prev[:,2]-0.5*σ₃₃*prev[:,3]-σ₃₄*prev[:,4]
    else
        return @. @views -0.5*σ₄₄*next-0.5*σ₄₄*prev[:,4]-σ₄₁*prev[:,1]-σ₄₂*prev[:,2]-σ₄₃*prev[:,3]
    end

end

PDE = SchrodingerPDEPolynomic((Ω,Ω),(C1,C2,C3,C4),F,N,T)