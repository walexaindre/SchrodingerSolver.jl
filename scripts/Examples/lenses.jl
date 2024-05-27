using GLMakie
using LinearAlgebra
using Printf

focallentoangle = (f) -> 2*atand(widthsensor/(2*f))


const d = 1000
const widthsensor = 18.0
const heightsensor = 13.5

const wmod = 17.3
const hmod = 13.0

const focalength = [8.0,12.0,25.0, 45.0, 60.0, 75.0, 150.0, 275.0, 300.0]
const anglesx = focallentoangle.(focalength)
const anglesy = 2*atand.(0.75*tand.(anglesx/2))

width = (ang) -> 2*d*tand(ang/2)
height = (ang) -> 2*d*tand(ang/2)


focalengthw = (ang) -> widthsensor/(2*tand(ang/2))
focalengthh = (ang) -> heightsensor/(2*tand(ang/2))

fig = Figure(size = (800, 800))
ax = Axis(fig[1,1],aspect=DataAspect())
for i in 1:length(anglesx)

    anglex = anglesx[i]
    angley = anglesy[i]
    widthx = width(anglex)
    widthy = width(angley)
    heightx = height(anglex)
    heighty = height(angley)
    medx = widthx/2
    medy = widthy/2
    r = Rect(-medx, -medy, widthx, widthy)

    diagonal = sqrt(wmod^2 + hmod^2)
    aov = 2*atand(diagonal/(2*focalength[i]))
    

    focalengthhx = focalengthw(anglex)


    fx = @sprintf("Focal Length: %.2f",focalengthhx)
    angx = @sprintf("%.2f",anglex)
    angy = @sprintf("%.2f",angley)
    aov = @sprintf("%.2f",aov)
    poly!(ax, r, strokewidth = 2,label = "Angle: $angx ° x $angy ° => $fx mm AOV: $aov °")
    
end

ax.yticks= WilkinsonTicks(6; k_min = 5)
ax.xticks= WilkinsonTicks(6; k_min = 5)


Legend(fig[1,2],fig.content[1])
fig
