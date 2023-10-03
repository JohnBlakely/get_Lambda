include("src/recfast.jl")
using .Recfast
p = Recfast.Params()
#Recfast.set_input_variables!(p)
Recfast.H_z_loc(p, 1)
