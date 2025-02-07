# Material Parameter List
# Circular Knit with elastomeric fibers
matPar = {
    "mu1"                   : 0.3,
    "alpha1"                : 3.5,
    "D1"                    : 0.0,
    "thickness"             : 0.17,           # textile thickness, mm
    "rho"                   : 9.134e-10,      # textile density
    "alpha_damping"         : 1.5,            # Damping of the material
    "PR"                    : 0.5,            # Poisson's ratio of material
    "matModel"              : 'OGD'           # Defined Material Model
    }

# mdb.models['Model-1'].materials['Knit'].hyperelastic.setValues(n=3, table=((
#     0.27692414,8.67196519349,-0.5338621818360234,7.628592692561689,0.8357421157647154,4.478503339924652,0,0,0), ))

# Material Parameter List
# # Circular Knit with elastomeric fibers from workstation
# matPar = {
#     "mu1"                   : 0.3,#*2.52,
#     "alpha1"                : 3.5,
#     "D1"                    : 0.0,
#     "thickness"             : 0.17,#06746534966163997,           # textile thickness, mm
#     "rho"                   : 9.134e-10,      # textile density
#     "alpha_damping"         : 1.5,            # Damping of the material
#     "PR"                    : 0.5,            # Poisson's ratio of material
#     "matModel"              : 'OGD'           # Defined Material Model
#     }