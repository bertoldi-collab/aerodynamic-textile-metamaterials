# Material Parameter List
# Nylon Heat Bondable Fabric with 1.5mm strips
matPar = {
    "C10"                   : 12,           # Isotropic stiffness
    "D0"                    : 0.001,          # Compressibility
    "thickness"             : 0.17,           # textile thickness, mm
    "rho"                   :  9.134e-10,     # textile density
    "alpha_damping"         : 1.5,            # Damping of the material
    "material_orientation"  : 0,              # Orientation of material in degrees
    "PR"                    : 1/(2*(12*0.001+1)),            # Poisson's ratio of material
    "matModel"              : 'NH'           # Defined Material Model
    }