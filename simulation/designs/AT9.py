# HD1 Sample: Created Sept 2021

# Name

name_part = 'AT9'
scale_factor = 1.8
# Unit-cell Parameters
t = 1.125
angle_base = 30  # Unit-cell orientation
theta_base = angle_base * np.pi/180.0
origin_base = [7, 0] # Location of UC tile origin (0,0) on the experiment map

# Finite Size Parameters
L_base = 13.765718*scale_factor
ep = np.array([-0.25,  0.25,  0.5,  0.75],dtype='float64')*L_base
nu = np.array([0.443,  0.1443,  0.0000,  0.1443],dtype='float64')*L_base
H_base = nu[0]

# Finite Size Parameters
L_finite = 75
H_finite = 150

