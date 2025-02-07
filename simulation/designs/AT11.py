# HD1 Sample: Created Sept 2021

# Name
name_part = 'AT11'
scale_factor = 1.73
# Unit-cell Parameters
t = 1.125
angle_base = 0  # Unit-cell orientation
theta_base = angle_base * np.pi/180.0
origin_base = [10.25, -0.01] # Location of UC tile origin (0,0) on the experiment map

# Finite Size Parameters
L_base = 13.765718*scale_factor
ep = np.array([0,  0.5,  0.5,  1],dtype='float64')*L_base
nu = np.array([0.5,  0,  0.,  0.],dtype='float64')*L_base
H_base = nu[0]


# Finite Size Parameters
L_finite = 75
H_finite = 150

