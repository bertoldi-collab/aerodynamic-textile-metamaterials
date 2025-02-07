# HD1 Sample: Created Sept 2021

# Name

name_part = 'AT7'
scale_factor = 1.754
# Unit-cell Parameters
ep = np.array([2.3541, 4, 6.8829, 10.88],dtype='float64')*scale_factor
nu = np.array([ 6.4679, 4.12, 0, 4.12],dtype='float64')*scale_factor
t = 1.125
angle_base = 55  # Unit-cell orientation
theta_base = angle_base * np.pi/180.0
origin_base = [2, 0] # Location of UC tile origin (0,0) on the experiment map

# Finite Size Parameters
L_base = 13.765718*scale_factor
H_base = nu[0]

# Finite Size Parameters
L_finite = 75
H_finite = 150

