import numpy as np

WorkerNumber=-1442178722970234786

jobN = "job_-1442178722970234786"

ep=[2.9061138071990587, 1.592078406635753, 11.445005430592593, 7.112396215048806]

nu=[18.18398879075348, 14.310774515809507, 1.7853200560265194, 17.125262221728587]

t=1.125

angle_base=41.55876163198384

L_base=15.037105215484656

H_base =18.18398879075348

disp_load=0.22

cyl_offset=0.2

name_part = 'HD1'
scale_factor = 1.485
# Unit-cell Parameters
t = 1.125
angle_base = 60  # Unit-cell orientation
theta_base = angle_base * np.pi/180.0
origin_base = [10.25, -0.01] # Location of UC tile origin (0,0) on the experiment map

# Finite Size Parameters
L_base = 13.765718*scale_factor
ep = np.array([0.5000,  0.5001,  1.000,  1],dtype='float64')*L_base
nu = np.array([0.8660,  0.8661,  0.0000,  0.000],dtype='float64')*L_base
H_base = nu[0]