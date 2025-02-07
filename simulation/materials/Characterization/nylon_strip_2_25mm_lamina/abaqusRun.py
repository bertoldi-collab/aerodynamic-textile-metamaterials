# Import Tool
import os
import sys
import numpy as np
import inspect

# Find the correct path name
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
folderdir1 = os.path.dirname(currentdir)
folderdir2 = os.path.dirname(folderdir1)
folderdir3 = os.path.dirname(folderdir2)
parentdir = os.path.dirname(folderdir3)

#########################################################
#                       IMPORTS                         #
#########################################################

# Import classes
execfile(parentdir+'\\abaqus\\abaqus_finite_strips_class.py')


# Import material parameters
























matPar = {
	"PR"                    : 0.3,            # Poisson's ratio of material
    "thickness"             : 0.17,           # textile thickness, mm
    "rho"                   :  9.134e-10,     # textile density
    "alpha_damping"         : 1.5,            # Damping of the material
    "material_orientation"  : 0,              # Orientation of material in degrees
    "DepVar"                : 1,              # Number of DepVar variables
    "matModel"              : 'LAM'           # Defined Material Model
    }

#########################################################
#                         MODEL                         #
#########################################################

# Setup a new model
Mdb()
model = mdb.models['Model-1']
name_part = 'Part-1'


#########################################################
#                       APPLICATION                     #
#########################################################
# Deformation gradient
DeformationF= 20
width = 75
n_strips = 12
t = 2.25
radius = (width-2.25*n_strips)/(2*n_strips)

# sample.print_finite('print_SVG','strip_basic')
cantilever_size = [75,35]
cc = 1
for i,angle in enumerate(np.linspace(0,90,2)):
	for j,length in enumerate(np.linspace(35,75,2)):
		sample = generate_strip(model,length,width,n_strips,t,0,name_part)
		model.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, 
					application=QUASI_STATIC, initialConditions=OFF, initialInc=0.1, maxNumInc=
					500, minInc=5e-12, name='step-'+name_part, nlgeom=ON, nohaf=OFF, previous=
					'Initial', timePeriod=2)
		matPar['material_orientation'] = angle+90
		sample.material_properties(matPar)
		sample.mesh(S3,TRI,t/3) # ElemCode (S3,S4R), ElemShape (TRI,QUAD,QUADDOM)
		sample.generate_RP('RefPoint-0')
		sample.assembly(sample.name_part+'-1')
		sample.amplitude('amp-'+name_part,2,2,20)    # Exp_fact 1 (linear ramp) up to 100 (Exponential)
		sample.boundary_conditions_gravity(-9806.65,'amp-'+name_part,'step-'+name_part)
		sample.create_job('job-strip-angle-'+str(int(angle))+'-'+str(int(length)),nCores = 1)
		sample.run_job()
		sample.post_process_gravity('Abaqus_Bending_'+str(int(cc)),length)
		cc = cc+1

# Unit-Cells
length = 97
sample = generate_strip(model,length,width,n_strips,t,0,name_part)

for i,angle in enumerate(np.linspace(0,90,7)):	
	model.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, 
					application=QUASI_STATIC, initialConditions=OFF, initialInc=0.1, maxNumInc=
					500, minInc=5e-12, name='step-'+name_part, nlgeom=ON, nohaf=OFF, previous=
					'Initial', timePeriod=1)
	disp,force = sample.read_csv('SS1-'+str(i+1)+'.csv')
	DeformationF = disp[-1]
	matPar['material_orientation'] = angle+90
	sample.material_properties(matPar)
	sample.mesh(S3,TRI,t/3) # ElemCode (S3,S4R), ElemShape (TRI,QUAD,QUADDOM)
	sample.generate_RP('RefPoint-0')
	sample.assembly(sample.name_part+'-1')
	sample.amplitude_exp('amp-'+name_part,disp,force)    # Exp_fact 1 (linear ramp) up to 100 (Exponential)
	sample.boundary_conditions(DeformationF,'amp-'+name_part,'step-'+name_part)
	sample.create_job('job-strip-angle-'+str(int(angle)),nCores = 1)
	sample.run_job()
	sample.post_process('Abaqus_Uniaxial_'+str(int(angle)))

