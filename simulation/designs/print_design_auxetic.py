# Import Tool
import os
import sys
import inspect
import numpy as np

# Import tooling
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)

#########################################################
#                       IMPORTS                         #
#########################################################

# Import classes
execfile(parentdir+'\\abaqus\\abaqus_unit_cell_class.py')
execfile(parentdir+'\\abaqus\\abaqus_finite_size_class.py')

# Import material parameters
execfile(parentdir+'\\materials\\LAM_nylon_strips_2-25.py')
woven_matPar = matPar
execfile(parentdir+'\\materials\\OGD_knit3.py')
knit_matPar = matPar

# Import Design
name_part ='AT13'
execfile(parentdir+'\\designs\\{}.py'.format(name_part))

#########################################################
#                         MODEL                         #
#########################################################

# Setup a new model
Mdb()
model = mdb.models['Model-1']
name_part = name_part
name_assem = name_part+'-finite'
name_step = 'Step-1'
name_print = name_part+'-print'

#########################################################
#                       SIMULATION                      #
#########################################################

model.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, 
                application=QUASI_STATIC, initialConditions=OFF, initialInc=0.1, maxNumInc=
                10000, minInc=5e-12, name=name_step, nlgeom=ON, nohaf=OFF, previous=
                'Initial', timePeriod=1)

#########################################################
#                       APPLICATION                     #
#########################################################
dimension=THREE_D

# Sample Design Parameters
# Offset parameter: To center unit-cells within a given sample
origin_base=[0.1,0]

# Width and Height of the sample
W_sample = 75
H_sample = 150

# Legacy Parameters
L_finite = W_sample
H_finite = H_sample

#Generate Unit-Cell
sampleUC = unit_cell(model,name_part)
x_offset = height/(2*np.tan(alpha))
area = np.abs(np.cross(np.array([width*2-x_offset*2,0]),np.array([0,height])))
sampleUC.generate_parallelogram(0,width*2-x_offset*2,height,0,dimension=THREE_D)

# Partition surface with re-entrant design
sampleUC.RE_partition(width,height,alpha,t)

# Apply material properties to the design (allows for easy removal of knit)
sampleUC.material_properties(matPar,sampleUC.model.parts[name_part].faces,matName='Woven')
sampleUC.material_properties(knit_matPar,sampleUC.knitFaces,matName='Knit')
sampleUC.model.parts[name_part].RemoveFaces(deleteCells=False, faceList=sampleUC.knitFaces)

UC_list = [sampleUC]
UC_boundingbox = [[0,0,0,W_sample,H_sample,origin_base[0], origin_base[1],0]] # [number-unit-cell,xbox1,ybox1,xbox2,ybox2,xorigin,yorigin,zorigin]
sampleFinite = generate_finite(model,UC_list,UC_boundingbox)
sampleFinite.geometry_uniaxial(W_sample,H_sample,name_assem,dimension=dimension)
sampleFinite.print_finite(name_part+'_temp',name_part+'_outer',instron_tabs=True,hook_tabs=False)

