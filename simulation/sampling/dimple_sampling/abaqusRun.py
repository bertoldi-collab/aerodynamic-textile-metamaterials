# Import Tool
import os
import sys
import inspect
import numpy as np

# Import tooling
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
folderdir1 = os.path.dirname(currentdir)
folderdir2 = os.path.dirname(folderdir1)
parentdir = os.path.dirname(folderdir2)

#########################################################
#                       IMPORTS                         #
#########################################################
for root, dirs, files in os.walk('.'):
    for filename in files:
        if filename.endswith('.lck'):
            filepath = os.path.join(root, filename)
            os.remove(filepath)
            #print(f"Deleted {filepath}")
# Import classes
execfile(parentdir+'\\abaqus\\abaqus_unit_cell_class.py')

# Import material parameters
execfile(parentdir+'\\materials\\LAM_nylon_strips_2-25.py')
woven_matPar = matPar
execfile(parentdir+'\\materials\\OGD_knit3.py')
knit_matPar = matPar

#########################################################
#                         MODEL                         #
#########################################################
# Exceute parameter code

# Setup a new model
Mdb()
model = mdb.models['Model-1']
name_part = 'sample'

#########################################################
#                       SIMULATION                      #
#########################################################

model.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, 
                application=QUASI_STATIC, initialConditions=OFF, initialInc=0.1, maxNumInc=
                1000000, minInc=5e-12, name='step-'+name_part, nlgeom=ON, nohaf=OFF, previous=
                'Initial', timePeriod=1)

#########################################################
#                       APPLICATION                     #
#########################################################
# Deformation gradient
DeformationF=[(0, 0),
        (0, disp_load)]

# Unit-Cell and RP Generation
sampleUC = unit_cell(model,name_part)
x_offset = height/(2*np.tan(alpha))
area = np.abs(np.cross(np.array([width*2-x_offset*2,0]),np.array([0,height])))
sampleUC.generate_parallelogram(0,width*2-x_offset*2,height,0,dimension=THREE_D)
sampleUC.generate_RP('RefPoint-1','RefPoint-2')

# Partition surface with re-entrant design
sampleUC.RE_dimple_tracker(width,height,alpha,t)

sampleUC.RE_partition(width,height,alpha,t)


V = 0
sampleUC.material_properties(woven_matPar,sampleUC.wovenFaces,matName='Woven',surface=TOP_SURFACE)
V = V + sampleUC.getVolume(sampleUC.wovenFaces,woven_matPar['thickness'])
if sampleUC.knitFaces == None:
    print('No Knit Faces Detected, full Woven Material.')
else:
    sampleUC.material_properties(knit_matPar,sampleUC.knitFaces,matName='Knit',surface=MIDDLE_SURFACE)
    V = V + sampleUC.getVolume(sampleUC.knitFaces,knit_matPar['thickness'])

sampleUC.mesh([S4R,S3],QUAD_DOMINATED,0.35,perturbation=False) # ElemCode [S4R,S3],[CPE4H,CPE3H], ElemShape (TRI,QUAD,QUAD_DOMINATED)
sampleUC.assembly(sampleUC.name_part+'-1')
sampleUC.PBC_cyl2(sampleUC.xlatticeVec,sampleUC.ylatticeVec,'RefPoint-1','RefPoint-2','corner_nodes','edge_nodes',radius=57.15/2,offset=cyl_offset)#(sampleUC.xlatticeVec[0])*12/(2*pi))
sampleUC.analytical_cylinder(offset=-cyl_offset)
sampleUC.int_contact(tangential=True)
sampleUC.amplitude('amp-'+name_part,step_name='step-'+name_part,num_steps = 100, exp_fact = 100,num_csv_TP = 50,num_odb_TP = 50)    # Exp_fact 1 (linear ramp) up to 100 (Exponential)
sampleUC.PBC_BC_cyl2(DeformationF,'amp-'+name_part,'step-'+name_part,origin_fix=False)

sampleUC.create_job(name_part,nCores = 4,initial_conditions=False,file_subroutine=True,path_subroutine=parentdir+'\\subroutines\\MPC_PBC_CYL_2.f',MPC=True,Corr_Normal=True)
sim_time = sampleUC.run_job()

# Outputs
np.savetxt('woven_density.txt', [sampleUC.getVolume(sampleUC.wovenFaces,1)/area])
np.savetxt('sim_time.txt', [sim_time])
sampleUC.output_E_P_cyl('stress',0.5,area)
sampleUC.output_E_D_cyl()
sampleUC.output_E_V_cyl('velocity')
sampleUC.output_E_R_cyl()
sampleUC.output_dimple(file_name='dimple.csv',step=0)

# Clean up folder
sampleUC.clean_up(retainFiles=['.sta'])
