from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import math
import numpy as np
import displayGroupMdbToolset as dgm
import displayGroupOdbToolset as dgo
import assembly


#########################################################
#                   ABAQUS CLASS                        #
#########################################################

#   Basic abaqus class that builds a strip finite size
#   that is defined by a set of geoemetric paramter and
#   designed for material parameter optimizations.

#   INPUTS:
#   model - abaqus model variables
#   length - length of the strips minus the radius
#   width - width of the whole sample
#   n_strips - number of strips across the whole sample
#   t - width of the strips in the center length
#   radius - radius of the segments at the top and bottom
#   name_part - name of the part to be tested

class generate_strip:
    def __init__(self,model,length,width,n_strips,t,radius,name_part,dim_3D=False):
        # Unpack Variables
        self.model = model
        self.length = length
        self.width = width
        self.n_strips = n_strips
        self.t = t
        self.radius = radius
        self.name_part = name_part
        self.dim_3D = dim_3D
        self.generate_geometry(dim_3D)

    def generate_geometry(self,dim_3D = False):
        name_sketch = 'sketch'

        if dim_3D:
            self.model.ConstrainedSketch(name=name_sketch, sheetSize=200.0)
            self.model.sketches[name_sketch].rectangle(point1=(0, 0), point2=
                (self.length,self.width))
            # Generate part
            self.model.Part(dimensionality=THREE_D, name=self.name_part, type=
                DEFORMABLE_BODY)
            self.model.parts[self.name_part].BaseSolidExtrude(depth=self.t,sketch=
                self.model.sketches[name_sketch])
            self.model.parts[self.name_part].Set(faces=self.model.parts[self.name_part].faces.findAt(((self.length,self.width/2,self.t/2), )), name='edge_bot')
            self.model.parts[self.name_part].Set(faces=self.model.parts[self.name_part].faces.findAt(((0,self.width/2,self.t/2), )), name='edge_top')
            

        else:        
            self.model.ConstrainedSketch(name=name_sketch, sheetSize=200.0)

            # Draw Arcs
            self.model.sketches[name_sketch].ArcByCenterEnds(center=(0.0, 
                self.radius), direction=COUNTERCLOCKWISE, point1=(0.0, 0.0), point2=(self.radius, self.radius))
            self.model.sketches[name_sketch].ArcByCenterEnds(center=(self.radius*2+self.t, 
                self.radius+self.length), direction=CLOCKWISE, point1=(self.radius+self.t, self.radius+self.length), point2=(self.radius*2+self.t, self.radius*2+self.length))
            self.model.sketches[name_sketch].ArcByCenterEnds(center=(self.radius*2+self.t, 
                self.radius), direction=CLOCKWISE, point1=(self.radius*2+self.t, 0.0), point2=(self.radius+self.t, self.radius))
            self.model.sketches[name_sketch].ArcByCenterEnds(center=(0.0, 
                self.radius+self.length), direction=COUNTERCLOCKWISE, point1=(self.radius, self.radius+self.length), point2=(0.0, self.radius*2+self.length))

            # Sides of strip
            self.model.sketches[name_sketch].Line(point1=(self.radius, self.radius), point2=
                (self.radius, self.radius+self.length))
            self.model.sketches[name_sketch].Line(point1=(self.radius+self.t, self.radius), point2=
                (self.radius+self.t, self.radius+self.length))

            # Top and Bottom segmrenets & Related Sets
            bot = self.model.sketches[name_sketch].Line(point1=(0.0, 0.0), point2=(
                self.radius*2+self.t, 0.0))
            top = self.model.sketches[name_sketch].Line(point1=(0.0,  self.radius*2+self.length), point2=(
                self.radius*2+self.t,  self.radius*2+self.length))
            botPoint = bot.pointOn+(0,)
            topPoint = top.pointOn+(0,)

            # Generate part
            self.model.Part(dimensionality=THREE_D, name=self.name_part, type=
                DEFORMABLE_BODY)
            self.model.parts[self.name_part].BaseShell(sketch=
                self.model.sketches[name_sketch])
            self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((botPoint, )), name='edge_bot')
            self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((topPoint, )), name='edge_top')

            
            del self.model.sketches[name_sketch]
            return

    def material_properties(self,woven_matPar):
        ################################################
        # Material and Properties
        ################################################
        self.woven_matPar = woven_matPar
        self.Ktransverse_woven = (5./6.)*(3770)*woven_matPar["thickness"]

        # Woven Material and Section
        self.model.Material(name='Woven')


        if woven_matPar['matModel'] == 'NH':
            self.model.materials['Woven'].Hyperelastic(
                materialType=ISOTROPIC, 
                testData=OFF, 
                type=NEO_HOOKE, 
                volumetricResponse=VOLUMETRIC_DATA, 
                table=((woven_matPar["C10"], woven_matPar["D0"]), ))
        # Ogden
        if woven_matPar['matModel']=='OGD':
            self.model.materials['Woven'].Hyperelastic(materialType=ISOTROPIC, 
            table=((matPar['mu1'], matPar['alpha1'], matPar['D1']), ), testData=OFF, type=OGDEN, 
            volumetricResponse=VOLUMETRIC_DATA)
        if woven_matPar['matModel']=='AB':
            self.model.materials['Woven'].Hyperelastic(materialType=ISOTROPIC,
                table=((matPar['mu1'], matPar['alpha1'], matPar['D1']), ), testData=OFF, type=ARRUDA_BOYCE, 
                volumetricResponse=VOLUMETRIC_DATA)
        if woven_matPar['matModel'] == 'ORT':
            self.model.materials['Woven'].Elastic(moduli=INSTANTANEOUS, table=((D1111, D1122, D2222, D1133, D2233, D3333, 
                D1212, D1313, D2323), ), type=ORTHOTROPIC)
        if woven_matPar['matModel'] == 'ANI':
            self.model.materials['Woven'].Elastic(table=((D1111, D1122, D2222, D1133, D2233, D3333,
            D1112, D2212, D3312, D1212, D1113, D2213, D3313, D1213, D1313, D1123, D2223, D3323, 
            D1223, D1323, D2323), ), type=ANISOTROPIC)
        if woven_matPar['matModel'] == 'LAM':
            self.model.materials['Woven'].Elastic(table=((E1, E2, v12, 
                G12, G13, G23), ), type=LAMINA)
        self.model.materials['Woven'].Density(table=((self.woven_matPar['rho'], ), ))
        self.model.materials['Woven'].Damping(alpha=self.woven_matPar['alpha_damping'])
        if "DepVar" in woven_matPar:
            self.model.materials['Woven'].Depvar(n=woven_matPar['DepVar'])

        # Shell Section Applied
        if self.dim_3D:
            self.model.HomogeneousSolidSection(material='Woven', name=
                'Section_Woven', thickness=None)
            self.model.parts[self.name_part].Set(cells=
                self.model.parts[self.name_part].cells.findAt(((self.length/2, self.width/2, self.t/2), 
                )), name='Set-Cell')
            self.model.parts[self.name_part].SectionAssignment(offset=0.0, 
                offsetField='', offsetType=MIDDLE_SURFACE, region=
                self.model.parts[self.name_part].sets['Set-Cell'], sectionName=
                'Section_Woven', thicknessAssignment=FROM_SECTION)
        else:
            self.model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
                integrationRule=SIMPSON, material='Woven', name='Section_Woven', 
                nodalThicknessField='', numIntPts=5, poissonDefinition=VALUE, poisson=self.woven_matPar['PR'],
                preIntegrate=OFF, temperature=GRADIENT, thickness=self.woven_matPar["thickness"], thicknessField='', 
                thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
            # Apply properties to faces
            wovenFaces= self.model.parts[self.name_part].faces
            self.model.parts[self.name_part].SectionAssignment(offset=0.0, 
                offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                faces=wovenFaces), sectionName='Section_Woven', thicknessAssignment=
                FROM_SECTION)
            if "K11" in woven_matPar:
                self.model.sections['Section_Woven'].TransverseShearShell(k11=self.woven_matPar['K11'], 
                    k12=self.woven_matPar['K11']/1.2, k22=self.woven_matPar['K11'])
                    # Apply material orientation
            self.model.parts[self.name_part].MaterialOrientation(
                additionalRotationField='', additionalRotationType=ROTATION_ANGLE, angle=self.woven_matPar["material_orientation"]
                , axis=AXIS_3, fieldName='', localCsys=
                None, orientationType=
                SYSTEM, region=Region(
                faces=wovenFaces))
        
        return

    def mesh_3D(self):
        self.model.parts[self.name_part].seedPart(deviationFactor=0.1, 
            minSizeFactor=0.1, size=self.t/4)
        self.model.parts[self.name_part].setElementType(elemTypes=(ElemType(
            elemCode=C3D8RH, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, 
            hourglassControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
            ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(
            self.model.parts[self.name_part].cells.findAt(((self.length/2, self.width/2, self.t/2), 
            )), ))
        self.model.parts[self.name_part].generateMesh()
        self.model.parts[self.name_part].Set(nodes= self.model.parts[self.name_part].sets['edge_top'].nodes, name='top_nodes')
        return

    def mesh(self,elemCode,elemShape,mesh_size_tile):
        self.model.parts[self.name_part].setElementType(elemTypes=(ElemType(
            elemCode=elemCode, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
            hourglassControl=DEFAULT), ElemType(elemCode=elemCode, elemLibrary=STANDARD)), 
            regions=(self.model.parts[self.name_part].faces, ))

        self.model.parts[self.name_part].setMeshControls(elemShape=elemShape, regions=(self.model.parts[self.name_part].faces))

        # Set seed sizes
        mesh_tol = 0.0000000000001
        self.model.parts[self.name_part].seedEdgeBySize(constraint=FINER, 
        deviationFactor=mesh_tol, edges=
        self.model.parts[self.name_part].edges,
            minSizeFactor=mesh_tol, size=mesh_size_tile)
        self.model.parts[self.name_part].seedPart(deviationFactor=mesh_tol, 
            minSizeFactor=mesh_tol, size=mesh_size_tile)
        
        self.model.parts[self.name_part].generateMesh()

        self.model.parts[self.name_part].Set(nodes= self.model.parts[self.name_part].sets['edge_top'].nodes, name='top_nodes')

        return

    
    def generate_RP(self,NameRef1):
        self.NameRef1 = NameRef1
        self.model.Part(
            dimensionality=THREE_D, name=self.NameRef1, type=DEFORMABLE_BODY)
        self.model.parts[self.NameRef1].ReferencePoint(point=(self.radius+self.t/2, self.radius*2+self.length, 0.0))
        self.model.rootAssembly.Instance(dependent=ON, name=self.NameRef1, part=self.model.parts[self.NameRef1])

        # Create set of reference points
        self.model.rootAssembly.Set(name=NameRef1, referencePoints=(
            self.model.rootAssembly.instances[NameRef1].referencePoints[1],))

    def assembly(self,name_assem):
        # Generate Assembly
        self.name_assem = name_assem
        self.model.rootAssembly.Instance(dependent=ON, name=self.name_assem, part=self.model.parts[self.name_part])
        return

    def amplitude_exp(self,amp,disp,force):
        self.disp = disp
        self.force = force
        num_timesteps = len(disp)
        all_t = np.linspace(0,1,num_timesteps)
        all_tuples = num_timesteps*[0]
        for i in range(num_timesteps):
            all_tuples[i] = (all_t[i],disp[i]/disp[-1])
        self.model.TabularAmplitude(data=all_tuples, name=
            amp, timeSpan=STEP)

        time_TP = len(all_t)*[0]
        for i in range(len(all_t)):
            time_TP[i] = ((all_t[i]),)

        self.model.TimePoint(name='TimePoints-1', points=tuple(time_TP))

        #Increase minimum number of attempts from 5 to 30
        self.model.steps['step-'+self.name_part].control.setValues(allowPropagation=OFF, 
            resetDefaultValues=OFF, timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0, 
            12.0, 25.0, 6.0, 3.0, 50.0))

        self.model.HistoryOutputRequest(createStepName='step-'+self.name_part, name=
            'H-Output-2', rebar=EXCLUDE, region=   
            self.model.rootAssembly.sets[self.NameRef1], sectionPoints=DEFAULT, 
            variables=('U1','U2','U3','RF1','RF2','RF3'),timePoint='TimePoints-1')
        
        self.model.fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'V', 'A', 'RF', 'CF', 'CSTRESS', 
            'CDISP', 'STH', 'SE'))

    def amplitude(self,amp,num_steps,exp_fact,num_TP):
        # Exp_fact is numbers from 1 (linear ramp) to 100+ (exponential)
        self.num_steps = num_steps
        self.exp_fact = exp_fact
        amp_data = self.exp_amp(num_steps,exp_fact)

        self.model.TabularAmplitude(data=amp_data, name=
            amp, timeSpan=STEP)

        #Increase minimum number of attempts from 5 to 30
        self.model.steps['step-'+self.name_part].control.setValues(allowPropagation=OFF, 
            resetDefaultValues=OFF, timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0, 
            12.0, 25.0, 6.0, 3.0, 50.0))

        self.model.HistoryOutputRequest(createStepName='step-'+self.name_part, name=
            'H-Output-2', rebar=EXCLUDE, region=   
            self.model.rootAssembly.sets[self.NameRef1], sectionPoints=DEFAULT, 
            variables=('U1','U2','U3','RF1','RF2'),timePoint='TimePoints-1')
        
        self.model.fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'V', 'A', 'RF', 'CF', 'CSTRESS', 
            'CDISP', 'STH', 'SE'))

        self.num_TP = num_TP
        amp_TP = np.linspace(0,1,num_TP)

        time_TP = num_TP*[0]
        if exp_fact == 1:
            for i,t in enumerate(np.linspace(0,1,num_TP)):
                time_TP[i] = ((float(t)),)
        else:
            for i in range(num_TP):
                    if i != 0:
                        time_TP[i] = (((np.log10(amp_TP[i]*((exp_fact)-1)+1))/np.log10(exp_fact)),)
                    else:
                        time_TP[i] = ((0),)
        self.model.TimePoint(name='TimePoints-1', points=tuple(time_TP))

    def exp_amp(self,num_timesteps,b):
        def exp_fxn(x,b):
            if b ==1:
                return 1
            else:
                return 1./(b-1) * (b**x - 1)

        all_t = np.linspace(0,1,num_timesteps)
        all_tuples = num_timesteps*[0]
        for i in range(num_timesteps):
            all_tuples[i] = (all_t[i],exp_fxn(all_t[i],b))
        return all_tuples

    def boundary_conditions(self,DeformationF,amp,step):
        REF1 = self.model.rootAssembly.sets[self.NameRef1]
        TopSet = self.name_assem+'.top_nodes'
        for Dim1 in [1,2]:
            self.model.Equation(name='topRP-'+str(Dim1),terms=(((1, TopSet, Dim1),(-1, self.NameRef1, Dim1),)))
        # Top BC
        if self.dim_3D:
            self.model.DisplacementBC(amplitude=amp, createStepName=step, 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-TOP', region=REF1, 
                u1=-DeformationF, 
                u2=UNSET, 
                u3=UNSET, 
                ur1=UNSET,
                ur2=UNSET,
                ur3=UNSET)
            # Bot BC
            self.model.DisplacementBC(amplitude=amp, createStepName=step, 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-BOT', region=self.model.rootAssembly.sets[self.name_assem+".edge_bot"], 
                u1=0, 
                u2=UNSET, 
                u3=UNSET, 
                ur1=UNSET,
                ur2=UNSET,
                ur3=UNSET)

        else:
            self.model.DisplacementBC(amplitude=amp, createStepName=step, 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-TOP', region=REF1, 
                u1=0, 
                u2=DeformationF, 
                u3=UNSET, 
                ur1=UNSET,
                ur2=UNSET,
                ur3=UNSET)

            # Bot BC
            self.model.DisplacementBC(amplitude=amp, createStepName=step, 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-BOT', region=self.model.rootAssembly.sets[self.name_assem+".edge_bot"], 
                u1=0, 
                u2=0, 
                u3=UNSET, 
                ur1=UNSET,
                ur2=UNSET,
                ur3=UNSET)

    def boundary_conditions_gravity(self,DeformationF,amp,step):
        REF1 = self.model.rootAssembly.sets[self.NameRef1]
        TopSet = self.name_assem+'.top_nodes'
        for Dim1 in [1,2,3]:
            self.model.Equation(name='topRP-'+str(Dim1),terms=(((1, TopSet, Dim1),(-1, self.NameRef1, Dim1),)))
        self.model.TabularAmplitude(data=((0.0, 0.0), (1.0, 1.0), (2.0, 
            1.0)), name='Amp-Gravity', smooth=SOLVER_DEFAULT, timeSpan=STEP)
        self.model.Gravity(comp3=DeformationF, createStepName=step, 
            distributionType=UNIFORM, field='', name='Load-1',amplitude='Amp-Gravity')
        # Top BC
        self.model.DisplacementBC(amplitude=amp, createStepName=step, 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-TOP', region=REF1, 
            u1=UNSET, 
            u2=UNSET, 
            u3=UNSET, 
            ur1=UNSET,
            ur2=UNSET,
            ur3=UNSET)

        # Bot BC
        self.model.DisplacementBC(amplitude=amp, createStepName=step, 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-BOT', region=self.model.rootAssembly.sets[self.name_assem+".edge_bot"], 
            u1=0, 
            u2=0, 
            u3=0, 
            ur1=0,
            ur2=0,
            ur3=0)

    def create_job(self,jobName,nCores,file_subroutine=None):
        # Regenerate assembly
        self.model.rootAssembly.regenerate()
        self.jobName = jobName
        self.file_subroutine = file_subroutine
        self.nCores = nCores
        # Create temporary job
        jobName_tmp = 'job_tmp'
        mdb.Job(
            name=jobName_tmp,
            numCpus=self.nCores, 
            numDomains=self.nCores, 
            explicitPrecision=DOUBLE_PLUS_PACK,
            nodalOutputPrecision=SINGLE, 
            userSubroutine=self.file_subroutine,
            atTime=None, 
            contactPrint=OFF, 
            description='', 
            echoPrint=OFF, 
            getMemoryFromAnalysis=True, 
            historyPrint=OFF, 
            memory=90, 
            memoryUnits=PERCENTAGE, 
            model='Model-1', 
            modelPrint=OFF, 
            multiprocessingMode=THREADS, 
            numGPUs=0, 
            parallelizationMethodExplicit=DOMAIN, 
            activateLoadBalancing=False, 
            queue=None, 
            resultsFormat=ODB, 
            scratch='', 
            type=ANALYSIS, 
            waitHours=0, 
            waitMinutes=0)

        # Update input file to include material orientations
        mdb.jobs[jobName_tmp].writeInput()

        oldFile = open(jobName_tmp + '.inp', 'r')
        newFile = open(self.jobName + '.inp', 'w+')
        
        # linesToSkip = 0
        for line in oldFile:
            if line.startswith('*Orientation,'):
                newFile.write(line)
            elif line.startswith('*Transverse Shear'):
                newFile.write('*Transverse Shear Stiffness\n')
            else:
                newFile.write(line)

        oldFile.close()
        newFile.close()
        print(self.jobName)
        # Create final job
        mdb.JobFromInputFile(
            name=self.jobName, 
            inputFileName=self.jobName+'.inp',
            numCpus=self.nCores,
            numDomains=self.nCores,
            explicitPrecision=DOUBLE_PLUS_PACK,
            nodalOutputPrecision=SINGLE, 
            userSubroutine=self.file_subroutine,
            atTime=None, 
            getMemoryFromAnalysis=True, 
            memory=90, 
            memoryUnits=PERCENTAGE, 
            multiprocessingMode=THREADS, 
            numGPUs=0, 
            parallelizationMethodExplicit=DOMAIN, 
            activateLoadBalancing=False, 
            queue=None, 
            resultsFormat=ODB, 
            scratch='', 
            type=ANALYSIS, 
            waitHours=0, 
            waitMinutes=0)
        del mdb.jobs['job_tmp']

    def run_job(self,waitForCompletion=True):
        mdb.jobs[self.jobName].submit(
            consistencyChecking=OFF)
        print('Submitted: {}'.format(self.jobName))
        if waitForCompletion:
            print('Waiting: {}'.format(self.jobName))
            mdb.jobs[self.jobName].waitForCompletion()
        self.clean_up()

    def clean_up(self,retainFiles=['.odb','.sta','.msg']):
        extensions = ['.odb','.sta','.dat','.com','.ipm','.log','.prt','.sim',
            '.msg','.lck','.inp']
        for e in retainFiles:
            extensions.remove(e)

        for e in extensions:
            try:
                os.remove(self.jobName + e)
            except:
                pass
        try:
            os.remove('job_tmp.inp')
        except:
            pass

    def post_process_gravity(self,fileName,length):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        # frame of pre-loading step
        num=0
        dispX=session.XYDataFromHistory(name='dispX',odb=odb,outputVariableName='Spatial displacement: U1 PI: REFPOINT-0 Node 1 in NSET REFPOINT-0')
        dispY=session.XYDataFromHistory(name='dispY',odb=odb,outputVariableName='Spatial displacement: U2 PI: REFPOINT-0 Node 1 in NSET REFPOINT-0')
        dispZ=session.XYDataFromHistory(name='dispZ',odb=odb,outputVariableName='Spatial displacement: U3 PI: REFPOINT-0 Node 1 in NSET REFPOINT-0')

        RES=np.zeros((len(dispX)-num,3))
        for i in range(num,len(dispX)):
            RES[i-num][0]=dispX[i][1]
            RES[i-num][1]=length+dispY[i][1]
            RES[i-num][2]=-dispZ[i][1]

        np.savetxt(fileName+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return

    def post_process(self,fileName):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        # frame of pre-loading step
        num=0
        disp=session.XYDataFromHistory(name='disp',odb=odb,outputVariableName='Spatial displacement: U2 PI: REFPOINT-0 Node 1 in NSET REFPOINT-0')
        force = session.XYDataFromHistory(name='force',odb=odb,outputVariableName='Reaction force: RF2 PI: REFPOINT-0 Node 1 in NSET REFPOINT-0')

        RES=np.zeros((len(disp)-num,2))
        for i in range(num,len(disp)):
            RES[i-num][0]=disp[i][1]
            RES[i-num][1]=force[i][1]

        np.savetxt(fileName+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return disp,force

    def print_finite(self,name_print,name_file):

        self.model.Part(name=name_print, objectToCopy=
            self.model.parts[self.name_part])

        ################################################
        # Creating base fabric
        ################################################
        p = self.model.parts[name_print]
        session.viewports['Viewport: 1'].setValues(displayedObject=p)
        session.viewports['Viewport: 1'].partDisplay.setValues(renderStyle=WIREFRAME)
        session.printToFile(fileName=name_file, format=SVG, canvasObjects=(
            session.viewports['Viewport: 1'], ))
        print('SVG file '+'printSVG'+' has been saved.')

    def read_csv(self,fileName):
        # file = open("data.csv")

        arr = np.loadtxt(fileName,delimiter=",", dtype=str)
        disp = []
        force = []
        for point in arr:
            disp.append(float(point[0]))
            force.append(float(point[1]))
        return disp,force