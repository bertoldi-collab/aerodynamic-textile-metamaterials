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
import csv
import math
import numpy as np
import assembly 
import mesh
import time
import os

#########################################################
#                   ABAQUS CLASS                        #
#########################################################

#   Basic abaqus class that builds a metamaterial unit-cell
#   that is defined by a 10-parameter geoemetric system. By
#   varying geometric parameters and 

class unit_cell:
    def __init__(self,model,name_part):
        # Cordinate based options instead of getSequence
        session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
        # Unpack Variables
        self.model = model
        self.name_part = name_part

    def generate_parallelogram(self,angle_base,L_base,H_base,x_offset,dimension=THREE_D,cut=False,P1 = [0,0],P2=[0,0]):
        self.angle_base = angle_base
        self.H_base = H_base
        self.L_base = L_base
        self.x_offset = x_offset
        self.dimension = dimension
        self.cut = cut
        ################################################
        # Outline Geometry of the unit-cell in local coordinates
        ################################################
        
        #counter clock-wise points of the unit cell fabric
        self.xrec = [0, self.L_base, self.L_base+self.x_offset,self.x_offset]
        self.yrec = [0, 0, self.H_base, self.H_base]

        ################################################
        # Unit-cell Rotations into global coordinate
        ################################################

        # Apply transforms to variables (angle and origin of UC)
        cc, ss = np.cos(self.angle_base* np.pi / 180.), np.sin(self.angle_base* np.pi / 180.)
        self.rotMat = np.array(((cc, -ss, 0), (ss, cc, 0), (0, 0, 1)))

        # Rotate the base UC and the unit-cell points based on angles
        for i,(xreci,yreci) in enumerate(zip(self.xrec,self.yrec)):
            self.xrec[i],self.yrec[i],dz = np.dot(self.rotMat,np.array([xreci,yreci,0]))

        # Define Lattice Vectors
        self.xlatticeVec = [self.xrec[1],self.yrec[1]]
        self.ylatticeVec = [self.xrec[3],self.yrec[3]]

        ################################################
        # Creating unit cell fabric
        ################################################
        s = self.model.ConstrainedSketch(name='base_textile', sheetSize=50)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        bot = s.Line(point1=(self.xrec[0], self.yrec[0]), point2=(self.xrec[1], self.yrec[1]))
        right = s.Line(point1=(self.xrec[1], self.yrec[1]), point2=(self.xrec[2], self.yrec[2]))
        top = s.Line(point1=(self.xrec[2], self.yrec[2]), point2=(self.xrec[3], self.yrec[3]))
        left = s.Line(point1=(self.xrec[3], self.yrec[3]), point2=(self.xrec[0], self.yrec[0]))

        botPoint = bot.pointOn+(0,)
        rightPoint = right.pointOn+(0,)
        topPoint = top.pointOn+(0,)
        leftPoint = left.pointOn+(0,)

        p = self.model.Part(name=self.name_part, dimensionality=self.dimension, type=DEFORMABLE_BODY)
        self.model.parts[self.name_part].BaseShell(sketch=s)
        s.unsetPrimaryObject()

        # Create Sets for later use
        e = self.model.parts[self.name_part].edges
        self.model.parts[self.name_part].Set(edges=e,  name='edgesAll')
        self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((botPoint, )), name='edge_bot')
        self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((rightPoint, )), name='edge_right')
        self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((topPoint, )), name='edge_top')
        self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((leftPoint, )), name='edge_left')
        self.v = self.model.parts[self.name_part].vertices
        self.v_pointOn = self.v.pointsOn

        del self.model.sketches['base_textile']

        self.ff = self.model.parts[self.name_part].faces

        ################################################
        # Cut Unit-Cell
        ################################################
        if cut:
            self.cut_ellipse(P1,P2)

        print('Successfully created homogenous paralellogram '+self.name_part)
    
    def generate_RP(self,NameRef1,NameRef2):
        self.NameRef1 = NameRef1
        self.NameRef2 = NameRef2
        self.model.Part(
            dimensionality=self.dimension, name=self.NameRef1, type=DEFORMABLE_BODY)
        self.model.parts[self.NameRef1].ReferencePoint(point=(self.xrec[1]/2, self.yrec[2]*2, 0.0))
        self.model.Part(
            dimensionality=self.dimension, name=self.NameRef2, type=DEFORMABLE_BODY)
        self.model.parts[self.NameRef2].ReferencePoint(point=(self.xrec[1]/2, self.yrec[2]*2, 0.0))
        self.model.rootAssembly.Instance(dependent=ON, name=self.NameRef1, part=self.model.parts[self.NameRef1])
        self.model.rootAssembly.Instance(dependent=ON, name=self.NameRef2, part=self.model.parts[self.NameRef2])
        # Create set of reference points
        self.model.rootAssembly.Set(name=NameRef1, referencePoints=(
            self.model.rootAssembly.instances[NameRef1].referencePoints[1],))
        self.model.rootAssembly.Set(name=NameRef2, referencePoints=(
            self.model.rootAssembly.instances[NameRef2].referencePoints[1],))

    def generate_RP_cyl(self,NameRef2):
        self.NameRef2 = NameRef2
        self.model.Part(
            dimensionality=self.dimension, name=self.NameRef2, type=DEFORMABLE_BODY)
        self.model.parts[self.NameRef2].ReferencePoint(point=(self.xrec[1]/2, self.yrec[2]*2, 0.0))
        self.model.rootAssembly.Instance(dependent=ON, name=self.NameRef2, part=self.model.parts[self.NameRef2])
        # Create set of reference points
        self.model.rootAssembly.Set(name=NameRef2, referencePoints=(
            self.model.rootAssembly.instances[NameRef2].referencePoints[1],))

    def RE_dimple_tracker(self,width,height,alpha,t):
        ################################################
        # Overlay Hexagonal Sketch and Partition
        ################################################
        s2 = self.model.ConstrainedSketch(gridSpacing=5, name='dimple_partition', 
            sheetSize=4*width)
        # Calculate x-offset based on geometry
        x_offset = height/(2*np.tan(alpha))


        dimple_line_1 = s2.Line(point1=(0,0), point2=(width*2,0))
        dimple_line_2 = s2.Line(point1=(0,0+height/2), point2=(width+width-x_offset,0+height/2))
        dimple_1_point = dimple_line_1.pointOn+(0,)
        dimple_2_point = dimple_line_2.pointOn+(0,)

        s2.setPrimaryObject(option=STANDALONE)
        self.model.parts[self.name_part].PartitionFaceBySketch(faces=
            self.model.parts[self.name_part].faces, sketch=s2)


        # Create Sets for later use
        e = self.model.parts[self.name_part].edges
        self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((dimple_1_point, )), name='dimple_1')
        self.model.parts[self.name_part].Set(edges=self.model.parts[self.name_part].edges.findAt((dimple_2_point, )), name='dimple_2')



    def RE_partition(self,width,height,alpha,t):
        self.width = width
        self.height = height
        self.alpha = alpha
        self.t = t

        # Calculate x-offset based on geometry
        x_offset = height/(2*np.tan(alpha))
        # Create connecting path to generate both Hexagon 1 and Hexagon 2
        self.hexPoints1x = [0,  width,  width-x_offset, width,0,x_offset,0]
        self.hexPoints1y = [0,0,height/2,height,height,height/2,0]

        self.hexPoints2x = [0+width-x_offset,  width+width-x_offset,  width-x_offset+width-x_offset, width+width-x_offset,0+width-x_offset,x_offset+width-x_offset,0+width-x_offset]
        self.hexPoints2y = [0+height/2,0+height/2,height/2+height/2,height+height/2,height+height/2,height/2+height/2,0+height/2]

        self.pattern_hypot = hypot(0,height)

        ################################################
        # Overlay Hexagonal Sketch and Partition
        ################################################
        s1 = self.model.ConstrainedSketch(gridSpacing=5, name='hexagon_tess', 
            sheetSize=4*width)

        s1.setPrimaryObject(option=STANDALONE)

        # Grid generation
        n_p = 2  # range from -n to n
        xp = n_p*2+1  # number of points along x-axis
        yp = n_p*2+1  # number of points along y-axis
        # Generate a grid of integers from -n to n
        x_values = np.linspace(-n_p, n_p, xp)
        y_values = np.linspace(-n_p, n_p, yp)
        x_grid, y_grid = np.meshgrid(x_values, y_values)
        
        # Flatten the grids
        xGrid = x_grid.flatten()
        yGrid = y_grid.flatten()
                
        #counter clock-wise points of the unit cell fabric
        xrec = [0, width*2-x_offset*2, width*2-x_offset*2,0]
        yrec = [0, 0, height, height]

        ################################################
        # Unit-cell Rotations into global coordinate
        ################################################

        # Define Lattice VectorsF
        xlatticeVec = [xrec[1],yrec[1]]
        ylatticeVec = [xrec[3],yrec[3]]

        count = 0
        for x_off,y_off in zip(xGrid,yGrid):   
            count = count+ 1
            L1x_off = xlatticeVec[0]*x_off
            L1y_off = xlatticeVec[1]*x_off
            L2x_off = ylatticeVec[0]*y_off
            L2y_off = ylatticeVec[1]*y_off
            # Create empty object lists
            hex1 = []#np.zeros(6,dtype=dict)
            hex2 = []#np.zeros(6,dtype=dict)
            #Generate Hexagon lines for 3x3 Unit-cell grid
            for j in range(0,6,1):
                if np.sqrt((self.hexPoints1x[j]+L1x_off+L2x_off-(self.hexPoints1x[j+1]+L1x_off+L2x_off))**2+(self.hexPoints1y[j]+L1y_off+L2y_off-(self.hexPoints1y[j+1]+L1y_off+L2y_off))**2)>0.01:
                    new_line = s1.Line(point1=(self.hexPoints1x[j]+L1x_off+L2x_off, self.hexPoints1y[j]+L1y_off+L2y_off), point2=(self.hexPoints1x[j+1]+L1x_off+L2x_off, self.hexPoints1y[j+1]+L1y_off+L2y_off))
                    hex1.append(new_line)
                if np.sqrt((self.hexPoints2x[j]+L1x_off+L2x_off-(self.hexPoints2x[j+1]+L1x_off+L2x_off))**2+(self.hexPoints2y[j]+L1y_off+L2y_off-(self.hexPoints2y[j+1]+L1y_off+L2y_off))**2)>0.01:
                    new_line = s1.Line(point1=(self.hexPoints2x[j]+L1x_off+L2x_off, self.hexPoints2y[j]+L1y_off+L2y_off), point2=(self.hexPoints2x[j+1]+L1x_off+L2x_off, self.hexPoints2y[j+1]+L1y_off+L2y_off))
                    hex2.append(new_line)

            n_p1 = len(s1.geometry)+2 # length of geometry preoffset

            
            s1.offset(distance=self.t, objectList=list(hex1), side=RIGHT)

            try:
                s1.offset(distance=self.t, objectList=list(hex2), side=RIGHT)
            except:
                pass
            n_p2 = len(s1.geometry)+1 # length of geometry preoffset
            # Delete Guidelines for both Hexagons
            s1.delete(objectList=list(hex1))
            s1.delete(objectList=list(hex2))
            # if count == 2:
            #     print(xbox)
            # print(count)
        try:
            self.model.parts[self.name_part].PartitionFaceBySketch(faces=
                self.model.parts[self.name_part].faces, sketch=s1)
        except:
            pass

        del self.model.sketches['hexagon_tess']
        self.TM_remove_faces(self.model.parts[self.name_part])

        self.TM_find_faces_RE(xrec,yrec)

        print("Textile Metamaterial pattern successfully partitioned on unit-cell.")

        return self.xrec,self.yrec


    def TM_remove_faces(self,part):
        tol = 0.0001
        for f in part.faces:
            if f.getSize(printResults=False) < tol:
                edges = f.getEdges()
                edgeTuple = ()
                for e in edges:
                    edgeTuple = edgeTuple+(part.edges[e],)
                part.RemoveRedundantEntities(edgeList=edgeTuple)

    def midpoint(self,p1, p2):
        return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2]

    def isBetween(self,a, b, c):
        if max(a[0],b[0]) >= c[0] and min(a[0],b[0]) <= c[0] and max(a[1],b[1]) >= c[1] and min(a[1],b[1]) <= c[1]:
            return True
        else:
            return False


    def TM_find_faces_RE(self,xrec,yrec):
        # Determine which face belongs to base, or to woven
        ff =self.model.parts[self.name_part].faces
        # Woven Faces
        self.wovenFaces = ff.findAt(((0,0,0),),((xrec[1],yrec[1],0),), ((xrec[2],yrec[2],0),), ((xrec[3],yrec[3],0),))
        # self.wovenFaces = ff.findAt((((0,0,0),(xrec[1],yrec[1],0),(xrec[2],yrec[2],0),(xrec[3],yrec[3],0),),))
        self.knitFaces = FaceArray([f for f in ff if f not in self.wovenFaces])

    def material_properties(self,matPar,faces,matName='Material',surface=MIDDLE_SURFACE):
        # Woven Material and Section
        self.model.Material(name=matName)
        Ktransverse = (5./6.)*(3770)*matPar["thickness"]

        # Damping and Density Parameters
        try:
            self.model.materials[matName].Density(table=((matPar['rho'], ), ))
            self.model.materials[matName].Damping(alpha=matPar['alpha_damping'])
        except:
            print("Warning: Density and Damping parameters not implemented.")

        ################################################
        # Material Model: Parameter Assignment
        ################################################

        # Holzapfel Gausser Ogden
        if matPar['matModel']=='HGO':
            # Repackage material properties
            nFibers = int(len([k for k in matPar.keys() if k.startswith('k')])/2)
            k1_vals = [matPar['k1{}'.format(i+1)] for i in range(nFibers)]
            k2_vals = [matPar['k2{}'.format(i+1)] for i in range(nFibers)]
            matPropVec = [matPar['C10'], matPar['D0']] +       \
            [val for pair in zip(k1_vals, k2_vals) for val in pair]
            self.model.materials[matName].UserMaterial(
                mechanicalConstants=tuple(matPropVec))

        if matPar['matModel']=='HGO_ANISO':
            # Repackage material properties
            nFibers = int(len([k for k in matPar.keys() if k.startswith('k')])/2)
            k1_vals = [matPar['k1{}'.format(i+1)] for i in range(nFibers)]
            k2_vals = [matPar['k2{}'.format(i+1)] for i in range(nFibers)]
            matPropVec = [matPar['C10'], matPar['D0']] +       \
            [val for pair in zip(k1_vals, k2_vals) for val in pair]
            self.model.materials[matName].Hyperelastic(anisotropicType=
            USER_DEFINED, materialType=ANISOTROPIC, moduliTimeScale=INSTANTANEOUS,localDirections=2, formulation=
            INVARIANT, properties=len(matPropVec), table=(tuple(matPropVec),))

        # Neohookean
        if matPar['matModel']=='NH':
            matPropVec = [matPar['C10'], matPar['D0']]
            self.model.materials[matName].Hyperelastic(
                materialType=ISOTROPIC, 
                testData=OFF, 
                type=NEO_HOOKE, 
                volumetricResponse=VOLUMETRIC_DATA, 
                table=((matPar['C10'], matPar['D0']), ))

        if matPar['matModel']=='LEO':
            self.model.materials['Woven'].Elastic(moduli=INSTANTANEOUS, table=((matPar['D1111'], matPar['D1122'], matPar['D2222'], matPar['D1133'], 
            matPar['D2233'], matPar['D3333'],matPar['D1212'], matPar['D1313'], matPar['D2323']), ), type=ORTHOTROPIC)

        if matPar['matModel'] == 'LAM':
            self.model.materials['Woven'].Elastic(table=((matPar['E1'], matPar['E2'], matPar['v12'], 
                matPar['G12'], matPar['G13'], matPar['G23']), ), type=LAMINA)
        
        # Ogden
        if matPar['matModel']=='OGD':
            self.model.materials[matName].Hyperelastic(materialType=ISOTROPIC, 
            table=((matPar['mu1'], matPar['alpha1'], matPar['D1']), ), testData=OFF, type=OGDEN, 
            volumetricResponse=VOLUMETRIC_DATA)

        # Arruda-Boyce
        if matPar['matModel']=='AB':
            self.model.materials[matName].Hyperelastic(materialType=ISOTROPIC,
                table=((matPar['mu1'], matPar['alpha1'], matPar['D1']), ), testData=OFF, type=ARRUDA_BOYCE, 
                volumetricResponse=VOLUMETRIC_DATA)

        if "DepVar" in matPar:
            self.model.materials[matName].Depvar(n=woven_matPar['DepVar'])

        ################################################
        # Section Generation
        ################################################

        if self.dimension == THREE_D:
            self.model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
                integrationRule=SIMPSON, material=matName, name='Section'+matName, 
                nodalThicknessField='', numIntPts=5, poissonDefinition=VALUE, poisson=matPar['PR'],
                preIntegrate=OFF, temperature=GRADIENT, thickness=matPar['thickness'], thicknessField='', 
                thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
        else:
            self.model.HomogeneousSolidSection(material=matName, name='Section'+matName
                , thickness=matPar['thickness'])

        if matPar['matModel']=='HGO' or matPar['matModel']=='HGO_ANISO' :
            self.model.sections['Section'+matName].TransverseShearShell(k11=Ktransverse,
                k12=Ktransverse/1.2, k22=Ktransverse)

        ################################################
        # Section Assignment
        ################################################
        self.model.parts[self.name_part].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=surface, region=Region(
            faces=faces), sectionName='Section'+matName, thicknessAssignment=
            FROM_SECTION)

        if "material_orientation" in matPar:
            self.model.parts[self.name_part].MaterialOrientation(
                additionalRotationField='', additionalRotationType=ROTATION_ANGLE, angle=matPar["material_orientation"]
                , axis=AXIS_3, fieldName='', localCsys=
                None, orientationType=
                SYSTEM, region=Region(
                faces=faces))
        
        print('Successfully created and assigned material properties for '+matName+' material.')

    def getVolume(self,faces,thickness):
        # Get the volume of a series of faces
        faces_singular = []
        [faces_singular.append(x) for x in faces if x not in faces_singular]
        V = 0
        faceList = []
        for f in faces_singular:
            if f.pointOn[0] != faceList:
                V = V + f.getSize(printResults=False)*thickness
                faceList.append(f.pointOn[0])
        return V

    def mesh(self,elemCode,elemShape,mesh_size_tile,perturbation=TRUE):
        self.mesh_size_tile = mesh_size_tile
        self.model.parts[self.name_part].setElementType(elemTypes=(ElemType(
            elemCode=elemCode[0], elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
            hourglassControl=DEFAULT), ElemType(elemCode=elemCode[1], elemLibrary=STANDARD)), 
            regions=(self.model.parts[self.name_part].faces, ))
        if elemShape == QUAD_DOMINATED:
            self.model.parts[self.name_part].setMeshControls(algorithm=ADVANCING_FRONT,
            allowMapped=False,elemShape=elemShape, regions=(self.model.parts[self.name_part].faces))
        else:
            self.model.parts[self.name_part].setMeshControls(elemShape=elemShape, regions=(self.model.parts[self.name_part].faces))
        # Set seed sizes
        mesh_tol = 1e-12
        mesh_dev = 1e-5
        if self.cut:
            for e in self.model.parts[self.name_part].edges:
                if e not in self.model.parts[self.name_part].sets['ellipse_face'].edges:
                    self.model.parts[self.name_part].seedEdgeBySize(constraint=FINER, 
                    deviationFactor=mesh_dev, edges=
                    (e,), minSizeFactor=mesh_tol, size=mesh_size_tile)
            for e in self.model.parts[self.name_part].sets['ellipse_face'].edges:
                self.model.parts[self.name_part].seedEdgeBySize(constraint=FIXED, edges=
                    (e,), size=0.6)
        else:
            self.model.parts[self.name_part].seedEdgeBySize(constraint=FINER, 
                    deviationFactor=mesh_tol, edges=
                    self.model.parts[self.name_part].edges, minSizeFactor=mesh_tol, size=mesh_size_tile)
        
        self.model.parts[self.name_part].seedPart(deviationFactor=mesh_dev, 
            minSizeFactor=mesh_tol, size=mesh_size_tile)
        self.model.parts[self.name_part].generateMesh(seedConstraintOverride=ON)

        # Create Node Sets
        self.model.parts[self.name_part].Set(nodes= self.model.parts[self.name_part].sets['edgesAll'].nodes, name='edge_nodes')
      	# Corner Set
        bb_size = mesh_size_tile/5.0
        n = self.model.parts[self.name_part].nodes
        self.model.parts[self.name_part].Set(name='cn_TL', nodes=n.getByBoundingBox(self.xrec[3]-bb_size, self.yrec[3]-bb_size, 0, self.xrec[3] + bb_size, self.yrec[3] + bb_size, 0))
        self.model.parts[self.name_part].Set(name='cn_TR', nodes=n.getByBoundingBox(self.xrec[2]-bb_size, self.yrec[2]-bb_size, 0, self.xrec[2] + bb_size, self.yrec[2] + bb_size, 0))
        self.model.parts[self.name_part].Set(name='cn_BR', nodes=n.getByBoundingBox(self.xrec[1]-bb_size, self.yrec[1]-bb_size, 0, self.xrec[1] + bb_size, self.yrec[1] + bb_size, 0))
        self.model.parts[self.name_part].Set(name='origin', nodes=n.getByBoundingBox(self.xrec[0]-bb_size, self.yrec[0]-bb_size, 0, self.xrec[0] + bb_size, self.yrec[0] + bb_size, 0))
        self.model.parts[self.name_part].Set(nodes= list((self.model.parts[self.name_part].sets['origin'].nodes,self.model.parts[self.name_part].sets['cn_TL'].nodes,self.model.parts[self.name_part].sets['cn_TR'].nodes,self.model.parts[self.name_part].sets['cn_BR'].nodes)), name='corner_nodes')

        if perturbation:
            for n in self.model.parts[self.name_part].nodes:
                if n not in self.model.parts[self.name_part].sets['edge_nodes'].nodes:
                    self.model.parts[self.name_part].editNode(
                        nodes=n, 
                        offset3=1e-6*np.random.randn())

    def assembly(self,name_assem):
        # Generate Assembly
        self.name_assem = name_assem
        self.model.rootAssembly.Instance(dependent=ON, name=self.name_assem, part=self.model.parts[self.name_part])

    
    def wrapMesh_custom(self,part,radius):
        offset = 0
        for n in part.nodes:
            [x,y,z] = n.coordinates
            x = x+1
            theta = x/radius+np.pi/4
            part.editNode(nodes=n,coordinate1=radius*np.cos(theta))
            part.editNode(nodes=n,coordinate2=radius*np.sin(theta))
            part.editNode(nodes=n,coordinate3=y)



    def PBC_cyl2(self,latticeVec1,latticeVec2,NameRef1,NameRef,NameCornerNodes,NameEdgeNodes,radius,offset = 0):
        # Wrap Part
        self.offset = offset
        self.radius = radius+offset
        radius = self.radius
        self.name_cyl = self.name_part+'-cyl'
        self.model.parts[self.name_part].PartFromMesh(copySets=True, name=self.name_cyl)
        self.model.parts[self.name_cyl].wrapMesh(radius=self.radius)
        # self.model.parts[self.name_part].PartFromMesh(copySets=True, name=self.name_cyl+'-2')
        # self.wrapMesh_custom(self.model.parts[self.name_cyl],self.radius)
        try:
            self.model.parts[self.name_cyl].sectionAssignments[0].setValues(
                offset=0.35, offsetField='', offsetType=SINGLE_VALUE)
            self.model.parts[self.name_cyl].sectionAssignments[1].setValues(
                offset=0, offsetField='', offsetType=SINGLE_VALUE)
        except:
            pass
        # Create Cylindrical Boundary Conditions
        del self.model.rootAssembly.features[self.name_part+'-1']
        self.model.rootAssembly.DatumCsysByDefault(CYLINDRICAL)
        self.model.rootAssembly.Instance(dependent=ON, name=self.name_cyl+'-1', 
            part=self.model.parts[self.name_cyl])

        self.model.rootAssembly.Set(name='nodes-all', nodes=
            self.model.rootAssembly.instances[self.name_cyl+'-1'].nodes)
        # NameEdgeNodes, name of the set for all edge nodes (including corner)
        ############################
        # Pre-Processing
        ############################
        # Nodes
        corner_nodes =  self.model.rootAssembly.sets[self.name_cyl+'-1.'+NameCornerNodes].nodes
        edge_nodes =    self.model.rootAssembly.sets[self.name_cyl+'-1.'+NameEdgeNodes].nodes
        # Create a set for the origin nodes
        for node in corner_nodes:
            if node in self.model.rootAssembly.sets[self.name_cyl+'-1.origin'].nodes:
                self.model.rootAssembly.Set(name='MPC-origin', nodes=self.NodeSequenceCyl(node))
                origin_node = self.model.rootAssembly.sets['MPC-origin'].nodes
        if self.dimension == THREE_D:
            dim_matrix = [4, 5, 6]
        elif self.dimension == TWO_D_PLANAR:
            dim_matrix = []
        else:
            raise BaseException('Error: Dimension not defined for PBC.')
        self.latticeVec1 = [0,latticeVec1[0],latticeVec1[1]]
        self.latticeVec2 = [0,latticeVec2[0],latticeVec2[1]]
        
        LatticeVec = (self.latticeVec1,self.latticeVec2)

        ############################
        # Corner Boundary Conditions
        ############################
        # Find the cylindrical coordinate lattice vectors
        # Note: lattice-vectors must be provided in standard PBC
        datum_id = self.model.rootAssembly.features['Datum csys-1'].id
        self.datum_id = datum_id
        usedNodes = []
        repConst = 1
        usedNodes.append(origin_node[0])
        # For every other corner node, tie it back to the origin.
        coordR_origin = origin_node[0].coordinates[0]
        coordPhi_origin = np.arctan(origin_node[0].coordinates[1]/origin_node[0].coordinates[0])    # X-coordinate: Always length dx with respect to origin
        coordZ_origin = origin_node[0].coordinates[2]    # Y-coordinate: Always length dy wtih respect to origin
        for node in corner_nodes:
            if node not in origin_node:
                coordR = node.coordinates[0]
                coordPhi = np.arctan(node.coordinates[1]/node.coordinates[0])    # X-coordinate: Always length dx with respect to origin
                coordZ = node.coordinates[2]    # Y-coordinate: Always length dy wtih respect to origin

                dR = coordR-coordR_origin
                dPhi = (coordPhi-coordPhi_origin)*radius
                dZ = coordZ-coordZ_origin
                # Create a set with current node to be used in equations
                self.origin_node = origin_node
                # self.model.rootAssembly.Set(name='MPC-origin-' + str(
                #     repConst), nodes=origin_node)
                self.model.rootAssembly.Set(name='MPC-corner-node-' + str(
                    repConst), nodes=self.NodeSequenceCyl(node))
                self.model.MultipointConstraint(controlPoint=
                    self.model.rootAssembly.sets[NameRef1], csys=None, mpcType=
                    USER_MPC, name='MPC-corner-' + str(repConst), surface=
                    self.model.rootAssembly.sets['MPC-corner-node-' + str(
                    repConst)], userMode=NODE_MODE_MPC
                    , userType=0)

                repConst = repConst + 1     # Increase integer for naming equation constraint
                usedNodes.append(node)

        ############################
        # Edge Boundary Conditions
        ############################
        repConst = 1
        for node1 in edge_nodes:
            if node1 not in corner_nodes:
                if node1 in usedNodes:
                    pass
                else:
                    stop = False
                    # Find Node1 Coordinates
                    coordR1 = node1.coordinates[0]
                    coordPhi1 =  np.arctan(node1.coordinates[1]/node1.coordinates[0])
                    coordZ1 = node1.coordinates[2]
                    for node2 in edge_nodes:
                        if node2 not in corner_nodes:
                            # Find Node2 Coordinates
                            coordR2 = node1.coordinates[0]
                            coordPhi2 =  np.arctan(node2.coordinates[1]/node2.coordinates[0])
                            coordZ2 = node2.coordinates[2]
                            # Find distance between nodes
                            dR = coordR2 - coordR1
                            dPhi = (coordPhi2 - coordPhi1) # X-Distance between nodes
                            dZ = coordZ2 - coordZ1  # Y-Distance between nodes
                            for vec in LatticeVec:
                                tol = 5e-5
                                if abs(dPhi*radius - vec[1]) < tol and abs(dZ - vec[2]) < tol:
                                    if stop:
                                        break
                                    else:
                                        # Correct combination found begin creating sets for use in equations constraints
                                        self.model.rootAssembly.Set(name='MPC-node-1-' + str(
                                            repConst), nodes=self.NodeSequenceCyl(node1))
                                        self.model.rootAssembly.Set(name='MPC-node-2-' + str(
                                            repConst), nodes=self.NodeSequenceCyl(node2))
                                        self.model.MultipointConstraint(controlPoint=
                                            self.model.rootAssembly.sets[NameRef1], csys=None, mpcType=
                                            USER_MPC, name='MPC-nodes-' + str(repConst), surface=
                                            self.model.rootAssembly.sets['MPC-node-2-' + str(
                                            repConst)], userMode=NODE_MODE_MPC
                                            , userType=0)
                                        # Remove used node from available list
                                        usedNodes.append(node1)
                                        usedNodes.append(node2)
                                        repConst = repConst + 1  # Increase integer for naming equation constraint
                                        # Don't look further, go to following node.
                                        stop = True
        
        if self.PBC_check(usedNodes,edge_nodes):
            raise BaseException("Error: PBC not implemeneted correctly, stopping python script.")

        print('Periodic Boundary Conditions successfully implemented')
        return 
    

    def PBC_planar(self,latticeVec1,latticeVec2,NameRef1,NameRef,NameCornerNodes,NameEdgeNodes):
        #Variables
        # latticeVec1 & latticeVec2 two dimensional list of x and y coordinates
        # exp: latticeVec1 = [x0,y0]
        # NameRef1 and NameRef, names of the reference points to which constraints will be tied
        # NameCornerNodes, name of the set for corner nodes
        # NameEdgeNodes, name of the set for all edge nodes (including corner)
        ############################
        # Pre-Processing
        ############################
        # Note: lattice-vectors must be provided in standard PBC
        LatticeVec = (latticeVec1,latticeVec2)
        # Nodes
        corner_nodes =  self.model.rootAssembly.sets[self.name_part+'-1.'+NameCornerNodes].nodes
        edge_nodes =    self.model.rootAssembly.sets[self.name_part+'-1.'+NameEdgeNodes].nodes
        # Create a set for the origin nodes
        for node in corner_nodes:
            if node.coordinates == (0,0,0):
                self.model.rootAssembly.Set(name='corner-origin', nodes=self.NodeSequence(node))
                origin_node = self.model.rootAssembly.sets['corner-origin'].nodes
        if self.dimension == THREE_D:
            dim_matrix = [3, 4, 5, 6]
        elif self.dimension == TWO_D_PLANAR:
            dim_matrix = []
        else:
            raise BaseException('Error: Dimension not defined for PBC.')

        ############################
        # Corner Boundary Conditions
        ############################
        usedNodes = []
        repConst = 1
        usedNodes.append(origin_node[0])
        # For every other corner node, tie it back to the origin.
        for node in corner_nodes:
            if node not in origin_node:
                coordX = node.coordinates[0]    # X-coordinate: Always length dx with respect to origin
                coordY = node.coordinates[1]    # Y-coordinate: Always length dy wtih respect to origin
                # Create a set with current node to be used in equations
                self.model.rootAssembly.Set(name='corner-node-' + str(
                    repConst), nodes=self.NodeSequence(node))
                for Dim1 in [1,2]:          # Create planar equations
                    self.model.Equation(name='corner-plrDOF-' + str(repConst) + '-' + str(Dim1),
                                                terms=((-1.0, 'corner-node-' + str(repConst), Dim1), (1.0, 'corner-origin', Dim1),
                                                        (coordX, NameRef1, Dim1), (coordY, NameRef, Dim1)))
                for Dim1 in dim_matrix:     # Create equality equations
                    self.model.Equation(name='corner-rotDOF-' + str(Dim1) + '-' + str(repConst),
                                                terms=((-1.0, 'corner-node-' + str(repConst), Dim1), (1.0, 'corner-origin', Dim1)))    
                repConst = repConst + 1     # Increase integer for naming equation constraint
                usedNodes.append(node)

        ############################
        # Edge Boundary Conditions
        ############################
        repConst = 1
        for node1 in edge_nodes:
            if node1 not in corner_nodes:
                if node1 in usedNodes:
                    pass
                else:
                    stop = False
                    # Find Node1 Coordinates
                    coordX1 = node1.coordinates[0]
                    coordY1 = node1.coordinates[1]
                    for node2 in edge_nodes:
                        if node2 not in corner_nodes:
                            # Find Node2 Coordinates
                            coordX2 = node2.coordinates[0]
                            coordY2 = node2.coordinates[1]
                            # Find distance between nodes
                            dx = coordX2 - coordX1  # X-Distance between nodes
                            dy = coordY2 - coordY1  # Y-Distance between nodes
                            for vec in LatticeVec:
                                if round(10000*(vec[0]-dx)) == 0 and round(10000*(vec[1]-dy)) == 0:
                                    if stop:
                                        break
                                    else:
                                        # Correct combination found begin creating sets for use in equations constraints
                                        self.model.rootAssembly.Set(name='Node-1-' + str(
                                            repConst), nodes=self.NodeSequence(node1))
                                        self.model.rootAssembly.Set(name='Node-2-' + str(
                                            repConst), nodes=self.NodeSequence(node2))
                                        for Dim1 in [1, 2]:
                                            self.model.Equation(name='node-plrDOF-' + str(repConst) + '-' + str(Dim1),
                                                                        terms=((1.0, 'Node-1-' + str(repConst), Dim1), (-1.0, 'Node-2-' + str(repConst), Dim1),
                                                                                (dx, NameRef1, Dim1), (dy, NameRef, Dim1)))
                                        for Dim1 in dim_matrix:
                                            self.model.Equation(name='node-rotDOF-' + str(repConst) + '-' + str(Dim1),
                                                                        terms=((-1.0, 'Node-1-' + str(repConst), Dim1), (1.0, 'Node-2-' + str(repConst), Dim1)))   
                                                                            # Remove used node from available list
                                        usedNodes.append(node1)
                                        usedNodes.append(node2)
                                        repConst = repConst + 1  # Increase integer for naming equation constraint
                                        # Don't look further, go to following node.
                                        stop = True
        
        if self.PBC_check(usedNodes,edge_nodes):
            raise BaseException("Error: PBC not implemeneted correctly, stopping python script.")

        print('Periodic Boundary Conditions successfully implemented')
        return 

    def PBC_check(self,used_nodes,edge_nodes):
        node_check = np.zeros(len(edge_nodes))

        for i,node in enumerate(edge_nodes):
            if node in used_nodes:
                node_check[i] = 1
            else:
                print('Node PBC: ')
                print(node)

        if 0 in node_check:
            return True
        else:
            return False

    def NodeSequence(self,node):
        # Example of Abaqus being dumb:
        # When creating a set, abaqus wants node objects, however
        # if given a node object as instancePart[2] it will fail
        # creation. Otherwise instancePart[2:3] is the correct version
        instancePart = self.model.rootAssembly.instances[self.name_part+'-1'].nodes
        return instancePart[node.label-1:node.label]
    
    def NodeSequenceCyl(self,node):
        # Example of Abaqus being dumb:
        # When creating a set, abaqus wants node objects, however
        # if given a node object as instancePart[2] it will fail
        # creation. Otherwise instancePart[2:3] is the correct version
        instancePart = self.model.rootAssembly.instances[self.name_cyl+'-1'].nodes
        return instancePart[node.label-1:node.label]
        
    def amplitude_exp(self,step_name,amp,disp,force):
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
            'H-Output-2-'+step_name, rebar=EXCLUDE, region=   
            self.model.rootAssembly.sets[self.NameRef1], sectionPoints=DEFAULT, 
            variables=('U1','U2','U3','RF1','RF2','RF3'),timePoint='TimePoints-1')

        self.model.HistoryOutputRequest(createStepName='step-'+self.name_part, name=
            'H-Output-3-'+step_name, rebar=EXCLUDE, region=   
            self.model.rootAssembly.sets[self.NameRef2], sectionPoints=DEFAULT, 
            variables=('U1','U2','U3','RF1','RF2','RF3'),timePoint='TimePoints-1')
        
        self.model.fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'V', 'A', 'RF', 'CF', 'CSTRESS', 
            'CDISP', 'STH', 'SE','COORD','IVOL'),timePoint='TimePoints-1')

    def amplitude(self,amp,num_steps,exp_fact,num_csv_TP,num_odb_TP,step_name='step-',start_point=0,end_point=1,amp_tuple=0,maxU=25,offset=0):
        # Exp_fact is numbers from 1 (linear ramp) to 100+ (exponential)
        #step_name = step_name+self.name_part
        self.num_steps = num_steps
        self.exp_fact = exp_fact

        if amp_tuple != 0:
            self.model.TabularAmplitude(data=amp_tuple, name=
                amp, timeSpan=STEP)
        else:
            amp_data = self.exp_amp(num_steps,exp_fact,start_point,end_point,step_name,offset)
            self.model.TabularAmplitude(data=amp_data, name=
                amp, timeSpan=STEP)

        #Increase minimum number of attempts from 5 to 30
        if maxU == maxU:
            self.model.steps[step_name].control.setValues(allowPropagation=OFF, 
                resetDefaultValues=OFF, timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0, 
                12.0, 25, 6.0, 3.0, 50.0))

        self.model.fieldOutputRequests['F-Output-1'].setValues(timePoint='TimePoints-ODB')

        self.model.HistoryOutputRequest(createStepName=step_name, name=
            'H-Output-2-'+step_name, rebar=EXCLUDE, region=   
            self.model.rootAssembly.sets[self.NameRef1], sectionPoints=DEFAULT, 
            variables=('U1','U2','U3','RF1','RF2','RF3'),timePoint='TimePoints-CSV')

        if not hasattr(self,'radius'):
            self.model.HistoryOutputRequest(createStepName=step_name, name=
                'H-Output-3-'+step_name, rebar=EXCLUDE, region=   
                self.model.rootAssembly.sets[self.NameRef2], sectionPoints=DEFAULT, 
                variables=('U1','U2','U3','RF1','RF2','RF3'),timePoint='TimePoints-CSV')
        
        self.model.fieldOutputRequests['F-Output-1'].setValues(variables=(
        'UR','S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'V', 'A', 'RF', 'CF', 'CSTRESS', 
        'CDISP', 'STH', 'SE','COORD','EVOL'),timePoint='TimePoints-ODB')

        self.num_csv_TP = num_csv_TP
        self.num_odb_TP = num_odb_TP
        amp_csv_TP = np.linspace(0,1,num_csv_TP)
        amp_odb_TP = np.linspace(0,1,num_odb_TP)

        time_csv_TP = num_csv_TP*[0]
        time_odb_TP = num_odb_TP*[0]
        if exp_fact == 1:
            for i,t in enumerate(np.linspace(0,1,num_csv_TP)):
                time_csv_TP[i] = ((float(t)),)
            for i,t in enumerate(np.linspace(0,1,num_odb_TP)):
                time_odb_TP[i] = ((float(t)),) 
        else:
            for i in range(num_csv_TP):
                    if i != 0:
                        time_csv_TP[i] = (((np.log10(amp_csv_TP[i]*((exp_fact)-1)+1))/np.log10(exp_fact)),)
                    else:
                        time_csv_TP[i] = ((0),)
            for i in range(num_odb_TP):
                    if i != 0:
                        time_odb_TP[i] = (((np.log10(amp_odb_TP[i]*((exp_fact)-1)+1))/np.log10(exp_fact)),)
                    else:
                        time_odb_TP[i] = ((0),)
        self.model.TimePoint(name='TimePoints-CSV', points=tuple(time_csv_TP))
        self.model.TimePoint(name='TimePoints-ODB', points=tuple(time_odb_TP))

    def exp_amp(self,num_timesteps,b,start_point,end_point,step_name,offset):
        # if step_name[-1] == '2':
        #     start_point=1

        def exp_fxn(x,b):
            if b ==1:
                return x
            else:
                return 1./(b-1) * (b**x - 1)

        all_t = np.linspace(start_point,end_point,num_timesteps)
        all_tuples = (num_timesteps)*[0]
        for i in range(num_timesteps):
            all_tuples[i] = (all_t[i],exp_fxn(all_t[i],b)+offset)
        return all_tuples

    def int_contact(self,intProp_name = 'IntProp-1',tangential=True):
        self.model.ContactProperty(intProp_name)
        if tangential:
            self.model.interactionProperties[intProp_name].TangentialBehavior(
                formulation=FRICTIONLESS)
        self.model.interactionProperties[intProp_name].NormalBehavior(
            allowSeparation=ON, constraintEnforcementMethod=AUGMENTED_LAGRANGE, 
            pressureOverclosure=HARD,clearanceAtZeroContactPressure=0,contactStiffnessScaleFactor=10)
        self.model.rootAssembly.Surface(name='m_Surf-3', side1Faces=
            self.model.rootAssembly.instances['cylinder-1'].faces)
        self.model.rootAssembly.Set(name='s_Set-125', nodes=
            self.model.rootAssembly.instances[self.name_cyl+'-1'].nodes)
        self.model.SurfaceToSurfaceContactStd(adjustMethod=OVERCLOSED, 
            clearanceRegion=None, createStepName='step-'+self.name_part, datumAxis=None, 
            initialClearance=OMIT, interactionProperty='IntProp-1', master=
            self.model.rootAssembly.surfaces['m_Surf-3'], name='Int-1', 
            slave=self.model.rootAssembly.sets['s_Set-125'], sliding=FINITE, 
            surfaceSmoothing=NONE, thickness=OFF, tied=OFF)
        mdb.models['Model-1'].interactions['Int-1'].setValues(adjustMethod=OVERCLOSED, 
            bondingSet=None, enforcement=NODE_TO_SURFACE, initialClearance=OMIT, 
            sliding=SMALL, smooth=0.2, supplementaryContact=SELECTIVE, 
            surfaceSmoothing=NONE, thickness=OFF, tied=OFF)
        
        mdb.models['Model-1'].interactionProperties['IntProp-1'].normalBehavior.setValues(
            allowSeparation=ON, clearanceAtZeroContactPressure=0.25, 
            constraintEnforcementMethod=PENALTY, contactStiffness=DEFAULT, 
            contactStiffnessScaleFactor=10.0, pressureOverclosure=HARD, 
            stiffnessBehavior=LINEAR)

    def analytical_cylinder(self,offset=0):
        # Generate an analytical cylinder to constrain PBC
        multiplier = 1
        height=multiplier*(np.abs(self.latticeVec1[2])+np.abs(self.latticeVec2[2]))
        #print('Height:{}'.format(height))
        tol = offset
        mdb.models['Model-1'].ConstrainedSketch(name='profile', sheetSize=200.0)
        mdb.models['Model-1'].sketches['profile'].ConstructionLine(point1=(0.0, 
            -200.0), point2=(0.0, 200.0))
        mdb.models['Model-1'].sketches['profile'].geometry.findAt((0.0, 0.0))
        mdb.models['Model-1'].sketches['profile'].FixedConstraint(entity=
            mdb.models['Model-1'].sketches['profile'].geometry.findAt((0.0, 0.0), 
            ))
        # print(self.radius+tol)
        # print(height)
        line = mdb.models['Model-1'].sketches['profile'].Line(point1=(self.radius+tol, 0), 
            point2=(self.radius+tol, -height))
        mdb.models['Model-1'].sketches['profile'].geometry.findAt((self.radius+tol, 0.0))
        mdb.models['Model-1'].sketches['profile'].VerticalConstraint(addUndoState=
            False, entity=line)
        mdb.models['Model-1'].Part(dimensionality=THREE_D, name='cylinder', type=
            ANALYTIC_RIGID_SURFACE)
        mdb.models['Model-1'].parts['cylinder'].AnalyticRigidSurfRevolve(sketch=
            mdb.models['Model-1'].sketches['profile'])
        del mdb.models['Model-1'].sketches['profile']

        self.model.rootAssembly.Instance(dependent=ON, name='cylinder-1', 
            part=self.model.parts['cylinder'])
        point_cyl = self.model.parts['cylinder'].faces[0].pointOn[0]
        self.model.parts['cylinder'].Surface(name='Surf-1', side1Faces=self.model.parts['cylinder'].faces.findAt((point_cyl, )))
        Cyl_RP = self.model.parts['cylinder'].ReferencePoint(point=(0,0,0))
        self.test = Cyl_RP
        self.model.rootAssembly.regenerate()
        self.model.rootAssembly.Set(name='Cyl_RP', referencePoints=(
            self.model.rootAssembly.instances['cylinder-1'].referencePoints[Cyl_RP.id], 
            ))

        mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(-90.0, 
            0.0, 0.0), axisPoint=(90.0, 0.0, 0.0), instanceList=('cylinder-1', ))
        self.model.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
            distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-Cyl', 
            region=self.model.rootAssembly.sets['Cyl_RP'], u1=SET, u2=SET, 
            u3=SET, ur1=SET, ur2=SET, ur3=SET)
            
    def PBC_BC_cyl2(self,DeformationF,amp,step,origin_fix=True,valfixed=OFF):
        REF1 = self.model.rootAssembly.sets[self.NameRef1]

        A = self.model.rootAssembly

        self.model.TabularAmplitude(data=((0.0, 1.0), (1.0, 1.0)), name=
            'Amp-2', smooth=SOLVER_DEFAULT, timeSpan=STEP)

        self.model.DisplacementBC(amplitude=amp, createStepName=step, 
            distributionType=UNIFORM, fieldName='', fixed=valfixed, name=
            'BC-REF1-'+step, region=REF1, 
            u1 = DeformationF[0][0],
            u2=  DeformationF[1][0], 
            u3=  DeformationF[1][1],  
            ur1=0,
            ur2=0,
            ur3=0)


        if origin_fix:
            self.model.DisplacementBC(amplitude=amp, createStepName=step, 
                distributionType=UNIFORM, fieldName='', fixed=valfixed, name=
                'BC-ENCASTRE', region=A.sets[self.name_cyl+'-1'+".origin"], 
                u1=0, 
                u2=0, 
                u3=UNSET, 
                ur1=UNSET,
                ur2=UNSET,
                ur3=UNSET)

    def PBC_BC(self,DeformationF,amp,step,origin_fix=True,valfixed=OFF):
        REF1 = self.model.rootAssembly.sets[self.NameRef1]
        REF2 = self.model.rootAssembly.sets[self.NameRef2]

        A = self.model.rootAssembly

        self.model.DisplacementBC(amplitude=amp, createStepName=step, 
            distributionType=UNIFORM, fieldName='', fixed=valfixed, localCsys=None, name=
            'BC-REF1-'+step, region=REF1, 
            u1=DeformationF[0][0], 
            u2=DeformationF[0][1], 
            u3=UNSET, 
            ur1=UNSET,
            ur2=UNSET,
            ur3=UNSET)

        self.model.DisplacementBC(amplitude=amp, createStepName=step, 
            distributionType=UNIFORM, fieldName='', fixed=valfixed, localCsys=None, name=
            'BC-REF2-'+step, region=REF2, 
            u1=DeformationF[1][0], 
            u2=DeformationF[1][1], 
            u3=UNSET, 
            ur1=UNSET,
            ur2=UNSET,
            ur3=UNSET)

        if origin_fix:
            self.model.DisplacementBC(amplitude=amp, createStepName=step, 
                distributionType=UNIFORM, fieldName='', fixed=valfixed, localCsys=None, name=
                'BC-ENCASTRE', region=A.sets[self.name_assem+".origin"], 
                u1=0, 
                u2=UNSET, 
                u3=0, 
                ur1=UNSET,
                ur2=UNSET,
                ur3=UNSET)

    def PBC_BC_Buckle(self,DeformationF,amp,step,origin_fix=True):
        REF1 = self.model.rootAssembly.sets[self.NameRef1]
        REF2 = self.model.rootAssembly.sets[self.NameRef2]

        A = self.model.rootAssembly

        self.model.DisplacementBC(createStepName=step, 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-REF1-'+step, region=REF1, 
            u1=DeformationF[0][0], 
            u2=DeformationF[0][1], 
            u3=UNSET, 
            ur1=UNSET,
            ur2=UNSET,
            ur3=UNSET)

        self.model.DisplacementBC(createStepName=step, 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-REF2-'+step, region=REF2, 
            u1=DeformationF[1][0], 
            u2=DeformationF[1][1], 
            u3=UNSET, 
            ur1=UNSET,
            ur2=UNSET,
            ur3=UNSET)

        if origin_fix:
            self.model.DisplacementBC(createStepName=step, 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-ENCASTRE', region=A.sets[self.name_assem+".origin"], 
                u1=0, 
                u2=0, 
                u3=0, 
                ur1=UNSET,
                ur2=UNSET,
                ur3=UNSET)
            
    def roughness_cyl(self):
        self.model.Part(dimensionality=THREE_D, name='RefPoint-3', type=
        DEFORMABLE_BODY)
        self.model.parts['RefPoint-3'].ReferencePoint(point=(self.radius, 0.0, 0.0))
        self.model.rootAssembly.Instance(dependent=ON, name='RefPoint-3-1', 
            part=self.model.parts['RefPoint-3'])
        self.model.rootAssembly.Set(name='RefPoint-3', referencePoints=(
            self.model.rootAssembly.instances['RefPoint-3-1'].referencePoints[1], 
            ))
        self.model.rootAssembly.Set(name='m_Set-46', referencePoints=(
            self.model.rootAssembly.instances['RefPoint-3-1'].referencePoints[1], 
            ))
        node_eq = ()
        n = self.model.rootAssembly.instances[self.name_cyl+'-1'].nodes
        for n_i,node in enumerate(self.model.rootAssembly.sets['nodes-all'].nodes):
            if node in self.model.rootAssembly.sets[self.name_cyl+'-1.'+'edge_bot'].nodes or node in self.model.rootAssembly.sets[self.name_cyl+'-1.'+'edge_left'].nodes:
                if node not in self.model.rootAssembly.sets[self.name_cyl+'-1.'+'edge_top'].nodes or node not in self.model.rootAssembly.sets[self.name_cyl+'-1.'+'edge_right'].nodes:
                    n_mask = n.sequenceFromLabels(labels=(node.label,))
                    self.model.rootAssembly.Set(nodes=(n_mask,), name='node_'+str(n_i))
                    node_eq = node_eq + ((2.0, 'node_'+str(n_i), 1, self.datum_id),)
            if node not in self.model.rootAssembly.sets[self.name_cyl+'-1.'+'edge_nodes'].nodes:
                    n_mask = n.sequenceFromLabels(labels=(node.label,))
                    self.model.rootAssembly.Set(nodes=(n_mask,), name='node_'+str(n_i))
                    node_eq = node_eq + ((1.0, 'node_'+str(n_i), 1, self.datum_id),)

        mdb.models['Model-1'].interactions['Int-1'].setValues(adjustMethod=OVERCLOSED, 
            bondingSet=None, datumAxis=None, enforcement=NODE_TO_SURFACE, 
            initialClearance=0, sliding=SMALL, smooth=0.2, supplementaryContact=
            SELECTIVE, thickness=OFF, tied=OFF)
        self.model.keywordBlock.synchVersions()
        self.model.keywordBlock.insert((len(self.model.keywordBlock.sieBlocks)-2),'\n*Node File, nset=nodes-all\nCOORD\n')
        self.model.keywordBlock.synchVersions()
        self.model.keywordBlock.insert((len(self.model.keywordBlock.sieBlocks)-2),'\n*Node File, nset=RefPoint-2\nU\n')

    def create_job(self,jobName,nCores,file_subroutine=False,path_subroutine='',initial_conditions=False,Gshear=0,MPC=False,Corr_Normal=False):
        # Regenerate assembly
        self.model.rootAssembly.regenerate()
        self.jobName = jobName
        self.file_subroutine = file_subroutine
        self.nCores = nCores
        # Create temporary job
        jobName_tmp = 'job_tmp'
        if self.file_subroutine == False:
            mdb.Job(
                name=jobName_tmp,
                numCpus=self.nCores, 
                numDomains=self.nCores, 
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
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
                # if linesToSkip > 0:
                #     linesToSkip -= 1
                #     continue

                if line.startswith('*Material,') and Gshear!=0:
                    newFile.write(line)
                    newFile.write('*Transverse Shear\n{},\n'.format(Gshear))
                elif line.startswith('0, MPC-corner') and MPC:
                    line = line[0:23]+'MPC-origin,'+line[23:-1]
                    newFile.write(line)
                else:
                    newFile.write(line)

            oldFile.close()
            newFile.close()

            # Create final job
            mdb.JobFromInputFile(
                name=self.jobName, 
                inputFileName=self.jobName+'.inp',
                numCpus=self.nCores,
                numDomains=self.nCores,
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
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
        else:
            mdb.Job(
                name=jobName_tmp,
                numCpus=self.nCores, 
                numDomains=self.nCores, 
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
                userSubroutine=path_subroutine,
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
            node_normal = False
            # linesToSkip = 0
            for line in oldFile:
                # if linesToSkip > 0:
                #     linesToSkip -= 1
                if node_normal:
                    if line.startswith('*Element'):
                        node_normal=False
                    else:
                        # try:
                        
                        # print(line_temp)
                        try:
                            line_temp = np.array(line.replace(',','').split(),dtype=float).tolist()
                            normal = np.array([line_temp[1],line_temp[2],0])/np.linalg.norm(np.array([line_temp[1],line_temp[2],0]))
                            # Convert NumPy array to list and add to numbers_list
                            line_temp.extend(normal.tolist())
                            line_temp[0] = int(line_temp[0])

                            # Convert the list back to a string with commas
                            line = ','.join(map(str, line_temp))
                            line = line+'\n'
 
                        except:
                            print('Zero normal encountered')
                            node_normal = False

                
                if line.startswith('*Orientation,'):
                    newFile.write(line)
                    # newFile.write(line[:-1] + ',local directions=2\n')
                elif line.startswith('*Transverse Shear'):
                    newFile.write('*Transverse Shear Stiffness\n')
                elif line =='*Node\n' and Corr_Normal:
                    node_normal = True
                    newFile.write(line)
                elif line.startswith('0, MPC-corner') and MPC:
                    line = line[0:22]+'MPC-origin,'+line[22:-1]+'\n'
                    newFile.write(line)
                elif line.startswith('0, MPC-node') and MPC:
                    line_temp = line.split(',')
                    node_temp = line_temp[1].split('-')
                    node_temp[2] = '1'
                    line_temp[-1] = '-'.join(node_temp)
                    line_temp.append('RefPoint-1')
                    line = ','.join(line_temp)+'\n'
                    newFile.write(line)
                else:
                    newFile.write(line)

            oldFile.close()
            newFile.close()

            # Create final job
            mdb.JobFromInputFile(
                name=self.jobName, 
                inputFileName=self.jobName+'.inp',
                numCpus=self.nCores,
                numDomains=self.nCores,
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
                userSubroutine=path_subroutine,
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

    def run_job(self,waitForCompletion=True,timeOut=6000,retainFiles=['.odb','.sta','.msg']):
        t1 = time.time()
        mdb.jobs[self.jobName].submit(
            consistencyChecking=OFF)
        print('Submitted: {}'.format(self.jobName))
        if waitForCompletion:
            try:
                print('Waiting: {}'.format(self.jobName))
                mdb.jobs[self.jobName].waitForCompletion(6000)
            except:
                print("Error: Job timed out")
        t2 = time.time()
        t_total = t2-t1
        time.sleep(2)
        
        # self.clean_up(retainFiles=retainFiles)
        return t_total

    def clean_up(self,retainFiles=['.odb','.sta','.msg']):
        extensions = ['.odb','.sta','.dat','.com','.ipm','.log','.prt','.sim',
            '.msg','.lck','.fil','inp','.stt','.mdl']
        for e in retainFiles:
            extensions.remove(e)
        # print(extensions)
        for e in extensions:
            try:
                os.remove(self.jobName + e)
            except:
                pass
        try:
            os.remove('job_tmp.inp')
        except:
            pass

    def output_U33_max(self,fileName):
        # Import the relavent data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)

        cc = 0
        U33max = np.zeros((len(odb.steps['step-'+self.name_part].frames)))
        for frame in odb.steps['step-'+self.name_part].frames:
            fieldU = frame.fieldOutputs['U']
            ndFieldU = fieldU.getSubset(region=odb.rootAssembly.instances[self.name_part.upper()+'-1'], position=NODAL)
            U33 = 0
            for value in ndFieldU.values:
                if U33 < abs(value.data[2]):
                    U33 = abs(value.data[2])
                    U33max[cc] = abs(value.data[2])
            cc = cc+1

        np.savetxt(str(fileName)+'.csv', U33max, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')

    def output_faces(self,file_name='faces.csv',cyl=False):
        # Import the relavent data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        if cyl:
            modified_rows = []
            n_elem = len(odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].elements)
            elem_array = np.zeros((n_elem,6))
            for i,elem in enumerate(odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].elements):
                np.append(elem_array,[[elem.label]+list(elem.connectivity)])
                elem_array[i,0] = elem.label
                elem_c = list(elem.connectivity)
                elem_array[i,1] = int(elem_c[0])-1
                elem_array[i,2] = int(elem_c[1])-1
                elem_array[i,3] = int(elem_c[2])-1
                if len(elem_c) == 4:
                    elem_array[i,4] = int(elem_c[3])-1
                    #modified_rows.append([elem_array[i,1],elem_array[i,2],elem_array[i,4]])
                else:
                    elem_array[i,4] = -1
                if 'WOVEN' in elem.sectionCategory.name:
                    elem_array[i,5] = 1
            elem_array = elem_array[elem_array[:, 0].argsort()]
        else:            
            n_elem = len(odb.rootAssembly.instances[self.name_part.upper()+'-1'].elements)
            elem_array = np.zeros((n_elem,5))
            for i,elem in enumerate(odb.rootAssembly.instances[self.name_part.upper()+'-1'].elements):
                np.append(elem_array,[[elem.label]+list(elem.connectivity)])
                elem_array[i,0] = elem.label
                elem_c = list(elem.connectivity)
                elem_array[i,1] = int(elem_c[0])-1
                elem_array[i,2] = int(elem_c[1])-1
                elem_array[i,3] = int(elem_c[2])-1
                if len(elem_c) == 4:
                    elem_array[i,4] = int(elem_c[3])-1
                else:
                    elem_array[i,4] = -1
                if 'WOVEN' in elem.sectionCategory.name:
                    elem_array[i,5] = 1
            elem_array = elem_array[elem_array[:, 0].argsort()]
        
        modified_rows = list(elem_array[:,1:5])
        modified_rows = []
        section_rows = []
        for elem in elem_array:
            if elem[4] == -1:
                modified_rows.append([elem[1],elem[2],elem[3]])
                section_rows.append([elem[5]])
            else:
                modified_rows.append([elem[1],elem[2],elem[4]])
                modified_rows.append([elem[2],elem[3],elem[4]])
                section_rows.append([elem[5]])
                section_rows.append([elem[5]])
        # output file paths
        output_csv_file = file_name


        # Write the modified data to the output CSV file
        with open(output_csv_file, mode='wb') as csv_file:
            csv_writer = csv.writer(csv_file)
            # Write the modified rows
            csv_writer.writerows(modified_rows)
        self.elem_connectivity = np.array(modified_rows)
        self.elem_section = np.array(section_rows)

        with open('sections.csv', mode='wb') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerows(self.elem_section)
    #     print('Post-processing UC is done. File saved.')
        print('Post-processing UC is done. File saved.')
        

    def output_PBC(self,file_name='PBC.csv',pos=-1,cyl=False,tol=0,offset=False):
        # Import the relavent data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        if pos == -1:
            nFrame_list = [len(odb.steps['step-'+self.name_part].frames) -1]
        elif pos == -2:
            nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
            print(nFrame_list)
        else:
            nFrame_list = [pos-1]
        V1_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_BR'].nodes[0].label
        V2_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_TL'].nodes[0].label
        n_PBC = (len(sampleUC.model.parts[sampleUC.name_part].sets['edgesAll'].nodes)-3)/2
        with open('integers.txt', 'w') as file:
            file.write("{} {}".format(V1_label, V2_label))
        
        elements = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].elements
        nodes = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodes
        
        for q,nFrame in enumerate(odb.steps['step-'+self.name_part].frames):
            modified_rows = []
            session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
            session.writeFieldReport(fileName='vertices_temp.csv', append=ON, 
                sortItem='Node Label', odb=odb, step=0, frame=nFrame, outputPosition=NODAL, 
                variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), (
                COMPONENT, 'COOR3'), )), ), stepFrame=SPECIFY)

            # Input and output file paths
            input_csv_file = 'vertices_temp.csv'
            output_csv_file = str(q)+'_'+file_name
            output_section_csv_file = 'sections.csv'
            fieldU = nFrame.fieldOutputs['COORD']
            fieldUR = nFrame.fieldOutputs['UR']
            # Corner Nodes
            node_origin = fieldU.getSubset(region=odb.rootAssembly.nodeSets['MPC-ORIGIN'], position=NODAL)
            node_corner_1 = fieldU.getSubset(region=odb.rootAssembly.nodeSets['MPC-CORNER-NODE-1'], position=NODAL)
            node_corner_2 = fieldU.getSubset(region=odb.rootAssembly.nodeSets['MPC-CORNER-NODE-2'], position=NODAL)
            node_corner_3 = fieldU.getSubset(region=odb.rootAssembly.nodeSets['MPC-CORNER-NODE-3'], position=NODAL)
            node_origin_R = fieldUR.getSubset(region=odb.rootAssembly.nodeSets['MPC-ORIGIN'], position=NODAL)
            node_corner_1_R = fieldUR.getSubset(region=odb.rootAssembly.nodeSets['MPC-CORNER-NODE-1'], position=NODAL)
            node_corner_2_R = fieldUR.getSubset(region=odb.rootAssembly.nodeSets['MPC-CORNER-NODE-2'], position=NODAL)
            node_corner_3_R = fieldUR.getSubset(region=odb.rootAssembly.nodeSets['MPC-CORNER-NODE-3'], position=NODAL)
            # return node_origin
            # Open the input CSV file for reading
            for [node,nodeR] in [[node_corner_1,node_corner_1_R],[node_corner_2,node_corner_2_R],[node_corner_3,node_corner_3_R]]:
                n1x, n1y, n1z = node_origin.values[0].data
                n1xR,n1yR,n1zR = node_origin_R.values[0].data
                n2x, n2y, n2z = node.values[0].data
                n2xR, n2yR, n2zR = nodeR.values[0].data
                modified_rows.append([n1x,n1y,n1z,n1xR,n1yR,n1zR,n2x,n2y,n2z,n2xR, n2yR, n2zR])
            for j in range(n_PBC):
                node_1 = fieldU.getSubset(region=odb.rootAssembly.nodeSets['MPC-NODE-1-'+str(j+1)], position=NODAL)
                node_2 = fieldU.getSubset(region=odb.rootAssembly.nodeSets['MPC-NODE-2-'+str(j+1)], position=NODAL)
                node_1_R = fieldUR.getSubset(region=odb.rootAssembly.nodeSets['MPC-NODE-1-'+str(j+1)], position=NODAL)
                node_2_R = fieldUR.getSubset(region=odb.rootAssembly.nodeSets['MPC-NODE-2-'+str(j+1)], position=NODAL)
                n1x, n1y, n1z = node_1.values[0].data
                n2x, n2y, n2z = node_2.values[0].data
                n1xR, n1yR, n1zR = node_1_R.values[0].data
                n2xR, n2yR, n2zR = node_2_R.values[0].data
                modified_rows.append([n1x,n1y,n1z,n1xR,n1yR,n1zR,n2x,n2y,n2z,n2xR,n2yR,n2zR])

            # Write the modified data to the output CSV file
            with open(output_csv_file, mode='wb') as csv_file:
                csv_writer = csv.writer(csv_file)
                
                # Write the header row
                #csv_writer.writerow(['Node Label', 'COORD-COOR1', 'COORD-COOR2', 'COORD-COOR3'])
                
                # Write the modified rows
                csv_writer.writerows(modified_rows)
            os.remove(input_csv_file)


    def output_P_D_double(self,file_name='dimple.csv',pressure=551.25e-6,pos=-1,cyl=False,tol=0,step=0):
        # Import the relavent data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        print(pos)
        if pos == -1:
            if step == 0:
                nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
                Frames = odb.steps['step-'+self.name_part].frames
            else:
                nFrame_list = list(range(0,len(odb.steps['step-pressure'].frames)))
                Frames = odb.steps['step-pressure'].frames
        elif pos == -2:
            nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
            print(nFrame_list)
        else:
            nFrame_list = [pos-1]
        origin_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['ORIGIN'].nodes[0].label
        V1_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_BR'].nodes[0].label
        V2_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_TL'].nodes[0].label
        with open('integers.txt', 'w') as file:
            file.write("{} {} {}".format(origin_label,V1_label, V2_label))
        print()
        #     # Create a list to store the modified rows
        modified_rows = []
        for nFrame,Frame in zip(nFrame_list,Frames):
            radius1 = []
            radius2 = []
            for Dimple in [1,2]:
                session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
                session.writeFieldReport(fileName='dimple_temp.csv', append=ON, 
                    sortItem='Node Label', odb=odb, step=step, frame=nFrame, outputPosition=NODAL, 
                    variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), (
                    COMPONENT, 'COOR3'), )), ), stepFrame=SPECIFY)
                fieldU = Frame.fieldOutputs['COORD']
                nodes = fieldU.getSubset(region=odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['DIMPLE_{}'.format(Dimple)], position=NODAL)
                # Input and output file paths
                input_csv_file = 'dimple_temp.csv'
                output_csv_file = file_name
                self.nodes = nodes

                
                # Iterate through the rows and extract the desired columns
                for node in nodes.values:
                    x,y,z = [float(node.dataDouble[0]), float(node.dataDouble[1]), float(node.dataDouble[2])]
                    r = np.sqrt(x**2+y**2)-self.radius+tol
                    theta = np.arctan2(y,x)
                    # modified_row = [pressure*(nFrame/nFrame_list[-1]),r]#[row[node_label_idx], row[coord1_idx], row[coord2_idx], row[coord3_idx]]
                    if Dimple ==1:
                        radius1.append(r)
                    else:
                        radius2.append(r)
                
            modified_rows.append([pressure*nFrame/nFrame_list[-1],np.max(radius1),np.max(radius2)])
            # Write the modified data to the output CSV file
        with open(output_csv_file, mode='wb') as csv_file:
            csv_writer = csv.writer(csv_file)
            
            # Write the header row
            #csv_writer.writerow(['Node Label', 'COORD-COOR1', 'COORD-COOR2', 'COORD-COOR3'])
            
            # Write the modified rows
            csv_writer.writerows(modified_rows)

            os.remove(input_csv_file)
        odb.close()
        print('Post-processing UC is done. File saved.')

    def output_P_D(self,file_name='dimple.csv',pressure=551.25e-6,pos=-1,cyl=False,tol=0,step=0):
        # Import the relavent data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        print(pos)
        if pos == -1:
            if step == 0:
                nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
                Frames = odb.steps['step-'+self.name_part].frames
            else:
                nFrame_list = list(range(0,len(odb.steps['step-pressure'].frames)))
                Frames = odb.steps['step-pressure'].frames
        elif pos == -2:
            nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
            print(nFrame_list)
        else:
            nFrame_list = [pos-1]
        V1_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_BR'].nodes[0].label
        V2_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_TL'].nodes[0].label
        with open('integers.txt', 'w') as file:
            file.write("{} {}".format(V1_label, V2_label))
        print()
        #     # Create a list to store the modified rows
        modified_rows = []
        for nFrame,Frame in zip(nFrame_list,Frames):
            radius1 = []
            radius2 = []
            for Dimple in [1,2]:
                session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
                session.writeFieldReport(fileName='dimple_temp.csv', append=ON, 
                    sortItem='Node Label', odb=odb, step=step, frame=nFrame, outputPosition=NODAL, 
                    variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), (
                    COMPONENT, 'COOR3'), )), ), stepFrame=SPECIFY)
                fieldU = Frame.fieldOutputs['COORD']
                nodes = fieldU.getSubset(region=odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['DIMPLE_{}'.format(Dimple)], position=NODAL)
                # Input and output file paths
                input_csv_file = 'dimple_temp.csv'
                output_csv_file = file_name
                self.nodes = nodes

                
                # Iterate through the rows and extract the desired columns
                for node in nodes.values:
                    x,y,z = [float(node.data[0]), float(node.data[1]), float(node.data[2])]
                    r = np.sqrt(x**2+y**2)-self.radius+tol
                    theta = np.arctan2(y,x)
                    # modified_row = [pressure*(nFrame/nFrame_list[-1]),r]#[row[node_label_idx], row[coord1_idx], row[coord2_idx], row[coord3_idx]]
                    if Dimple ==1:
                        radius1.append(r)
                    else:
                        radius2.append(r)
                
            modified_rows.append([pressure*nFrame/nFrame_list[-1],np.max(radius1),np.max(radius2)])
            # Write the modified data to the output CSV file
        with open(output_csv_file, mode='wb') as csv_file:
            csv_writer = csv.writer(csv_file)
            
            # Write the header row
            #csv_writer.writerow(['Node Label', 'COORD-COOR1', 'COORD-COOR2', 'COORD-COOR3'])
            
            # Write the modified rows
            csv_writer.writerows(modified_rows)

            os.remove(input_csv_file)
        odb.close()
        print('Post-processing UC is done. File saved.')


    def output_dimple(self,file_name='dimple.csv',pos=-1,cyl=False,tol=0,step=0):
        # Import the relavent data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        print(pos)
        if pos == -1:
            if step == 0:
                nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
                Frames = odb.steps['step-'+self.name_part].frames
            else:
                nFrame_list = list(range(0,len(odb.steps['step-pressure'].frames)))
                Frames = odb.steps['step-pressure'].frames
        elif pos == -2:
            nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
            print(nFrame_list)
        else:
            nFrame_list = [pos-1]
        V1_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_BR'].nodes[0].label
        V2_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_TL'].nodes[0].label
        with open('integers.txt', 'w') as file:
            file.write("{} {}".format(V1_label, V2_label))
        print()
        for nFrame,Frame in zip(nFrame_list,Frames):
            for Dimple in [1,2]:
                session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
                session.writeFieldReport(fileName='dimple_temp.csv', append=ON, 
                    sortItem='Node Label', odb=odb, step=step, frame=nFrame, outputPosition=NODAL, 
                    variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), (
                    COMPONENT, 'COOR3'), )), ), stepFrame=SPECIFY)
                fieldU = Frame.fieldOutputs['COORD']
                nodes = fieldU.getSubset(region=odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['DIMPLE_{}'.format(Dimple)], position=NODAL)
                # Input and output file paths
                input_csv_file = 'dimple_temp.csv'
                output_csv_file = str(nFrame)+'_'+str(Dimple)+'_'+file_name
                self.nodes = nodes

                #     # Create a list to store the modified rows
                modified_rows = []
                
                # Iterate through the rows and extract the desired columns
                for node in nodes.values:
                    x,y,z = [float(node.data[0]), float(node.data[1]), float(node.data[2])]
                    r = np.sqrt(x**2+y**2)-self.radius+tol
                    theta = np.arctan2(y,x)
                    modified_row = [theta*np.sqrt(x**2+y**2), z, r]#[row[node_label_idx], row[coord1_idx], row[coord2_idx], row[coord3_idx]]
                    modified_rows.append(modified_row)

                # Write the modified data to the output CSV file
                with open(output_csv_file, mode='wb') as csv_file:
                    csv_writer = csv.writer(csv_file)
                    
                    # Write the header row
                    #csv_writer.writerow(['Node Label', 'COORD-COOR1', 'COORD-COOR2', 'COORD-COOR3'])
                    
                    # Write the modified rows
                    csv_writer.writerows(modified_rows)

            os.remove(input_csv_file)
        odb.close()
        print('Post-processing UC is done. File saved.')

    def output_vertices(self,file_name='vertices.csv',pos=-1,cyl=False,tol=0):
        # Import the relavent data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        print(pos)
        if pos == -1:
            nFrame_list = [len(odb.steps['step-'+self.name_part].frames) -1]
        elif pos == -2:
            nFrame_list = list(range(0,len(odb.steps['step-'+self.name_part].frames)))
            print(nFrame_list)
        else:
            nFrame_list = [pos-1]
        origin_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['ORIGIN'].nodes[0].label
        V1_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_BR'].nodes[0].label
        V2_label = odb.rootAssembly.instances[self.name_cyl.upper()+'-1'].nodeSets['CN_TL'].nodes[0].label
        with open('integers.txt', 'w') as file:
            file.write("{} {} {}".format(origin_label,V1_label, V2_label))

        for nFrame in nFrame_list:
            session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
            session.writeFieldReport(fileName='vertices_temp.csv', append=ON, 
                sortItem='Node Label', odb=odb, step=0, frame=nFrame, outputPosition=NODAL, 
                variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), (
                COMPONENT, 'COOR3'), )), ), stepFrame=SPECIFY)
            
            # Input and output file paths
            input_csv_file = 'vertices_temp.csv'
            output_csv_file = str(nFrame)+'_'+file_name

            # Open the input CSV file for reading
            with open(input_csv_file, mode='r') as csv_file:
                csv_reader = csv.reader(csv_file)
                
                # # Skip the first two rows (headers)
                # next(csv_reader)
                # next(csv_reader)
                
                # Get the column indices for the desired columns
                header = next(csv_reader)
                node_label_idx = header.index('    Node Label')
                coord1_idx = header.index('   COORD-COOR1')
                coord2_idx = header.index('   COORD-COOR2')
                coord3_idx = header.index('   COORD-COOR3')

                # Create a list to store the modified rows
                modified_rows = []
                
                # Iterate through the rows and extract the desired columns
                for i,row in enumerate(csv_reader):
                        if cyl:
                            if i > 1:
                                    # try:
                                    # print(row)
                                    x,y,z = [float(row[coord1_idx]), float(row[coord2_idx]), float(row[coord3_idx])]
                                    # except:
                                    #     pass
                                    r = np.sqrt(x**2+y**2)-self.radius+tol
                                    # print(self.radius)
                                    theta = np.arctan2(y,x)
                                    modified_row = [theta*np.sqrt(x**2+y**2), z, r]#[row[node_label_idx], row[coord1_idx], row[coord2_idx], row[coord3_idx]]
                                    # modified_row = [r*np.cos(theta), z, r*np.sin(theta)]#[row[node_label_idx], row[coord1_idx], row[coord2_idx], row[coord3_idx]]
                                    modified_rows.append(modified_row)
                        else:
                            if i > 1:
                                modified_row = [row[coord1_idx], row[coord2_idx], row[coord3_idx]]#[row[node_label_idx], row[coord1_idx], row[coord2_idx], row[coord3_idx]]
                                modified_rows.append(modified_row)

            # Write the modified data to the output CSV file
            with open(output_csv_file, mode='wb') as csv_file:
                csv_writer = csv.writer(csv_file)
                
                # Write the header row
                #csv_writer.writerow(['Node Label', 'COORD-COOR1', 'COORD-COOR2', 'COORD-COOR3'])
                
                # Write the modified rows
                csv_writer.writerows(modified_rows)
            time.sleep(1)
            os.remove(input_csv_file)
        print('Post-processing UC is done. File saved.')

    def output_F_sig(self,fileName,PR,volume):
            # Import the relavent data
            session.journalOptions.setValues(replayGeometry=COORDINATE,
                                        recoverGeometry=COORDINATE )
                                        
            NaNValue = False
            #------------------------------------------------------------------------
            project = self.jobName
            # Post-processing
            odb = openOdb(project + '.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=odb)

            num=0

            F11_R=session.XYDataFromHistory(name='F11',odb=odb,outputVariableName='Spatial displacement: U1 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

            F21_R=session.XYDataFromHistory(name='F21',odb=odb,outputVariableName='Spatial displacement: U2 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

            S11_R=session.XYDataFromHistory(name='S11',odb=odb,outputVariableName='Reaction force: RF1 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

            S21_R=session.XYDataFromHistory(name='S21',odb=odb,outputVariableName='Reaction force: RF2 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

            F12_R=session.XYDataFromHistory(name='F12',odb=odb,outputVariableName='Spatial displacement: U1 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

            F22_R=session.XYDataFromHistory(name='F22',odb=odb,outputVariableName='Spatial displacement: U2 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

            S12_R=session.XYDataFromHistory(name='S12',odb=odb,outputVariableName='Reaction force: RF1 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

            S22_R=session.XYDataFromHistory(name='S22',odb=odb,outputVariableName='Reaction force: RF2 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

            RES=np.zeros((len(S22_R)-num,10))

            for i in range(num,len(S22_R)):
                if i == 0:
                    RES[i-num][0]=0
                    RES[i-num][1]=0
                    RES[i-num][2]=0
                    RES[i-num][3]=0
                    RES[i-num][4]=0
                    RES[i-num][5]=0
                    RES[i-num][6]=0
                    RES[i-num][7]=0
                    RES[i-num][8]=0
                    RES[i-num][9]=0
                else:
                    F11=F11_R[i][1]+1
                    F12=F12_R[i][1]
                    F21=F21_R[i][1]
                    F22=F22_R[i][1]+1
                    F33 = pow((F11*F22-F12*F21),(-PR/(1-PR)))
                    F = np.array([[F11,F12,0],[F21,F22,0],[0,0,F33]])
                    # Compute the singular value decomposition (SVD)
                    U, s, V = np.linalg.svd(F)
                    # Compute the polar decomposition of A
                    R = np.dot(U, V)
                    V = np.dot(F,np.transpose(R)) 
                    eL = self.matrix_log(V)
                    ep11 = eL[0][0]
                    ep12 = eL[0][1]
                    ep21 = eL[1][0]
                    ep22 = eL[1][1]
                    ep33 = eL[2][2]
                    S11=S11_R[i][1]
                    S12=S12_R[i][1] 
                    S21=S21_R[i][1]
                    S22=S22_R[i][1]
                    RES[i-num][0]=F11-1
                    RES[i-num][1]=F12
                    RES[i-num][2]=F21
                    RES[i-num][3]=F22-1
                    RES[i-num][4]=F33-1
                    volume_n = F33*pow((F11*F22-F12*F21),(-PR/(1-PR)))
                    P = np.array([[S11/volume,S12/volume,0],[S21/volume,S22/volume,0],[0,0,0]])
                    sig = (1/np.linalg.det(F))*(np.dot(F,np.transpose(P)))
                    RES[i-num][5]=sig[0][0]
                    RES[i-num][6]=sig[0][1]
                    RES[i-num][7]=sig[1][0]
                    RES[i-num][8]=sig[1][1]
                    RES[i-num][9]=0
                    if S11 != S11 or S12 != S12 or S21 != S21 or S22 != S22:
                        NaNValue = True

            np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
            odb.close()
            print('Post-processing UC is done. File saved.')
            return NaNValue,F22

    def output_ep_sig(self,fileName,PR,volume):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
                                    
        NaNValue = False
        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)

        num=0

        F11_R=session.XYDataFromHistory(name='F11',odb=odb,outputVariableName='Spatial displacement: U1 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        F21_R=session.XYDataFromHistory(name='F21',odb=odb,outputVariableName='Spatial displacement: U2 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        S11_R=session.XYDataFromHistory(name='S11',odb=odb,outputVariableName='Reaction force: RF1 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        S21_R=session.XYDataFromHistory(name='S21',odb=odb,outputVariableName='Reaction force: RF2 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        F12_R=session.XYDataFromHistory(name='F12',odb=odb,outputVariableName='Spatial displacement: U1 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        F22_R=session.XYDataFromHistory(name='F22',odb=odb,outputVariableName='Spatial displacement: U2 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        S12_R=session.XYDataFromHistory(name='S12',odb=odb,outputVariableName='Reaction force: RF1 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        S22_R=session.XYDataFromHistory(name='S22',odb=odb,outputVariableName='Reaction force: RF2 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        RES=np.zeros((len(S22_R)-num,10))

        for i in range(num,len(S22_R)):
            if i == 0:
                RES[i-num][0]=0
                RES[i-num][1]=0
                RES[i-num][2]=0
                RES[i-num][3]=0
                RES[i-num][4]=0
                RES[i-num][5]=0
                RES[i-num][6]=0
                RES[i-num][7]=0
                RES[i-num][8]=0
                RES[i-num][9]=0
            else:
                F11=F11_R[i][1]+1
                F12=F12_R[i][1]
                F21=F21_R[i][1]
                F22=F22_R[i][1]+1
                F33 = pow((F11*F22-F12*F21),(-PR/(1-PR)))
                F = np.array([[F11,F12,0],[F21,F22,0],[0,0,F33]])
                # Compute the singular value decomposition (SVD)
                U, s, V = np.linalg.svd(F)
                # Compute the polar decomposition of A
                R = np.dot(U, V)
                V = np.dot(F,np.transpose(R)) 
                eL = self.matrix_log(V)
                ep11 = eL[0][0]
                ep12 = eL[0][1]
                ep21 = eL[1][0]
                ep22 = eL[1][1]
                ep33 = eL[2][2]
                S11=S11_R[i][1]
                S12=S12_R[i][1] 
                S21=S21_R[i][1]
                S22=S22_R[i][1]
                RES[i-num][0]=ep11
                RES[i-num][1]=ep12
                RES[i-num][2]=ep21
                RES[i-num][3]=ep22
                RES[i-num][4]=ep33
                volume_n = F33*pow((F11*F22-F12*F21),(-PR/(1-PR)))
                P = np.array([[S11/volume,S12/volume,0],[S21/volume,S22/volume,0],[0,0,0]])
                sig = (1/np.linalg.det(F))*(np.dot(F,np.transpose(P)))
                RES[i-num][5]=sig[0][0]
                RES[i-num][6]=sig[0][1]
                RES[i-num][7]=sig[1][0]
                RES[i-num][8]=sig[1][1]
                RES[i-num][9]=0
                if S11 != S11 or S12 != S12 or S21 != S21 or S22 != S22:
                    NaNValue = True

        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return NaNValue,F22
    
    def matrix_log(self,matrix):
        """
        Compute the natural logarithm of a matrix using eigenvalue decomposition.
        
        Args:
        matrix (numpy.ndarray): The matrix to take the natural logarithm of.
        
        Returns:
        numpy.ndarray: The natural logarithm of the matrix.
        """
        # Compute the eigenvalue decomposition of the matrix
        eigvals, eigvecs = np.linalg.eig(matrix)
        
        # Take the natural logarithm of the eigenvalues
        log_eigvals = np.log(eigvals)
        
        # Construct the matrix logarithm from the logarithm of the eigenvalues and eigenvectors
        matrix_log = np.dot(eigvecs, np.dot(np.diag(log_eigvals), np.linalg.inv(eigvecs)))
        
        return matrix_log
    
    def output_E_D_cyl(self,fileName='depth'):
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

        F22_R=session.XYDataFromHistory(name='F31',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())
        
        F22_R=[]
        R_R = []
        cc_frame = 0
        for frame in odb.steps['step-'+self.name_part].frames:
            fieldCD = frame.fieldOutputs['COORD']
            ndFieldCD = fieldCD.getSubset(region=odb.rootAssembly.instances[self.name_cyl.upper()+'-1'], position=NODAL)
            fieldU = frame.fieldOutputs['U']
            ndFieldU = fieldU.getSubset(region=odb.rootAssembly.nodeSets['REFPOINT-1'], position=NODAL)
            D_max = 0
            cc = 0
            for i,value in enumerate(ndFieldCD.values):
                # print(value)
                depth = np.sqrt((value.data[0])**2+(value.data[1])**2)-(self.radius)
                if D_max < depth:
                    D_max = depth
                cc = cc+1
            cc_frame = cc_frame +1
            R_R.append(D_max)
            self.test = ndFieldU.values
            F22_R.append(ndFieldU.values[0].data[2])

        RES=np.zeros((len(F22_R)-num,2))
        for i in range(num,len(F22_R)):
            E=F22_R[i]
            R=R_R[i]
            RES[i-num][0]=E
            RES[i-num][1]=R


        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return


    def output_E_D_cyl_double(self,fileName='depth'):
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

        F22_R=session.XYDataFromHistory(name='F31',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())
        
        F22_R=[]
        R_R = []
        cc_frame = 0
        for frame in odb.steps['step-'+self.name_part].frames:
            fieldCD = frame.fieldOutputs['COORD']
            ndFieldCD = fieldCD.getSubset(region=odb.rootAssembly.instances[self.name_cyl.upper()+'-1'], position=NODAL)
            fieldU = frame.fieldOutputs['U']
            ndFieldU = fieldU.getSubset(region=odb.rootAssembly.nodeSets['REFPOINT-1'], position=NODAL)
            D_max = 0
            cc = 0
            for i,value in enumerate(ndFieldCD.values):
                # print(value)
                depth = np.sqrt((value.dataDouble[0])**2+(value.dataDouble[1])**2)-(self.radius)
                if D_max < depth:
                    D_max = depth
                cc = cc+1
            cc_frame = cc_frame +1
            R_R.append(D_max)
            self.test = ndFieldU.values
            F22_R.append(ndFieldU.values[0].dataDouble[2])

        RES=np.zeros((len(F22_R)-num,2))
        for i in range(num,len(F22_R)):
            E=F22_R[i]
            R=R_R[i]
            RES[i-num][0]=E
            RES[i-num][1]=R


        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return

    def output_E_R_cyl_double(self,fileName='roughness'):
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

        F22_R= []#session.XYDataFromHistory(name='F31',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())
        R_R = []
        cc_frame = 0
        for frame in odb.steps['step-'+self.name_part].frames:
            fieldU = frame.fieldOutputs['COORD']
            ndFieldU = fieldU.getSubset(region=odb.rootAssembly.instances[self.name_cyl.upper()+'-1'], position=NODAL)
            fieldRP = frame.fieldOutputs['U']
            ndFieldRP = fieldU.getSubset(region=odb.rootAssembly.nodeSets['REFPOINT-1'], position=NODAL)
            R_nodes = 0
            cc = 0
            for i,value in enumerate(ndFieldU.values):
                # print(value)
                R_nodes = R_nodes+(np.sqrt((value.dataDouble[0])**2+(value.dataDouble[1])**2)-(self.radius))
                cc = cc+1
            cc_frame = cc_frame +1
            F22_R.append(ndFieldRP.values[0].dataDouble[2])
            R_R.append(R_nodes/cc)

        RES=np.zeros((len(F22_R)-num,2))
        for i in range(num,len(F22_R)):
            E=F22_R[i]
            R=R_R[i]
            RES[i-num][0]=E
            RES[i-num][1]=R


        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return


    def output_E_R_cyl(self,fileName='roughness'):
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

        F22_R= []#session.XYDataFromHistory(name='F31',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())
        R_R = []
        cc_frame = 0
        for frame in odb.steps['step-'+self.name_part].frames:
            fieldU = frame.fieldOutputs['COORD']
            ndFieldU = fieldU.getSubset(region=odb.rootAssembly.instances[self.name_cyl.upper()+'-1'], position=NODAL)
            fieldRP = frame.fieldOutputs['U']
            ndFieldRP = fieldU.getSubset(region=odb.rootAssembly.nodeSets['REFPOINT-1'], position=NODAL)
            R_nodes = 0
            cc = 0
            for i,value in enumerate(ndFieldU.values):
                # print(value)
                R_nodes = R_nodes+(np.sqrt((value.data[0])**2+(value.data[1])**2)-(self.radius))
                cc = cc+1
            cc_frame = cc_frame +1
            F22_R.append(ndFieldRP.values[0].data[2])
            R_R.append(R_nodes/cc)

        RES=np.zeros((len(F22_R)-num,2))
        for i in range(num,len(F22_R)):
            E=F22_R[i]
            R=R_R[i]
            RES[i-num][0]=E
            RES[i-num][1]=R


        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return

    def output_E_V_cyl_double(self,fileName):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
                                    
        NaNValue = False
        ################################################

        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        # frame of pre-loading step
        num=0

        F22_R=session.XYDataFromHistory(name='F22',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        ALLKE_R=session.XYDataFromHistory(name='ALLKE',odb=odb,outputVariableName='Kinetic energy: ALLKE for Whole Model')
        RES=np.zeros((len(ALLKE_R)-num,2))

        for i in range(num,len(ALLKE_R)):
            F22=F22_R[i][1]
            ALLKE=ALLKE_R[i][1]

            RES[i-num][0]=F22
            RES[i-num][1]=ALLKE


        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return NaNValue,F22
    
    def output_E_V_cyl(self,fileName):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
                                    
        NaNValue = False
        ################################################

        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        # frame of pre-loading step
        num=0

        F22_R=session.XYDataFromHistory(name='F22',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        ALLKE_R=session.XYDataFromHistory(name='ALLKE',odb=odb,outputVariableName='Kinetic energy: ALLKE for Whole Model')
        RES=np.zeros((len(ALLKE_R)-num,2))

        for i in range(num,len(ALLKE_R)):
            F22=F22_R[i][1]
            ALLKE=ALLKE_R[i][1]

            RES[i-num][0]=F22
            RES[i-num][1]=ALLKE


        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return NaNValue,F22

    def output_E_P_cyl(self,fileName,PR,volume):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
                                    
        NaNValue = False
        ################################################
        a_1 = [self.xrec[1]-self.xrec[0], self.yrec[1]-self.yrec[0]]
        b_1 = [1, 0]
        theta_1 = np.arccos((np.dot(a_1, b_1)) / (np.linalg.norm(a_1) * np.linalg.norm(b_1)))
        ac_1 = a_1
        ac_1 = np.append(ac_1, 0)
        bc_1 = b_1
        bc_1 = np.append(bc_1, 0)
        crossab_1 = np.cross(ac_1, bc_1)
        if crossab_1[2] < 0:
            theta_1 = 2*np.pi - theta_1
        a_2 = [self.xrec[3]-self.xrec[0], self.yrec[3]-self.yrec[0]]
        b_2 = [1, 0]
        theta_2 = np.arccos((np.dot(a_2, b_2)) / (np.linalg.norm(a_2) * np.linalg.norm(b_2)))
        # print(theta_1)
        # print(theta_2)
        ac_2 = a_2
        ac_2 = np.append(ac_2, 0)
        bc_2 = b_2
        bc_2 = np.append(bc_2, 0)
        crossab_2 = np.cross(ac_2, bc_2)
        if crossab_2[2] < 0:
            theta_2 = 2*np.pi - theta_2
            # print(theta_2)
        d1 = np.sqrt((self.xrec[1]-self.xrec[0])**2+(self.yrec[1]-self.yrec[0])**2) #top/bot length
        d2 = np.sqrt((self.xrec[3]-self.xrec[0])**2+(self.yrec[3]-self.yrec[0])**2) #left/right length
        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        # frame of pre-loading step
        num=0


        F22_R=session.XYDataFromHistory(name='F31',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())


        S22_R=session.XYDataFromHistory(name='S31',odb=odb,outputVariableName='Reaction force: RF3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())


        RES=np.zeros((len(S22_R)-num,2))

        for i in range(num,len(S22_R)):
            F22=F22_R[i][1]
            S22=S22_R[i][1]

            RES[i-num][0]=F22
            RES[i-num][1]=S22/volume

            if S22 != S22:
                NaNValue = True

        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return NaNValue,F22

    def output_E_P_cyl_double(self,fileName,PR,volume):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
                                    
        NaNValue = False
        ################################################
        a_1 = [self.xrec[1]-self.xrec[0], self.yrec[1]-self.yrec[0]]
        b_1 = [1, 0]
        theta_1 = np.arccos((np.dot(a_1, b_1)) / (np.linalg.norm(a_1) * np.linalg.norm(b_1)))
        ac_1 = a_1
        ac_1 = np.append(ac_1, 0)
        bc_1 = b_1
        bc_1 = np.append(bc_1, 0)
        crossab_1 = np.cross(ac_1, bc_1)
        if crossab_1[2] < 0:
            theta_1 = 2*np.pi - theta_1
        a_2 = [self.xrec[3]-self.xrec[0], self.yrec[3]-self.yrec[0]]
        b_2 = [1, 0]
        theta_2 = np.arccos((np.dot(a_2, b_2)) / (np.linalg.norm(a_2) * np.linalg.norm(b_2)))
        # print(theta_1)
        # print(theta_2)
        ac_2 = a_2
        ac_2 = np.append(ac_2, 0)
        bc_2 = b_2
        bc_2 = np.append(bc_2, 0)
        crossab_2 = np.cross(ac_2, bc_2)
        if crossab_2[2] < 0:
            theta_2 = 2*np.pi - theta_2
            # print(theta_2)
        d1 = np.sqrt((self.xrec[1]-self.xrec[0])**2+(self.yrec[1]-self.yrec[0])**2) #top/bot length
        d2 = np.sqrt((self.xrec[3]-self.xrec[0])**2+(self.yrec[3]-self.yrec[0])**2) #left/right length
        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        # frame of pre-loading step
        num=0

        F22_R=session.XYDataFromHistory(name='F31',odb=odb,outputVariableName='Spatial displacement: U3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())
        S22_R=session.XYDataFromHistory(name='S31',odb=odb,outputVariableName='Reaction force: RF3 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())


        RES=np.zeros((len(S22_R)-num,2))

        for i in range(num,len(S22_R)):
            F22=F22_R[i][1]
            S22=S22_R[i][1]
            RES[i-num][0]=F22
            RES[i-num][1]=S22/volume
            if S22 != S22:
                NaNValue = True

        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return NaNValue,F22
        
    def read_csv(self,fileName):
        # file = open("data.csv")

        arr = np.loadtxt(fileName,delimiter=",", dtype=str)
        disp = []
        force = []
        for point in arr:
            disp.append(float(point[0]))
            force.append(float(point[1]))
        return disp,force



    def output_F_P(self,fileName,PR,volume):
        # Import the relavent data
        session.journalOptions.setValues(replayGeometry=COORDINATE,
                                    recoverGeometry=COORDINATE )
                                    
        NaNValue = False
        ################################################
        a_1 = [self.xrec[1]-self.xrec[0], self.yrec[1]-self.yrec[0]]
        b_1 = [1, 0]
        theta_1 = np.arccos((np.dot(a_1, b_1)) / (np.linalg.norm(a_1) * np.linalg.norm(b_1)))
        ac_1 = a_1
        ac_1 = np.append(ac_1, 0)
        bc_1 = b_1
        bc_1 = np.append(bc_1, 0)
        crossab_1 = np.cross(ac_1, bc_1)
        if crossab_1[2] < 0:
            theta_1 = 2*np.pi - theta_1
        a_2 = [self.xrec[3]-self.xrec[0], self.yrec[3]-self.yrec[0]]
        b_2 = [1, 0]
        theta_2 = np.arccos((np.dot(a_2, b_2)) / (np.linalg.norm(a_2) * np.linalg.norm(b_2)))
        # print(theta_1)
        # print(theta_2)
        ac_2 = a_2
        ac_2 = np.append(ac_2, 0)
        bc_2 = b_2
        bc_2 = np.append(bc_2, 0)
        crossab_2 = np.cross(ac_2, bc_2)
        if crossab_2[2] < 0:
            theta_2 = 2*np.pi - theta_2
            # print(theta_2)
        d1 = np.sqrt((self.xrec[1]-self.xrec[0])**2+(self.yrec[1]-self.yrec[0])**2) #top/bot length
        d2 = np.sqrt((self.xrec[3]-self.xrec[0])**2+(self.yrec[3]-self.yrec[0])**2) #left/right length
        #------------------------------------------------------------------------
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        # frame of pre-loading step
        num=0

        F11_R=session.XYDataFromHistory(name='F11',odb=odb,outputVariableName='Spatial displacement: U1 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        F21_R=session.XYDataFromHistory(name='F21',odb=odb,outputVariableName='Spatial displacement: U2 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        S11_R=session.XYDataFromHistory(name='S11',odb=odb,outputVariableName='Reaction force: RF1 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        S21_R=session.XYDataFromHistory(name='S21',odb=odb,outputVariableName='Reaction force: RF2 PI: '+self.NameRef1.upper()+' Node 1 in NSET '+self.NameRef1.upper())

        F12_R=session.XYDataFromHistory(name='F12',odb=odb,outputVariableName='Spatial displacement: U1 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        F22_R=session.XYDataFromHistory(name='F22',odb=odb,outputVariableName='Spatial displacement: U2 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        S12_R=session.XYDataFromHistory(name='S12',odb=odb,outputVariableName='Reaction force: RF1 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        S22_R=session.XYDataFromHistory(name='S22',odb=odb,outputVariableName='Reaction force: RF2 PI: '+self.NameRef2.upper()+' Node 1 in NSET '+self.NameRef2.upper())

        RES=np.zeros((len(S22_R)-num,10))

        for i in range(num,len(S22_R)):
            F11=F11_R[i][1]+1
            F12=F12_R[i][1]
            F21=F21_R[i][1]
            F22=F22_R[i][1]+1
            F33 = pow((F11*F22-F12*F21),(-PR/(1-PR)))
            S11=S11_R[i][1]
            S12=S12_R[i][1] 
            S21=S21_R[i][1]
            S22=S22_R[i][1]
            RES[i-num][0]=F11
            RES[i-num][1]=F12
            RES[i-num][2]=F21
            RES[i-num][3]=F22
            RES[i-num][4]=F33
            RES[i-num][5]=S11/volume
            RES[i-num][6]=S12/volume
            RES[i-num][7]=S21/volume
            RES[i-num][8]=S22/volume
            RES[i-num][9]=0
            if S11 != S11 or S12 != S12 or S21 != S21 or S22 != S22:
                NaNValue = True

        np.savetxt(str(fileName)+'.csv', RES, fmt='%6.8f',delimiter=',')
        odb.close()
        print('Post-processing UC is done. File saved.')
        return NaNValue,F22

        
    def inLine(self,A,B,C):
        # Initial Edge cases Check
        if (B[0]-A[0]) != 0:
            slopeAB = (B[1]-A[1])/(B[0]-A[0])
        elif (C[0]-A[0]) != 0:
            print('inLine Edge Case #1: Line AB slope infinite, C lies outside that line.')
            return False
        else:
            print('inLine Edge Case #2: All three points lie along line of infinite slope.')
            return True

        if  (C[0]-A[0]) != 0:
            slopeAC = (C[1]-A[1])/(C[0]-A[0])   # Slope of point AC
        else:
            print('inLine Edge Case #3: Line AB has finite, Line AC has infinite slope.')
            return False

        tol = 1e-2  # Tolerance of the inLine check
        if abs(slopeAB-slopeAC) < tol:
            return True
        else:
            return False

    def output_Poisson(self, fileName,flip=False):
        # Import the relevant data
        project = self.jobName
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)

        # Initialize area sums
        initial_area_sum = 0

        # Get the first and last frames
        first_frame = odb.steps['step-'+self.name_part].frames[0]

        # Function to calculate the area of a triangular element
        def calculate_triangle_area(node1, node2, node3):
            x1, y1, z1 = node1.data[0], node1.data[1], node1.data[2]
            x2, y2, z2 = node2.data[0], node2.data[1], node2.data[2]
            x3, y3, z3 = node3.data[0], node3.data[1], node3.data[2]
            vector1 = np.array([x2 - x1, y2 - y1, z2 - z1])
            vector2 = np.array([x3 - x1, y3 - y1, z3 - z1])
            cross_product = np.cross(vector1, vector2)
            area = 0.5 * np.linalg.norm(cross_product)
            if area > 80:
                print('Area is too large')
                print(node1.data)
                print(node2.data)
                print(node3.data)       
            return area
              
        # Iterate through each element and calculate the area for the first frame
        fieldCOORD = first_frame.fieldOutputs['COORD']
        ndFieldCOORD = fieldCOORD.getSubset(region=odb.rootAssembly.instances[self.name_part.upper()+'-1'], position=NODAL)
        for elem in odb.rootAssembly.instances[self.name_part.upper()+'-1'].elements:
            self.elems = elem
            connectivity = elem.connectivity
            if len(connectivity) == 3:  # Triangular element

                self.first_frame = ndFieldCOORD
                node1 = ndFieldCOORD.values[connectivity[0]-1]
                node2 = ndFieldCOORD.values[connectivity[1]-1]
                node3 = ndFieldCOORD.values[connectivity[2]-1]
                initial_area_sum += calculate_triangle_area(node1, node2, node3)
            elif len(connectivity) == 4:  # Quadrilateral element, split into two triangles
                node1 = ndFieldCOORD.values[connectivity[0]-1]
                node2 = ndFieldCOORD.values[connectivity[1]-1]
                node3 = ndFieldCOORD.values[connectivity[2]-1]
                node4 = ndFieldCOORD.values[connectivity[3]-1]
                initial_area_sum += calculate_triangle_area(node1, node2, node3)
                initial_area_sum += calculate_triangle_area(node1, node3, node4)

        poisson_list = []
        poisson_list_trad = []
        # Generate frame list:
        strain_list_RP1 = []
        strain_list_RP2 = []

        for i,frame in enumerate(odb.steps['step-'+self.name_part].frames):
            if i != 0:
                # Initialize area sums
                final_area_sum = 0
                # Iterate through each element and calculate the area for the last frame
                fieldCOORD = frame.fieldOutputs['COORD']
                ndFieldCOORD = fieldCOORD.getSubset(region=odb.rootAssembly.instances[self.name_part.upper()+'-1'], position=NODAL)

                fieldU = frame.fieldOutputs['U']
                ndFieldU1 = fieldU.getSubset(region=odb.rootAssembly.nodeSets['REFPOINT-1'], position=NODAL)
                ndFieldU2 = fieldU.getSubset(region=odb.rootAssembly.nodeSets['REFPOINT-2'], position=NODAL)
                for elem in odb.rootAssembly.instances[self.name_part.upper()+'-1'].elements:
                    connectivity = elem.connectivity
                    if len(connectivity) == 3:  # Triangular element
                        node1 = ndFieldCOORD.values[connectivity[0]-1]
                        node2 = ndFieldCOORD.values[connectivity[1]-1]
                        node3 = ndFieldCOORD.values[connectivity[2]-1]
                        final_area_sum += calculate_triangle_area(node1, node2, node3)
                    elif len(connectivity) == 4:  # Quadrilateral element, split into two triangles
                        node1 = ndFieldCOORD.values[connectivity[0]-1]
                        node2 = ndFieldCOORD.values[connectivity[1]-1]
                        node3 = ndFieldCOORD.values[connectivity[2]-1]
                        node4 = ndFieldCOORD.values[connectivity[3]-1]
                        final_area_sum += calculate_triangle_area(node1, node2, node3)
                        final_area_sum += calculate_triangle_area(node1, node3, node4)

                # Compare the initial and final area sums using new poisson ratio method
                area_change = final_area_sum - initial_area_sum 
                if flip:
                    e2 = ndFieldU1.values[0].data[0]
                else:
                    e2 = ndFieldU2.values[0].data[1]
                e1 = (final_area_sum/(initial_area_sum*(1+e2)))-1
                try:
                    poisson_ratio = -e1/e2
                except:
                    poisson_ratio = np.nan


                # Compare the poisson ratio method with the traditional method
                if flip:
                    e1 = ndFieldU2.values[0].data[1]
                    e2 = ndFieldU1.values[0].data[0]
                else:
                    e1 = ndFieldU1.values[0].data[0]
                    e2 = ndFieldU2.values[0].data[1]
                    
                try:
                    poisson_ratio_trad = -e1/e2
                except:
                    poisson_ratio_trad = np.nan

                # Append the Poisson ratio to the list
                poisson_list.append(poisson_ratio)
                poisson_list_trad.append(poisson_ratio_trad) 
                strain_list_RP2.append(e2) 

        # Save the Poisson ratio and strain list to a CSV file using numpy
        min_length = min(len(strain_list_RP2), len(poisson_list), len(poisson_list_trad))
        data = np.column_stack((strain_list_RP2[:min_length], poisson_list[:min_length], poisson_list_trad[:min_length]))
        np.savetxt(str(fileName) + '.csv', data, delimiter=',', header='Strain,Poisson_Ratio,Poisson_Ratio_Trad', comments='')
        odb.close()
