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

#   Basic abaqus class that builds a metamaterial finite size
#   sample based on a provided unit cell. 
# 
#   Note: For Abaqus use only, this class is not intended for standard python use 
#
#   INPUTS:
#   model - abaqus model variables
#   UC_list - list of UC classes that will be used in the design
#   UC_boundingbox - list of bounding boxes for the printing sample

class generate_finite:
    def __init__(self,model,UC_list,UC_boundingbox):
        # Unpack Variables
        self.model = model
        self.UC_list = UC_list
        self.UC_boundingbox = UC_boundingbox

    def generate_grid(self):
        UC_grid_temp = []
        UC_grid =[]
        for box in self.UC_boundingbox:
            L_oversize = hypot(box[3],box[4])
            UC = self.UC_list[box[0]]
            for i in range(int(-2*L_oversize/UC.L_base),int(2*L_oversize/UC.L_base),1):
                for j in range(int(-2*L_oversize/UC.pattern_hypot),int(2*L_oversize/UC.pattern_hypot),1):
                    xtemp = box[5]+UC.xlatticeVec[0]*i+UC.ylatticeVec[0]*j
                    ytemp = box[6]+UC.xlatticeVec[1]*i+UC.ylatticeVec[1]*j
                    ztemp = 0
                    UC_grid_temp.append([box[0],xtemp,ytemp,ztemp])
            # Check if all the points are inside the bounding box
            pt = np.zeros(2)
            rect = np.zeros(4)
            tol = hypot(UC.L_base,UC.pattern_hypot)*2.2
            rect[0] = box[1]-tol
            rect[1] = box[2]-tol
            rect[2] = box[3]+tol
            rect[3] = box[4]+tol
            cc = 0
            for point in UC_grid_temp:
                cc= cc+1
                pt[0] = point[1]
                pt[1] = point[2]
                if self.inRectangle(pt,rect[0],rect[1],rect[2],rect[3]):
                    UC_grid.append(point)

        return UC_grid

    def inRectangle(self,point,x0,y0,x1,y1):
        if x0 <= point[0] and x1 >= point[0] and y0 <= point[1] and y1 >= point[1]:
            return True
        else:
            return False



    def geometry_uniaxial(self,width,height,name_assem,dimension=THREE_D):
        # Store variables
        self.height = height
        self.width = width
        self.name_assem = name_assem
        self.dimension=dimension
        self.UC_grid = self.generate_grid()        
        #counter clock-wise points of the unit cell fabric
        self.xrec = [0, self.width, self.width,0]
        self.yrec = [0, 0, self.height, self.height]

        ################################################
        # Creating base fabric
        ################################################
        s = self.model.ConstrainedSketch(name='base_textile', sheetSize=50)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        bot = s.Line(point1=(self.xrec[0], self.yrec[0]), point2=(self.xrec[1], self.yrec[1]))
        right = s.Line(point1=(self.xrec[1], self.yrec[1]), point2=(self.xrec[2], self.yrec[2]))
        top = s.Line(point1=(self.xrec[2], self.yrec[2]), point2=(self.xrec[3], self.yrec[3]))
        left = s.Line(point1=(self.xrec[3], self.yrec[3]), point2=(self.xrec[0], self.yrec[0]))
        p = self.model.Part(name=self.name_assem, dimensionality=self.dimension, type=DEFORMABLE_BODY)
        self.model.parts[self.name_assem].BaseShell(sketch=s)
        s.unsetPrimaryObject()

        del self.model.sketches['base_textile']
        e = self.model.parts[self.name_assem].edges
        self.model.parts[self.name_assem].Set(edges=e,  name='edgesAll')
        vectorList = ['edge_bot','edge_right','edge_top','edge_left']
        for edge in self.model.parts[self.name_assem].edges:
            edgeList = []
            coords = edge.pointOn
            C = [coords[0][0],coords[0][1]]
            for i in range(4):
                A = [self.xrec[i],self.yrec[i]]
                if i == 3:
                    B = [self.xrec[0],self.yrec[0]]
                else:
                    B = [self.xrec[i+1],self.yrec[i+1]]
                if self.isBetween(A, B, C):
                    edgeList.append(e.findAt((edge.pointOn[0],), ))
                    self.model.parts[self.name_assem].Set(edges=edgeList,  name=vectorList[i])

        ################################################
        # CCreate cutter
        ################################################
        s = self.model.ConstrainedSketch(name='cut_sketch', sheetSize=50)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        P = 4
        bot = s.Line(point1=(self.xrec[0], self.yrec[0]), point2=(self.xrec[1], self.yrec[1]))
        right = s.Line(point1=(self.xrec[1], self.yrec[1]), point2=(self.xrec[2], self.yrec[2]))
        top = s.Line(point1=(self.xrec[2], self.yrec[2]), point2=(self.xrec[3], self.yrec[3]))
        left = s.Line(point1=(self.xrec[3], self.yrec[3]), point2=(self.xrec[0], self.yrec[0]))
        bot = s.Line(point1=(self.xrec[0]-L_finite*P, self.yrec[0]-H_finite*P), point2=(self.xrec[1]+L_finite*P, self.yrec[1]-H_finite*P))
        right = s.Line(point1=(self.xrec[1]+L_finite*P, self.yrec[1]-H_finite*P), point2=(self.xrec[2]+L_finite*P, self.yrec[2]+H_finite*P))
        top = s.Line(point1=(self.xrec[2]+L_finite*P, self.yrec[2]+H_finite*P), point2=(self.xrec[3]-L_finite*P, self.yrec[3]+H_finite*P))
        left = s.Line(point1=(self.xrec[3]-L_finite*P, self.yrec[3]+H_finite*P), point2=(self.xrec[0]-L_finite*P, self.yrec[0]-H_finite*P))
        p = self.model.Part(name='cut_temp', dimensionality=self.dimension, type=DEFORMABLE_BODY)
        self.model.parts['cut_temp'].BaseShell(sketch=s)
        s.unsetPrimaryObject()

        del self.model.sketches['cut_sketch']
        s = self.model.ConstrainedSketch(name='partition', sheetSize=50)
        UC_assem_list = []

        for i,indx in enumerate(self.UC_grid):
            pos = [indx[1]+origin_base[0],indx[2]+origin_base[1],indx[3]]
            UC = self.UC_list[indx[0]]
            UC_assem = self.duplicate_UC(UC,pos,i)
            UC_assem_list.append(UC_assem)

        tmerge = np.zeros(len(self.model.rootAssembly.instances.keys()),dtype=dict)

        for j in range(len(self.model.rootAssembly.instances.keys())):
            tmerge[j] = self.model.rootAssembly.instances[self.model.rootAssembly.instances.keys()[j]]

        self.model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
            instances=list(tmerge), 
            name='merge_temp', originalInstances=DELETE)
        self.model.rootAssembly.Instance(dependent=ON, name='cut_temp-1', 
            part=self.model.parts['cut_temp'])
        self.model.rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
            self.model.rootAssembly.instances['cut_temp-1'], ), 
            instanceToBeCut=self.model.rootAssembly.instances['merge_temp-1'], 
            name=self.name_assem, originalInstances=DELETE)

        del self.model.parts['cut_temp']
        del self.model.parts['merge_temp']
        self.TM_remove_faces(self.model.parts[self.name_assem])


    def geometry_biaxial(self,width,height,name_assem,dimension=THREE_D):
        # Store variables
        self.height = height
        self.width = width
        self.name_assem = name_assem
        self.dimension=dimension
        self.UC_grid = self.generate_grid()        
        #counter clock-wise points of the unit cell fabric
        self.xrec = [0, self.width, self.width,0]
        self.yrec = [0, 0, self.height, self.height]

        ################################################
        # Creating base fabric
        ################################################
        s = self.model.ConstrainedSketch(name='base_textile', sheetSize=50)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        bot = s.Line(point1=(self.xrec[0], self.yrec[0]), point2=(self.xrec[1], self.yrec[1]))
        right = s.Line(point1=(self.xrec[1], self.yrec[1]), point2=(self.xrec[2], self.yrec[2]))
        top = s.Line(point1=(self.xrec[2], self.yrec[2]), point2=(self.xrec[3], self.yrec[3]))
        left = s.Line(point1=(self.xrec[3], self.yrec[3]), point2=(self.xrec[0], self.yrec[0]))
        p = self.model.Part(name=self.name_assem, dimensionality=self.dimension, type=DEFORMABLE_BODY)
        self.model.parts[self.name_assem].BaseShell(sketch=s)
        s.unsetPrimaryObject()

        del self.model.sketches['base_textile']
        e = self.model.parts[self.name_assem].edges
        self.model.parts[self.name_assem].Set(edges=e,  name='edgesAll')
        vectorList = ['edge_bot','edge_right','edge_top','edge_left']
        for edge in self.model.parts[self.name_assem].edges:
            edgeList = []
            coords = edge.pointOn
            C = [coords[0][0],coords[0][1]]
            for i in range(4):
                A = [self.xrec[i],self.yrec[i]]
                if i == 3:
                    B = [self.xrec[0],self.yrec[0]]
                else:
                    B = [self.xrec[i+1],self.yrec[i+1]]
                if self.isBetween(A, B, C):
                    edgeList.append(e.findAt((edge.pointOn[0],), ))
                    self.model.parts[self.name_assem].Set(edges=edgeList,  name=vectorList[i])

        ################################################
        # CCreate cutter
        ################################################
        s = self.model.ConstrainedSketch(name='cut_sketch', sheetSize=50)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        P = 4
        bot = s.Line(point1=(self.xrec[0], self.yrec[0]), point2=(self.xrec[1], self.yrec[1]))
        right = s.Line(point1=(self.xrec[1], self.yrec[1]), point2=(self.xrec[2], self.yrec[2]))
        top = s.Line(point1=(self.xrec[2], self.yrec[2]), point2=(self.xrec[3], self.yrec[3]))
        left = s.Line(point1=(self.xrec[3], self.yrec[3]), point2=(self.xrec[0], self.yrec[0]))
        bot = s.Line(point1=(self.xrec[0]-L_finite*P, self.yrec[0]-H_finite*P), point2=(self.xrec[1]+L_finite*P, self.yrec[1]-H_finite*P))
        right = s.Line(point1=(self.xrec[1]+L_finite*P, self.yrec[1]-H_finite*P), point2=(self.xrec[2]+L_finite*P, self.yrec[2]+H_finite*P))
        top = s.Line(point1=(self.xrec[2]+L_finite*P, self.yrec[2]+H_finite*P), point2=(self.xrec[3]-L_finite*P, self.yrec[3]+H_finite*P))
        left = s.Line(point1=(self.xrec[3]-L_finite*P, self.yrec[3]+H_finite*P), point2=(self.xrec[0]-L_finite*P, self.yrec[0]-H_finite*P))
        p = self.model.Part(name='cut_temp', dimensionality=self.dimension, type=DEFORMABLE_BODY)
        self.model.parts['cut_temp'].BaseShell(sketch=s)
        s.unsetPrimaryObject()

        del self.model.sketches['cut_sketch']
        s = self.model.ConstrainedSketch(name='partition', sheetSize=50)
        UC_assem_list = []

        for i,indx in enumerate(self.UC_grid):
            pos = [indx[1]+origin_base[0],indx[2]+origin_base[1],indx[3]]
            UC = self.UC_list[indx[0]]
            UC_assem = self.duplicate_UC(UC,pos,i)
            UC_assem_list.append(UC_assem)

        tmerge = np.zeros(len(self.model.rootAssembly.instances.keys()),dtype=dict)

        for j in range(len(self.model.rootAssembly.instances.keys())):
            tmerge[j] = self.model.rootAssembly.instances[self.model.rootAssembly.instances.keys()[j]]

        self.model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
            instances=list(tmerge), 
            name='merge_temp', originalInstances=DELETE)
        self.model.rootAssembly.Instance(dependent=ON, name='cut_temp-1', 
            part=self.model.parts['cut_temp'])
        self.model.rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
            self.model.rootAssembly.instances['cut_temp-1'], ), 
            instanceToBeCut=self.model.rootAssembly.instances['merge_temp-1'], 
            name=self.name_assem, originalInstances=DELETE)

        del self.model.parts['cut_temp']
        del self.model.parts['merge_temp']
        self.TM_remove_faces(self.model.parts[self.name_assem])

    def TM_remove_faces(self,part):
        tol = 0.0001
        for f in part.faces:
            if f.getSize(printResults=False) < tol:
                edges = f.getEdges()
                edgeTuple = ()
                for e in edges:
                    edgeTuple = edgeTuple+(part.edges[e],)
                try:
                    part.RemoveRedundantEntities(edgeList=edgeTuple)
                except:
                    print("No redundent entries found.")

    def duplicate_UC(self,UC,pos,i):
        UC_assem = self.model.rootAssembly.Instance(dependent=ON, name='UC-'+str(i), part=
             self.model.parts[UC.name_part])
        self.model.rootAssembly.translate(instanceList=('UC-'+str(i), ), vector=(
            pos[0], pos[1], 0.0))
        return UC_assem


    def isBetween(self,a, b, c):
        if max(a[0],b[0]) > c[0] and min(a[0],b[0]) < c[0] and max(a[1],b[1]) > c[1] and min(a[1],b[1]) < c[1]:
            return True
        else:
            return False
        
    def print_finite(self,name_print,name_file,instron_tabs=True,hook_tabs=False):

        self.model.Part(name=name_print, objectToCopy=
            self.model.parts[self.name_assem])

        p = self.model.parts[name_print]

        clampW = self.width
        clampH = 40

        #counter clock-wise points of the base fabric
        xrec = [0, self.width, self.width, 0]
        yrec = [0, 0, -2.5, -2.5]

        N_top_hooks = 16.0
        N_bot_hooks = 12.0
        N_hooks = [N_bot_hooks,N_top_hooks]
        radius_hook = 1.1/2
        for i,Nh in enumerate(N_hooks):
            # print('Width = {}'.format(self.width))
            d_hook = self.width/(Nh)
            # print('D_hook = {}'.format(d_hook))
            s = self.model.ConstrainedSketch(name='hook_clampSketch_'+str(i), sheetSize=50)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=STANDALONE)
            s.Line(point1=(xrec[0], yrec[0]), point2=(xrec[1], yrec[1]))
            s.Line(point1=(xrec[1], yrec[1]), point2=(xrec[2], yrec[2]))
            s.Line(point1=(xrec[2], yrec[2]), point2=(xrec[3], yrec[3]))
            s.Line(point1=(xrec[3], yrec[3]), point2=(xrec[0], yrec[0]))


            s.CircleByCenterPerimeter(center=(
                d_hook/2-1.5, yrec[3]/2), point1=(d_hook/2-1.5, yrec[3]/2+radius_hook))
            s.CircleByCenterPerimeter(center=(
                d_hook/2+1.5, yrec[3]/2), point1=(d_hook/2+1.5, yrec[3]/2+radius_hook))
            s.linearPattern(angle1=0.0, angle2=
                90.0, geomList=(s.geometry.findAt(
                (d_hook/2+1.5, yrec[3]/2+radius_hook), ), 
                s.geometry.findAt((d_hook/2-1.5, yrec[3]/2+radius_hook), 
                )), number1=int(Nh), number2=1, spacing1=d_hook, 
                spacing2=0, vertexList=())
            p = self.model.Part(name='hook_clamp_'+str(i), dimensionality=self.dimension, type=DEFORMABLE_BODY)
            
            self.model.parts['hook_clamp_'+str(i)].BaseShell(sketch=s)
            s.unsetPrimaryObject()
            del s

        #counter clock-wise points of the base fabric
        xrec = [0, clampW, clampW, 0]
        yrec = [0, 0, -clampH, -clampH]

        ################################################
        # Creating base fabric
        ################################################

        s = self.model.ConstrainedSketch(name='clampSketch', sheetSize=50)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=STANDALONE)
        s.Line(point1=(xrec[0], yrec[0]), point2=(xrec[1], yrec[1]))
        s.Line(point1=(xrec[1], yrec[1]), point2=(xrec[2], yrec[2]))
        s.Line(point1=(xrec[2], yrec[2]), point2=(xrec[3], yrec[3]))
        s.Line(point1=(xrec[3], yrec[3]), point2=(xrec[0], yrec[0]))

        p = self.model.Part(name='clamp', dimensionality=self.dimension, type=DEFORMABLE_BODY)
        self.model.parts['clamp'].BaseShell(sketch=s)
        s.unsetPrimaryObject()

        s1 = self.model.ConstrainedSketch(name='tabSketch', sheetSize=50)
        tab_size = 10
        s1.rectangle(point1=(0, 0), 
            point2=(tab_size, tab_size))

        p = self.model.Part(name='tab', dimensionality=self.dimension, type=DEFORMABLE_BODY)
        self.model.parts['tab'].BaseShell(sketch=s1)
        s1.unsetPrimaryObject()
    
        del self.model.sketches['clampSketch']
        del self.model.sketches['tabSketch']

        ### Prepare Unit Cell ###
        if instron_tabs:
            self.model.rootAssembly.Instance(dependent=ON, name='clampBottom', part=self.model.parts['clamp'])
            self.model.rootAssembly.Instance(dependent=ON, name='clampTop', part=self.model.parts['clamp'])
            self.model.rootAssembly.Instance(dependent=ON, name='tabBL', part=self.model.parts['tab'])
            self.model.rootAssembly.Instance(dependent=ON, name='tabBR', part=self.model.parts['tab'])
            self.model.rootAssembly.Instance(dependent=ON, name='tabTL', part=self.model.parts['tab'])
            self.model.rootAssembly.Instance(dependent=ON, name='tabTR', part=self.model.parts['tab'])
            self.model.rootAssembly.Instance(dependent=ON, name='printTemp', part=self.model.parts[name_print])
            # Translate Tile
            self.model.rootAssembly.translate(instanceList=('clampTop', ), 
            vector=(0,self.height+abs(yrec[3]),0))
            self.model.rootAssembly.translate(instanceList=('tabBL', ), 
            vector=(-tab_size,-clampH/2-tab_size,0))
            self.model.rootAssembly.translate(instanceList=('tabBR', ), 
            vector=(clampW,-clampH/2-tab_size,0))
            self.model.rootAssembly.translate(instanceList=('tabTL', ), 
            vector=(-tab_size,self.height+clampH/2,0))
            self.model.rootAssembly.translate(instanceList=('tabTR', ), 
            vector=(clampW,self.height+clampH/2,0))

            self.model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
                instances=(self.model.rootAssembly.instances['clampBottom'], 
                self.model.rootAssembly.instances['clampTop'], 
                self.model.rootAssembly.instances['printTemp'],
                self.model.rootAssembly.instances['tabBL'],
                self.model.rootAssembly.instances['tabBR'],
                self.model.rootAssembly.instances['tabTL'],
                self.model.rootAssembly.instances['tabTR']), name=
                'printSVG_sample', originalInstances=DELETE)

            
            del self.model.parts['clamp']
            del self.model.parts['tab']
            del self.model.parts[name_print]
            del self.model.rootAssembly.instances['printSVG_sample-1']
            p = self.model.parts['printSVG_sample']
        elif hook_tabs:
            self.model.rootAssembly.Instance(dependent=ON, name='clampBottom', part=self.model.parts['hook_clamp_0'])
            self.model.rootAssembly.Instance(dependent=ON, name='clampTop', part=self.model.parts['hook_clamp_1'])
            self.model.rootAssembly.Instance(dependent=ON, name='printTemp', part=self.model.parts[name_print])
            # Translate Tile
            self.model.rootAssembly.translate(instanceList=('clampTop', ), 
            vector=(0,self.height+2.5,0))

            self.model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
                instances=(self.model.rootAssembly.instances['clampBottom'], 
                self.model.rootAssembly.instances['printTemp'],
                self.model.rootAssembly.instances['clampTop']), name=
                'printSVG_sample', originalInstances=DELETE)
            
            del self.model.parts['clamp']
            del self.model.parts['tab']
            del self.model.parts[name_print]
            del self.model.rootAssembly.instances['printSVG_sample-1']
            p = self.model.parts['printSVG_sample']
            
        else:
            p = self.model.parts[name_print]

        session.viewports['Viewport: 1'].setValues(displayedObject=p)
        session.viewports['Viewport: 1'].partDisplay.setValues(renderStyle=WIREFRAME)
        session.printToFile(fileName=name_file, format=SVG, canvasObjects=(
            session.viewports['Viewport: 1'], ))
        print('SVG file '+'printSVG_sample-1'+' has been saved.')