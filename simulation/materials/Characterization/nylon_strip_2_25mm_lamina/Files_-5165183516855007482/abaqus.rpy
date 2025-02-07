# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2021 replay file
# Internal Version: 2020_03_06-09.50.37 167380
# Run by david on Tue Apr 16 15:44:58 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=58.6115989685059, 
    height=90.4732894897461)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='job-strip-angle-0-75.odb')
#: Model: C:/Users/david/Documents/GitHub/R001-Textile-Metamaterials/materials/Characterization/nylon_strip_2_25mm_lamina/Files_-5165183516855007482/job-strip-angle-0-75.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       4
#: Number of Node Sets:          5
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=219.238, 
    farPlane=340.013, width=156.861, height=97.7426, cameraPosition=(141.833, 
    185.228, 145.026), cameraUpVector=(-0.778526, -0.331652, 0.532827), 
    cameraTarget=(-0.987626, 42.4069, 2.20567))
session.viewports['Viewport: 1'].view.setValues(nearPlane=226.763, 
    farPlane=296.021, width=162.245, height=101.098, cameraPosition=(251.034, 
    93.1999, -5.25585), cameraUpVector=(-0.348162, -0.134475, 0.927739), 
    cameraTarget=(11.6079, 31.7922, -15.1282))
session.viewports['Viewport: 1'].view.setValues(nearPlane=223.599, 
    farPlane=305.259, width=159.982, height=99.6875, cameraPosition=(246.308, 
    69.303, 55.0783), cameraUpVector=(-0.601923, 0.174913, 0.779162), 
    cameraTarget=(11.3544, 30.5105, -11.8922))
session.viewports['Viewport: 1'].view.setValues(nearPlane=221.029, 
    farPlane=317.371, width=158.143, height=98.5416, cameraPosition=(197.365, 
    183.848, 51.9984), cameraUpVector=(-0.604297, -0.127125, 0.786552), 
    cameraTarget=(8.19741, 37.8991, -12.0909))
