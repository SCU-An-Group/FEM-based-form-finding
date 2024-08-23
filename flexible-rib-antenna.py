# -*- coding: mbcs -*-
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
import numpy as np
import csv
import os
from os import path
import regionToolset
import math

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# execfile("parameter_for_spatial.py")
DampingFactor = 1e-8


theta = 0. * pi/180.0
m = 16#肋个数
angle = 180.0/m
F0 = 225.0#焦距
D = 500.0#口径
d = 50.0#中心轮毂直径
n = 40#生成抛物线点个数，越大抛物面精度越高
x = np.linspace((d/2.0*cos(theta)), sqrt(D**2/4.0 - d**2/4.0*sin(theta)**2), n)
y = 1./(4.*F0) * x**2
z = 1./(4.*F0) * (x**2+y**2)
x1 = x*cos(angle * pi/180.0)
y1 = 1./(4.*F0) *x**2
z1 = x*sin(angle * pi/180.0)
x2 = x*cos(-angle * pi/180.0)
y2 = 1./(4.*F0) *x**2
z2 = x*sin(-angle * pi/180.0)




# Parameters

#cable
density_c = 5.4e-9
Ec = 40000.0;
nu_c = 0.35
alpha_c = -3.0e-06;
area_c = 1
pretension_c_out = 1;
pretension_c_in = 0.1
delta_T_c_out = pretension_c_out/area_c/Ec/alpha_c;
delta_T_c_in = pretension_c_in/area_c/Ec/alpha_c;

#membrane
density_m = 2.16e-9
Em = 1.7;
nu_m = 0.3;
thickness_m = 0.08
alpha_m = 1e-05;
pretension_m = 0.15;
delta_T_m = pretension_m * (1-nu_m)/Em/alpha_m

#rib
density_r = 7.86e-9
Er = 206000
nu_r = 0.3
thickness_r = 0.1
theta_rib = 30
R_rib = 15



LargeNumber = 999999999999.0


Mdb()
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Membrane', type=
    DEFORMABLE_BODY)

#create point
name_point1 = []
name_point2 = []
for i in range(0,n):
    name_point1.append('P1_' + str(i+1))
    name_point2.append('P2_' + str(i+1))


for i in range(len(x)):
    name_point1[i] = mdb.models['Model-1'].parts['Membrane'].DatumPointByCoordinate(coords=(x1[i], y1[i], z1[i]))
    name_point2[i] = mdb.models['Model-1'].parts['Membrane'].DatumPointByCoordinate(coords=(x2[i], y2[i], z2[i]))

#create line
name_line1 = []
name_line2 = []

for i in range(0,n-1):
    name_line1.append('L1_' + str(i+1))
    name_line2.append('L2_' + str(i+1))
    
name_line3 = []
for i in range(0,n):
    name_line3.append('L3_' + str(i+1))


for i in range(len(x)-1):
    name_line1[i] = mdb.models['Model-1'].parts['Membrane'].WirePolyLine(mergeType=IMPRINT, meshable=
        ON, points=((mdb.models['Model-1'].parts['Membrane'].datums[name_point1[i].id], 
            mdb.models['Model-1'].parts['Membrane'].datums[name_point1[i+1].id])))
    name_line2[i] = mdb.models['Model-1'].parts['Membrane'].WirePolyLine(mergeType=IMPRINT, meshable=
        ON, points=((mdb.models['Model-1'].parts['Membrane'].datums[name_point2[i].id], 
            mdb.models['Model-1'].parts['Membrane'].datums[name_point2[i+1].id])))

for i in range(len(x)):
    name_line3[i] = mdb.models['Model-1'].parts['Membrane'].WirePolyLine(mergeType=IMPRINT, meshable=
        ON, points=((mdb.models['Model-1'].parts['Membrane'].datums[name_point1[i].id], 
            mdb.models['Model-1'].parts['Membrane'].datums[name_point2[i].id])))

#create faces
for i in range(len(x)-1):
    setname = 'rectangle-' + str(i)
    edges=mdb.models['Model-1'].parts['Membrane'].edges.getByBoundingBox(xMin=x1[i], yMin=-LargeNumber, zMin=-LargeNumber,
                                                             xMax=x2[i+1], yMax=LargeNumber, zMax=LargeNumber) 
    mdb.models['Model-1'].parts['Membrane'].Set(edges=edges,name=setname)
    mdb.models['Model-1'].parts['Membrane'].CoverEdges(edgeList=edges,tryAnalytical=True)


#repalce faces
setname = 'mem_allfaces' 
faces=mdb.models['Model-1'].parts['Membrane'].faces.getByBoundingBox(xMin=x1[0],yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=x2[n-1],yMax=LargeNumber,zMax=LargeNumber) 
# mdb.models['Model-1'].parts['Membrane'].Set(faces=faces,name=setname)
# mdb.models['Model-1'].parts['Membrane'].ReplaceFaces(faceList=faces, stitch=True)



# create rib
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='rib', type=
    DEFORMABLE_BODY)


# create rib point
name_point1 = []

for i in range(0,n):
    name_point1.append('P1_' + str(i+1))



for i in range(len(x)):
    name_point1[i] = mdb.models['Model-1'].parts['rib'].DatumPointByCoordinate(coords=(x1[i], y1[i], z1[i]))


# create rib line
name_line1 = []


for i in range(0,n-1):
    name_line1.append('L1_' + str(i+1))

    

for i in range(len(x)-1):
    name_line1[i] = mdb.models['Model-1'].parts['rib'].WirePolyLine(mergeType=IMPRINT, meshable=
        ON, points=((mdb.models['Model-1'].parts['rib'].datums[name_point1[i].id], 
            mdb.models['Model-1'].parts['rib'].datums[name_point1[i+1].id])))

mdb.models['Model-1'].parts['rib'].DatumAxisByTwoPoint(point1=
    mdb.models['Model-1'].parts['rib'].datums[1], point2=
    mdb.models['Model-1'].parts['rib'].datums[n])
mdb.models['Model-1'].ConstrainedSketch(gridSpacing=11.76, name='__profile__', 
    sheetSize=470.53)
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(
    -235.265, 0.0), point2=(235.265, 0.0))
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    -235.265), point2=(0.0, 235.265))
mdb.models['Model-1'].parts['rib'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(R_rib*cos(theta_rib/2), 
    R_rib*sin(theta_rib/2)), direction=COUNTERCLOCKWISE, point1=(0.0, 0.0), point2=(0.0, 2*R_rib*sin(theta_rib/2)))
mdb.models['Model-1'].parts['rib'].ShellSweep(path=
    mdb.models['Model-1'].parts['rib'].edges, 
    profile=mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=
    RIGHT, sketchUpEdge=
    mdb.models['Model-1'].parts['rib'].datums[2*n])
del mdb.models['Model-1'].sketches['__profile__']



#creat cable
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='cable', type=
    DEFORMABLE_BODY)


name_point1 = ['P1_1','P1_2']
name_point2 = ['P2_1','P2_2']

name_point1[0] = mdb.models['Model-1'].parts['cable'].DatumPointByCoordinate(coords=(x1[0], y1[0], z1[0]))
name_point1[1] = mdb.models['Model-1'].parts['cable'].DatumPointByCoordinate(coords=(x1[n-1], y1[n-1], z1[n-1]))
name_point2[0] = mdb.models['Model-1'].parts['cable'].DatumPointByCoordinate(coords=(x2[0], y2[0], z2[0]))
name_point2[1] = mdb.models['Model-1'].parts['cable'].DatumPointByCoordinate(coords=(x2[n-1], y2[n-1], z2[n-1]))

name_line1 = ['L1_1','L1_2']
name_line1[0] = mdb.models['Model-1'].parts['cable'].WirePolyLine(mergeType=IMPRINT, meshable=
    ON, points=((mdb.models['Model-1'].parts['cable'].datums[name_point1[0].id], 
        mdb.models['Model-1'].parts['cable'].datums[name_point2[0].id])))
name_line1[1] = mdb.models['Model-1'].parts['cable'].WirePolyLine(mergeType=IMPRINT, meshable=
    ON, points=((mdb.models['Model-1'].parts['cable'].datums[name_point1[1].id], 
        mdb.models['Model-1'].parts['cable'].datums[name_point2[1].id])))




# create sets

setname = 'AllRib'
rib_edges=mdb.models['Model-1'].parts['rib'].faces.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) 
mdb.models['Model-1'].parts['rib'].Set(faces=rib_edges,name=setname)

setname = 'AllMem'
mem_faces=mdb.models['Model-1'].parts['Membrane'].faces.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) 
mdb.models['Model-1'].parts['Membrane'].Set(faces=mem_faces,name=setname)

setname = 'AllCable'
cable_edges=mdb.models['Model-1'].parts['cable'].edges.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) 
mdb.models['Model-1'].parts['cable'].Set(edges=cable_edges,name=setname)


# Materials
mdb.models['Model-1'].Material(name='Material-membrane')
mdb.models['Model-1'].materials['Material-membrane'].Elastic(table=((Em, 
    nu_m), ))
mdb.models['Model-1'].materials['Material-membrane'].Expansion(table=((alpha_m, 
    ), ))

mdb.models['Model-1'].Material(name='Material-cables')
mdb.models['Model-1'].materials['Material-cables'].Elastic(table=((Ec, 
    nu_c), ))
mdb.models['Model-1'].materials['Material-cables'].Expansion(table=((alpha_c, 
    ), ))
mdb.models['Model-1'].materials['Material-cables'].elastic.setValues(
    noCompression=ON)

mdb.models['Model-1'].Material(name='Material-ribs')
mdb.models['Model-1'].materials['Material-ribs'].Elastic(table=((Er, 
    nu_r), ))


#properties

#membrane
mdb.models['Model-1'].MembraneSection(material='Material-membrane', name=
    'Section-2', poissonDefinition=DEFAULT, thickness=thickness_m, thicknessField='', 
    thicknessType=UNIFORM)
mdb.models['Model-1'].parts['Membrane'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mem_faces), sectionName=
    'Section-2', thicknessAssignment=FROM_SECTION)


mdb.models['Model-1'].parts['Membrane'].DatumCsysByThreePoints(coordSysType=
    CYLINDRICAL, name='Datum csys-1', origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 
    0.0), point2=(0.0, 0.0, 1.0))


#cable
mdb.models['Model-1'].TrussSection(area=area_c, material='Material-cables', name=
    'Section-cables')

mdb.models['Model-1'].parts['cable'].Set(edges=
    mdb.models['Model-1'].parts['cable'].edges, 
    name='AllCable')
mdb.models['Model-1'].parts['cable'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['cable'].sets['AllCable'], sectionName=
    'Section-cables', thicknessAssignment=FROM_SECTION)


#rib_shell
mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    material='Material-ribs', name='Section-ribs', nodalThicknessField='', 
    poissonDefinition=DEFAULT, preIntegrate=ON, thickness=thickness_r, thicknessField=
    '', thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
mdb.models['Model-1'].parts['rib'].SectionAssignment(offset=0.0, offsetField=''
    , offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['rib'].sets['AllRib'], sectionName='Section-ribs', 
    thicknessAssignment=FROM_SECTION)

#assembly

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Membrane-1', 
    part=mdb.models['Model-1'].parts['Membrane'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='cable-1', part=
    mdb.models['Model-1'].parts['cable'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='rib-1', part=
    mdb.models['Model-1'].parts['rib'])
mdb.models['Model-1'].rootAssembly.RadialInstancePattern(axis=(0.0, 1.0, 0.0), 
    instanceList=('Membrane-1', 'cable-1', 'rib-1'), number=16, point=(0.0, 
    0.0, 0.0), totalAngle=360.0)

#step
mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None, 
    continueDampingFactors=False, initialInc=0.0001, name='Step-1', nlgeom=ON, 
    previous='Initial', stabilizationMagnitude=DampingFactor, stabilizationMethod=
    DAMPING_FACTOR)
mdb.models['Model-1'].steps['Step-1'].setValues(minInc=1e-16)
mdb.models['Model-1'].steps['Step-1'].setValues(maxNumInc=10000)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'COORD'))
mdb.models['Model-1'].steps['Step-1'].Restart(frequency=0, numberIntervals=1, 
    overlay=OFF, timeMarks=OFF)

# #mesh


mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
    minSizeFactor=0.1, regions=(
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'], ), size=5.0)
mdb.models['Model-1'].rootAssembly.setMeshControls(elemShape=QUAD, regions=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].faces, )

mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=M3D4R, elemLibrary=STANDARD), ElemType(elemCode=M3D3, 
    elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=Region(mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].faces), )

for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
        minSizeFactor=0.1, regions=(
        mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-'+ str(i)], ), size=5.0)
    mdb.models['Model-1'].rootAssembly.setMeshControls(elemShape=QUAD, regions=
        mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-'+ str(i)].faces, )
    mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
        elemCode=M3D4R, elemLibrary=STANDARD), ElemType(elemCode=M3D3, 
        elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=Region(mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-'+ str(i)].faces), )

Region(faces=mem_faces)


mdb.models['Model-1'].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=
    mdb.models['Model-1'].rootAssembly.instances['cable-1'].edges, number=11)
mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=T3D2, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].rootAssembly.instances['cable-1'].edges, 
    ))

for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=
        mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-'+ str(i)].edges, number=11)
    mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
        elemCode=T3D2, elemLibrary=STANDARD), ), regions=
        (mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-'+ str(i)].edges, ))


mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
    minSizeFactor=0.1, regions=(
    mdb.models['Model-1'].rootAssembly.instances['rib-1'], ), size=5.0)
mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=S4R, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].faces, 
    ))

for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
        minSizeFactor=0.1, regions=(
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)], ), size=5.0)
    mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
        elemCode=S4R, elemLibrary=STANDARD), ), regions=
        (mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].faces, ))


mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-6'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-7'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-8'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-9'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-10'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-11'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-12'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-13'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-14'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-15'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-16'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-6'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-7'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-8'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-9'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-10'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-11'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-12'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-13'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-14'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-15'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-16'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-6'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-7'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-8'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-9'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-10'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-11'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-12'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-13'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-14'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-15'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-16']))



mdb.models['Model-1'].rootAssembly.Set(name='Membraneallnodes-1', nodes=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )
for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='Membraneallnodes-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-'+ str(i)].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )

mdb.models['Model-1'].rootAssembly.Set(name='riballnodes-1', nodes=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )
for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='riballnodes-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )

mdb.models['Model-1'].rootAssembly.Set(name='cableallnodes-1', nodes=
    mdb.models['Model-1'].rootAssembly.instances['cable-1'].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )
for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='cableallnodes-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-'+ str(i)].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )



mdb.models['Model-1'].rootAssembly.Set(name='inner_cable-1', nodes=
    mdb.models['Model-1'].rootAssembly.instances['cable-1'].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=y[0],zMax=LargeNumber) )
for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='inner_cable-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-'+ str(i)].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=y[0],zMax=LargeNumber) )

mdb.models['Model-1'].rootAssembly.Set(name='outer_cable-1', nodes=
    mdb.models['Model-1'].rootAssembly.instances['cable-1'].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )
for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='outer_cable-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-'+ str(i)].nodes.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )

mdb.models['Model-1'].rootAssembly.Set(name='fixed point of rib-1', nodes=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[0:2]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[82:85])

for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='fixed point of rib-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].nodes[0:2]+\
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].nodes[82:85])



mdb.models['Model-1'].rootAssembly.Surface(name='Surf-rib', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-2'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-3'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-4'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-5'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-6'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-7'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-8'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-9'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-10'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-11'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-12'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-13'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-14'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-15'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-16'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-membrane', 
    side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-2'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-3'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-4'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-5'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-6'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-7'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-8'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-9'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-10'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-11'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-12'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-13'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-14'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-15'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-16'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber))






mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.sets['riballnodes-1'], 
    name='Tie1-1', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-1']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Tie1-1'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.sets['riballnodes-2'], 
    name='Tie1-2', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-1']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Tie1-2'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)

for i in range(2,m):
    mdb.models['Model-1'].Tie(adjust=ON, master=
        mdb.models['Model-1'].rootAssembly.sets['riballnodes-'+ str(i)], 
        name='Tie'+ str(i) + '-1', positionToleranceMethod=COMPUTED, slave=
        mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-'+ str(i)]
        , thickness=ON, tieRotations=ON)
    mdb.models['Model-1'].constraints['Tie'+ str(i) + '-1'].setValues(positionTolerance=
        0.001, positionToleranceMethod=SPECIFIED)

for i in range(2,m):
    mdb.models['Model-1'].Tie(adjust=ON, master=
        mdb.models['Model-1'].rootAssembly.sets['riballnodes-'+ str(i+1)], 
        name='Tie'+ str(i) + '-2', positionToleranceMethod=COMPUTED, slave=
        mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-'+ str(i)]
        , thickness=ON, tieRotations=ON)
    mdb.models['Model-1'].constraints['Tie'+ str(i) + '-2'].setValues(positionTolerance=
        0.001, positionToleranceMethod=SPECIFIED)    

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.sets['riballnodes-16'], 
    name='Tie16-1', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-16']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Tie16-1'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.sets['riballnodes-1'], 
    name='Tie16-2', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-16']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Tie16-2'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)



mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-1'], 
    name='Tie17-1', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.sets['cableallnodes-1']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Tie17-1'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)

for i in range(2,m+1):
    mdb.models['Model-1'].Tie(adjust=ON, master=
        mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-'+ str(i)], 
        name='Tie17-'+ str(i), positionToleranceMethod=COMPUTED, slave=
        mdb.models['Model-1'].rootAssembly.sets['cableallnodes-'+ str(i)]
        , thickness=ON, tieRotations=ON)
    mdb.models['Model-1'].constraints['Tie17-'+ str(i)].setValues(positionTolerance=
        0.001, positionToleranceMethod=SPECIFIED)

#load
mdb.models['Model-1'].rootAssembly.regenerate()

for i in range(1,m+1):
    mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
        name='BC-'+ str(i), region=
        mdb.models['Model-1'].rootAssembly.sets['fixed point of rib-'+ str(i)])



for i in range(1,m+1):
    mdb.models['Model-1'].Temperature(createStepName='Initial', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(0.0, ), name='Membrane Predefined Field-'+ str(i), region=
        mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-'+ str(i)])
    mdb.models['Model-1'].predefinedFields['Membrane Predefined Field-'+ str(i)].setValuesInStep(
        magnitudes=(-delta_T_m, ), stepName='Step-1')

for i in range(1,m+1):
    mdb.models['Model-1'].Temperature(createStepName='Initial', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(0.0, ), name='innercable Predefined Field-'+ str(i), region=
        mdb.models['Model-1'].rootAssembly.sets['inner_cable-'+ str(i)])
    mdb.models['Model-1'].predefinedFields['innercable Predefined Field-'+ str(i)].setValuesInStep(
        magnitudes=(-delta_T_c_in, ), stepName='Step-1')

for i in range(1,m+1):
    mdb.models['Model-1'].Temperature(createStepName='Initial', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(0.0, ), name='outercable Predefined Field-'+ str(i), region=
        mdb.models['Model-1'].rootAssembly.sets['outer_cable-'+ str(i)])
    mdb.models['Model-1'].predefinedFields['outercable Predefined Field-'+ str(i)].setValuesInStep(
        magnitudes=(-delta_T_c_out, ), stepName='Step-1')


# job 
jobName = 'Iteration-1' 
mdb.saveAs(pathName=jobName)
mdb.Job(model='Model-1', name=jobName)

#### Check jobName only if it does not exist in the folder

if path.exists(jobName + '.odb') == True:
    print(jobName + " already exists")
else:
    mdb.jobs[jobName].writeInput(consistencyChecking=OFF) # write an inp file
    os.system('abaqus job=' + jobName + ' int')
    # mdb.jobs[jobName].submit()
    # mdb.jobs[jobName].waitForCompletion()

#### Delete .lck file when it exists
if path.exists(jobName + '.lck') == True:
    os.remove(jobName + '.lck')
else:
    print('...')


#---------------------------------------------------------------
# Iteration starts here
#---------------------------------------------------------------

n = 100 # iteration times
results=np.zeros((n,2))
for ind in range(1, n):

    itrModelName = 'Iteration-' + str(ind+1)
    
    mdb.Model(name=itrModelName, objectToCopy=mdb.models['Model-1'])
    
    ODB = openOdb(path = 'Iteration-' + str(ind) + '.odb')


    firstFrame = ODB.steps['Step-1'].frames[-1]
    coordinates = firstFrame.fieldOutputs['COORD']


     # Set new coordinates for membrane
    MembraneNodes=ODB.rootAssembly.nodeSets['MEMBRANEALLNODES-1']
    MembraneCoordField = coordinates.getSubset(region=MembraneNodes)
    MembraneCoord = MembraneCoordField.values
    NewCoord_Membrane = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['Membrane-1'].nodes), 3))
    
    for i in range(len(mdb.models[itrModelName].rootAssembly.instances['Membrane-1'].nodes)):
        NewCoord_Membrane[i][0]=np.float64(MembraneCoord[i].data[0]) 
        NewCoord_Membrane[i][1]=np.float64(MembraneCoord[i].data[1]) 
        NewCoord_Membrane[i][2]=np.float64(MembraneCoord[i].data[2])
    
#SET THE NEW COORDINATES
    mdb.models[itrModelName].rootAssembly.editNode(
        nodes=mdb.models[itrModelName].rootAssembly.instances['Membrane-1'].nodes,
        coordinates=NewCoord_Membrane)

    for j in range(2,m+1):
        MembraneNodes=ODB.rootAssembly.nodeSets['MEMBRANEALLNODES-'+ str(j)]
        MembraneCoordField = coordinates.getSubset(region=MembraneNodes)
        MembraneCoord = MembraneCoordField.values
        NewCoord_Membrane = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['Membrane-1-rad-'+ str(j)].nodes), 3))
        
        for i in range(len(mdb.models[itrModelName].rootAssembly.instances['Membrane-1-rad-'+ str(j)].nodes)):
            NewCoord_Membrane[i][0]=np.float64(MembraneCoord[i].data[0]) 
            NewCoord_Membrane[i][1]=np.float64(MembraneCoord[i].data[1]) 
            NewCoord_Membrane[i][2]=np.float64(MembraneCoord[i].data[2])
        
        mdb.models[itrModelName].rootAssembly.editNode(
            nodes=mdb.models[itrModelName].rootAssembly.instances['Membrane-1-rad-'+ str(j)].nodes,
            coordinates=NewCoord_Membrane)


        mdb.models[itrModelName].rootAssembly.regenerate()



    # Set new coordinates for rib
    RibNodes=ODB.rootAssembly.nodeSets['RIBALLNODES-1']
    RibCoordField = coordinates.getSubset(region=RibNodes)
    RibCoord = RibCoordField.values
    NewCoord_Rib = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['rib-1'].nodes), 3))
    
    for i in range(len(mdb.models[itrModelName].rootAssembly.instances['rib-1'].nodes)):
        NewCoord_Rib[i][0]=np.float64(RibCoord[i].data[0])
        NewCoord_Rib[i][1]=np.float64(RibCoord[i].data[1])
        NewCoord_Rib[i][2]=np.float64(RibCoord[i].data[2])
    
    #SET THE NEW COORDINATES
    mdb.models[itrModelName].rootAssembly.editNode(
        nodes=mdb.models[itrModelName].rootAssembly.instances['rib-1'].nodes,
        coordinates=NewCoord_Rib)
    
    mdb.models[itrModelName].rootAssembly.regenerate()

    for j in range(2,m+1):
        RibNodes=ODB.rootAssembly.nodeSets['RIBALLNODES-'+ str(j)]
        RibCoordField = coordinates.getSubset(region=RibNodes)
        RibCoord = RibCoordField.values
        NewCoord_Rib = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-'+ str(j)].nodes), 3))
        
        for i in range(len(mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-'+ str(j)].nodes)):
            NewCoord_Rib[i][0] = np.float64(RibCoord[i].data[0])
            NewCoord_Rib[i][1] = np.float64(RibCoord[i].data[1])
            NewCoord_Rib[i][2] = np.float64(RibCoord[i].data[2])
        
        mdb.models[itrModelName].rootAssembly.editNode(
            nodes=mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-'+ str(j)].nodes,
            coordinates=NewCoord_Rib)
    
    mdb.models[itrModelName].rootAssembly.regenerate()


    
    mdb.models[itrModelName].InitialState(createStepName='Initial', endIncrement=
    STEP_END, endStep=LAST_STEP, fileName='Iteration-' + str(ind), instances=(
    mdb.models[itrModelName].rootAssembly.instances['rib-1'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-16'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-15'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-14'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-13'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-12'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-11'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-10'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-9'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-8'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-7'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-6'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-5'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-4'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-3'], 
    mdb.models[itrModelName].rootAssembly.instances['rib-1-rad-2']), name=
    'Predefined Field-rib', updateReferenceConfiguration=ON)


    mdb.models[itrModelName].ContactProperty('IntProp-1')
    mdb.models[itrModelName].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=ROUGH)
    mdb.models[itrModelName].interactionProperties['IntProp-1'].NormalBehavior(
        allowSeparation=OFF, constraintEnforcementMethod=DEFAULT, 
        pressureOverclosure=HARD)
    mdb.models[itrModelName].SurfaceToSurfaceContactStd(adjustMethod=NONE, 
        clearanceRegion=None, createStepName='Step-1', datumAxis=None, 
        initialClearance=OMIT, interactionProperty='IntProp-1', master=
        mdb.models[itrModelName].rootAssembly.surfaces['Surf-rib'], name='Int-1', 
        slave=mdb.models[itrModelName].rootAssembly.surfaces['Surf-membrane'], 
        sliding=FINITE, thickness=ON)

# Set new coordinates for cable
    CableNodes=ODB.rootAssembly.nodeSets['CABLEALLNODES-1']
    CableCoordField = coordinates.getSubset(region=CableNodes)
    CableCoord = CableCoordField.values
    NewCoord_Cable = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['cable-1'].nodes), 3))
    
    for i in range(len(mdb.models[itrModelName].rootAssembly.instances['cable-1'].nodes)):
        NewCoord_Cable[i][0]=np.float64(CableCoord[i].data[0])
        NewCoord_Cable[i][1]=np.float64(CableCoord[i].data[1])
        NewCoord_Cable[i][2]=np.float64(CableCoord[i].data[2])
    
#SET THE NEW COORDINATES
    mdb.models[itrModelName].rootAssembly.editNode(
        nodes=mdb.models[itrModelName].rootAssembly.instances['cable-1'].nodes,
        coordinates=NewCoord_Cable)
      

    for j in range(2,m+1):
        CableNodes=ODB.rootAssembly.nodeSets['CABLEALLNODES-'+ str(j)]
        CableCoordField = coordinates.getSubset(region=CableNodes)
        CableCoord = CableCoordField.values
        NewCoord_Cable = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['cable-1-rad-'+ str(j)].nodes), 3))
        
        for i in range(len(mdb.models[itrModelName].rootAssembly.instances['cable-1-rad-'+ str(j)].nodes)):
            NewCoord_Cable[i][0]=np.float64(CableCoord[i].data[0])
            NewCoord_Cable[i][1]=np.float64(CableCoord[i].data[1])
            NewCoord_Cable[i][2]=np.float64(CableCoord[i].data[2])
        
        mdb.models[itrModelName].rootAssembly.editNode(
            nodes=mdb.models[itrModelName].rootAssembly.instances['cable-1-rad-'+ str(j)].nodes,
            coordinates=NewCoord_Cable)


    mdb.models[itrModelName].steps['Step-1'].setValues(initialInc=1e-6)


    mdb.models[itrModelName].rootAssembly.regenerate()
    

    a = mdb.models[itrModelName].rootAssembly

    x = []
    y = []
    z = []
    Area = []
    normal=[]

    for i in range(len(a.instances['Membrane-1'].elements)):

        elements = a.instances['Membrane-1'].elements[i:i+1]
        element = a.instances['Membrane-1']
        region = regionToolset.Region(elements=elements)
        properties = a.getMassProperties(regions=region)
        area = properties['area']
        areaCentroid = properties['areaCentroid']

        x.append(areaCentroid[0])
        y.append(areaCentroid[2])
        z.append(areaCentroid[1])
        Area.append(area)
        normal.append(atan(-element.elementFaces[i].getNormal()[0]/element.elementFaces[i].getNormal()[1]))


    prediction_z = []
    Area_projected = []

    for i in range(len(z)):
        prediction_z.append((1/900.0)*(x[i]**2 + y[i]**2))
        Area_projected.append(Area[i]*cos(normal[i]))


    error = []
    for i in range(len(z)):
        error.append(z[i] - prediction_z[i])
     
    totalArea = sum(Area_projected)

    squaredError = []
    absError = []

    for i in range(len(z)):
        squaredError.append(error[i]*error[i]*Area_projected[i])

    for val in error:
        # squaredError.append(val * val)#target-prediction之差平方 
        absError.append(abs(val))#误差绝对值

    targetDeviation = []
    targetMean = sum(z) / len(z)#target平均值
    for val in z:
        targetDeviation.append((val - targetMean) * (val - targetMean))

    results[ind-1]=(ind,sqrt(sum(squaredError) / totalArea))


    ODB.close()


    # job 
    jobName = 'Iteration-' + str(ind+1) 
    mdb.Job(model=itrModelName, name=jobName)
    # mdb.saveAs(pathName=jobName)#no cae
    # mdb.jobs[jobName].setValues(nodalOutputPrecision=FULL)


    if path.exists(jobName + '.odb') == True:
        print(jobName + " already exists")
    else:
        mdb.jobs[jobName].writeInput(consistencyChecking=OFF) # write an inp file
        os.system('abaqus job=' + jobName + ' int')
        # mdb.jobs[jobName].submit()
        # mdb.jobs[jobName].waitForCompletion()
    
    #### Delete .lck file when it exists
    if path.exists(jobName + '.lck') == True:
        os.remove(jobName + '.lck')
    else:
        print("...")


    ind = ind + 1    

csvFile = open("rmse-new-projected.csv","wb")
writer = csv.writer(csvFile)
writer.writerow(["x","y"])

for value in results:
    writer.writerow(value)
csvFile.close()
