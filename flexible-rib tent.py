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

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

DampingFactor = 1e-8

m = 6
theta = 0. * pi/180.0
F0 = 450.0
D = 2000.0
d = 0.0
h = d**2/16.0/F0*(1+sin(theta)**2)
n = 100
x = np.linspace((d/2.0*cos(theta)), sqrt(D**2/4.0 - d**2/4.0*sin(theta)**2), n)
y = 1./(4.*F0) * x**2
z = 1./(4.*F0) * (x**2+y**2)
x1 = x*cos(30. * pi/180.0)
y1 = x
z1 = x*sin(30. * pi/180.0)
x2 = x*cos(-30. * pi/180.0)
y2 = x
z2 = x*sin(-30. * pi/180.0)




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
pretension_m = 0.1;
delta_T_m = pretension_m * (1-nu_m)/Em/alpha_m

#rib
density_r = 7.86e-9
Er = 2000
nu_r = 0.3
thickness_r = 2



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

for i in range(len(x)-1):
    name_line3[i] = mdb.models['Model-1'].parts['Membrane'].WirePolyLine(mergeType=IMPRINT, meshable=
        ON, points=((mdb.models['Model-1'].parts['Membrane'].datums[name_point1[i+1].id], 
            mdb.models['Model-1'].parts['Membrane'].datums[name_point2[i+1].id])))

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

name_point1 = []
name_point2 = ['P2_1']

for i in range(0,n):
    name_point1.append('P1_' + str(i+1))


name_point2[0] = mdb.models['Model-1'].parts['rib'].DatumPointByCoordinate(coords=(x2[2], y2[2], z2[2]))
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

mdb.models['Model-1'].parts['rib'].DatumPlaneByThreePoints(isDependent=False, 
    point1=mdb.models['Model-1'].parts['rib'].datums[2], point2=
    mdb.models['Model-1'].parts['rib'].datums[4], point3=
    mdb.models['Model-1'].parts['rib'].datums[1])
mdb.models['Model-1'].parts['rib'].DatumAxisByNormalToPlane(isDependent=False, 
    plane=mdb.models['Model-1'].parts['rib'].datums[201], point=
    mdb.models['Model-1'].parts['rib'].datums[2])
mdb.models['Model-1'].ConstrainedSketch(gridSpacing=111.9, name='__profile__', 
    sheetSize=4476.16)
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(
    -2238.08, 0.0), point2=(2238.08, 0.0))
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    -2238.08), point2=(0.0, 2238.08))
mdb.models['Model-1'].parts['rib'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
    -10*sin(22.2073425 * pi/180.0), 10*cos(22.2073425 * pi/180.0)))
mdb.models['Model-1'].parts['rib'].ShellSweep(flipSweepDirection=ON, path=
    mdb.models['Model-1'].parts['rib'].edges
    , profile=mdb.models['Model-1'].sketches['__profile__'], profileNormal=ON, 
    sketchOrientation=RIGHT, sketchUpEdge=
    mdb.models['Model-1'].parts['rib'].datums[202])
del mdb.models['Model-1'].sketches['__profile__']



#creat cable
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='cable', type=
    DEFORMABLE_BODY)


name_point1 = ['P1_1']
name_point2 = ['P2_1']

name_point1[0] = mdb.models['Model-1'].parts['cable'].DatumPointByCoordinate(coords=(x1[n-1], y1[n-1], z1[n-1]))

name_point2[0] = mdb.models['Model-1'].parts['cable'].DatumPointByCoordinate(coords=(x2[n-1], y2[n-1], z2[n-1]))

name_line1 = ['L1_1']
name_line1[0] = mdb.models['Model-1'].parts['cable'].WirePolyLine(mergeType=IMPRINT, meshable=
    ON, points=((mdb.models['Model-1'].parts['cable'].datums[name_point1[0].id], 
        mdb.models['Model-1'].parts['cable'].datums[name_point2[0].id])))




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
    CYLINDRICAL, name='Datum csys-1', origin=(0.0, 0.0, 0.0), point1=
    mdb.models['Model-1'].parts['Membrane'].InterestingPoint(
    mdb.models['Model-1'].parts['Membrane'].edges.findAt((866.025404, 1000.0, 
    250.0), ), MIDDLE), point2=
    mdb.models['Model-1'].parts['Membrane'].datums[2*n])

mdb.models['Model-1'].parts['Membrane'].MaterialOrientation(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, fieldName='', localCsys=
    mdb.models['Model-1'].parts['Membrane'].datums[698], orientationType=SYSTEM
    , region=Region(faces=mdb.models['Model-1'].parts['Membrane'].faces))

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
    instanceList=('Membrane-1', 'cable-1', 'rib-1'), number=6, point=(0.0, 
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
mdb.models['Model-1'].rootAssembly.makeDependent(instances=(
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-6'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-6'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-6']))


mdb.models['Model-1'].parts['Membrane'].setMeshControls(elemShape=QUAD, regions=
    mdb.models['Model-1'].parts['Membrane'].faces)
mdb.models['Model-1'].parts['Membrane'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=50.0)
mdb.models['Model-1'].parts['Membrane'].setElementType(elemTypes=(ElemType(
    elemCode=M3D4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF), ElemType(
    elemCode=M3D3, elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=(
    mdb.models['Model-1'].parts['Membrane'].faces, ))
mdb.models['Model-1'].parts['Membrane'].generateMesh()

mdb.models['Model-1'].parts['rib'].seedEdgeByNumber(constraint=FINER, edges=mdb.models['Model-1'].parts['rib'].edges, number=3)
mdb.models['Model-1'].parts['rib'].generateMesh()
mdb.models['Model-1'].parts['rib'].setElementType(elemTypes=(ElemType(
    elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    hourglassControl=ENHANCED), ElemType(elemCode=S3, elemLibrary=STANDARD)), 
    regions=(mdb.models['Model-1'].parts['rib'].faces, ))

mdb.models['Model-1'].parts['cable'].seedEdgeByNumber(constraint=FINER, edges=
    mdb.models['Model-1'].parts['cable'].edges, number=20)
mdb.models['Model-1'].parts['cable'].generateMesh()
mdb.models['Model-1'].parts['cable'].setElementType(elemTypes=(ElemType(
    elemCode=T3D2, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].parts['cable'].edges, ))



mdb.models['Model-1'].rootAssembly.makeIndependent(instances=(
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-6'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-6'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-2'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-3'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-4'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-5'], 
    mdb.models['Model-1'].rootAssembly.instances['cable-1-rad-6']))


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
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[0:1]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[199:200])


for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='fixed point of rib-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].nodes[0:1]+\
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].nodes[199:200])


mdb.models['Model-1'].rootAssembly.Set(name='RIBEDGE', nodes=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[0:1]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[3:4]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[5:6]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[7:8]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[9:10]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[11:12]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[13:14]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[15:16]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[17:18]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[19:20]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[21:22]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[23:24]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[25:26]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[27:28]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[29:30]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[31:32]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[33:34]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[35:36]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[37:38]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[39:40]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[41:42]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[43:44]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[45:46]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[47:48]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[49:50]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[51:52]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[53:54]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[55:56]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[57:58]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[59:60]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[61:62]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[63:64]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[65:66]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[67:68]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[69:70]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[71:72]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[73:74]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[75:76]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[77:78]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[79:80]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[81:82]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[83:84]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[85:86]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[87:88]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[89:90]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[91:92]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[93:94]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[95:96]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[97:98]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[99:100]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[101:102]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[103:104]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[105:106]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[107:108]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[109:110]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[111:112]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[113:114]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[115:116]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[117:118]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[119:120]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[121:122]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[123:124]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[125:126]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[127:128]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[129:130]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[131:132]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[133:134]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[135:136]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[137:138]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[139:140]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[141:142]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[143:144]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[145:146]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[147:148]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[149:150]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[151:152]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[153:154]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[155:156]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[157:158]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[159:160]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[161:162]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[163:164]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[165:166]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[167:168]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[169:170]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[171:172]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[173:174]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[175:176]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[177:178]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[179:180]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[181:182]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[183:184]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[185:186]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[187:188]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[189:190]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[191:192]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[193:194]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[195:196]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[197:198]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[199:200]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[206:208]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[212:214]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[218:220]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[224:226]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[230:232]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[236:238]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[242:244]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[248:250]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[254:256]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[260:262]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[266:268]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[272:274]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[278:280]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[284:286]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[290:292]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[296:298]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[302:304]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[308:310]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[314:316]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[320:322]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[326:328]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[332:334]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[338:340]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[344:346]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[350:352]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[356:358]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[362:364]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[368:370]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[374:376]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[380:382]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[386:388]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[392:394]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[398:400]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[404:406]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[410:412]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[416:418]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[422:424]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[428:430]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[434:436]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[440:442]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[446:448]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[452:454]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[458:460]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[464:466]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[470:472]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[476:478]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[482:484]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[488:490]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[494:496]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[500:502]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[506:508]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[512:514]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[518:520]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[524:526]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[530:532]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[536:538]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[542:544]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[548:550]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[554:556]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[560:562]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[566:568]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[572:574]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[578:580]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[584:586]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[590:592]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[596:598]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[602:604]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[608:610]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[614:616]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[620:622]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[626:628]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[632:634]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[638:640]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[644:646]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[650:652]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[656:658]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[662:664]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[668:670]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[674:676]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[680:682]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[686:688]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[692:694]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[698:700]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[704:706]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[710:712]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[716:718]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[722:724]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[728:730]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[734:736]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[740:742]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[746:748]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[752:754]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[758:760]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[764:766]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[770:772]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[776:778]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[782:784]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[788:790]+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[794:796]) 


mdb.models['Model-1'].rootAssembly.Surface(name='Surf-membrane-1', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )
for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-membrane-'+ str(i), side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1-rad-'+str(i)].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) )    



mdb.models['Model-1'].rootAssembly.Surface(name='Surf-rib-1', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-2'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-rib-2', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-2'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-3'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-rib-3', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-3'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-4'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-rib-4', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-4'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-5'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-rib-5', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-5'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-6'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-rib-6', side12Elements=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
                                                            xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber)+\
    mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-6'].elements.getByBoundingBox(xMin=-LargeNumber,yMin=y[1],zMin=-LargeNumber,
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
    mdb.models['Model-1'].rootAssembly.sets['riballnodes-6'], 
    name='Tie6-1', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-6']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Tie6-1'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.sets['riballnodes-1'], 
    name='Tie6-2', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.sets['Membraneallnodes-6']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Tie6-2'].setValues(positionTolerance=
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


mdb.models['Model-1'].constraints['Tie1-2'].swapSurfaces()
mdb.models['Model-1'].constraints['Tie2-2'].swapSurfaces()
mdb.models['Model-1'].constraints['Tie3-2'].swapSurfaces()
mdb.models['Model-1'].constraints['Tie4-2'].swapSurfaces()
mdb.models['Model-1'].constraints['Tie5-2'].swapSurfaces()
mdb.models['Model-1'].constraints['Tie6-2'].swapSurfaces()    

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
    -1, endStep=1, fileName='Iteration-' + str(ind), instances=(
    mdb.models[itrModelName].rootAssembly.instances['rib-1'],  
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

    for i in range(1,m+1):
        mdb.models[itrModelName].SurfaceToSurfaceContactStd(adjustMethod=NONE, 
            clearanceRegion=None, createStepName='Step-1', datumAxis=None, 
            initialClearance=OMIT, interactionProperty='IntProp-1', master=
            mdb.models[itrModelName].rootAssembly.surfaces['Surf-membrane-'+str(i)], name=
            'Int-'+str(i), slave=mdb.models[itrModelName].rootAssembly.surfaces['Surf-rib-'+str(i)], 
            sliding=FINITE, thickness=OFF)


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
    
    ODB.close()


    mdb.models[itrModelName].keywordBlock.synchVersions(storeNodesAndElements=
        False)
    mdb.models[itrModelName].keywordBlock.replace(212, '\n')
    mdb.models[itrModelName].keywordBlock.replace(436, '\n')
    mdb.models[itrModelName].keywordBlock.synchVersions(storeNodesAndElements=
        False)
    mdb.models[itrModelName].keywordBlock.replace(642, '\n')
    mdb.models[itrModelName].keywordBlock.replace(850, '\n')
    mdb.models[itrModelName].keywordBlock.synchVersions(storeNodesAndElements=
        False)
    mdb.models[itrModelName].keywordBlock.replace(1056, '\n')
    mdb.models[itrModelName].keywordBlock.replace(1264, '\n')


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
