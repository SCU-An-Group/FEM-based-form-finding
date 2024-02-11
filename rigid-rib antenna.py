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

DampingFactor = 1e-10


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
pretension_m = 0.1;
delta_T_m = pretension_m * (1-nu_m)/Em/alpha_m

#rib
density_r = 7.86e-9
Er = 206000
nu_r = 0.3



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


# create rib
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='rib', type=
    DEFORMABLE_BODY)

mdb.models['Model-1'].parts['rib'].ReferencePoint(point=(0.0, 0.0, 0.0))
mdb.models['Model-1'].parts['rib'].Set(name='Set-RP', referencePoints=(
    mdb.models['Model-1'].parts['rib'].referencePoints[1], ))


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
rib_edges=mdb.models['Model-1'].parts['rib'].edges.getByBoundingBox(xMin=-LargeNumber,yMin=-LargeNumber,zMin=-LargeNumber,
                                                             xMax=LargeNumber,yMax=LargeNumber,zMax=LargeNumber) 
mdb.models['Model-1'].parts['rib'].Set(edges=rib_edges,name=setname)

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
mdb.models['Model-1'].materials['Material-membrane'].expansion.setValues(table=
    ((alpha_m, alpha_m, 0.0), ), type=ORTHOTROPIC)

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
    CYLINDRICAL, line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0), name=
    'Datum csys-1', origin=(0.0, 0.0, 0.0))
mdb.models['Model-1'].parts['Membrane'].MaterialOrientation(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, fieldName='', localCsys=
    mdb.models['Model-1'].parts['Membrane'].datums[279], orientationType=SYSTEM
    , region=mdb.models['Model-1'].parts['Membrane'].sets['AllMem'])


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


#rib 
rib_radius = 15
theta_np = np.linspace(0. * pi/180.0, -56. * pi/180.0, 100)
rib_thickness = 0.1
points = []
for theta in theta_np:
    points.append((rib_radius * np.cos(theta)-15, rib_radius * np.sin(theta), rib_thickness))
mdb.models['Model-1'].ArbitraryProfile(name='C-section', table=(points))

mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
    DURING_ANALYSIS, material='Material-ribs', name='Section-3', poissonRatio=
    0.0, profile='C-section', temperatureVar=LINEAR)
mdb.models['Model-1'].parts['rib'].SectionAssignment(offset=0.0, offsetField=''
    , offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['rib'].sets['AllRib'], sectionName='Section-3', 
    thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['rib'].assignBeamSectionOrientation(method=
    N1_COSINES, n1=(0.0, 0.0, -1.0), region=
    mdb.models['Model-1'].parts['rib'].sets['AllRib'])
mdb.models['Model-1'].sections['Section-3'].TransverseShearBeam(k13=None, k23=
    None, scfDefinition=COMPUTED, slendernessCompensation=COMPUTED)



# #assembly

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Membrane-1', 
    part=mdb.models['Model-1'].parts['Membrane'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='cable-1', part=
    mdb.models['Model-1'].parts['cable'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='rib-1', part=
    mdb.models['Model-1'].parts['rib'])
mdb.models['Model-1'].rootAssembly.RadialInstancePattern(axis=(0.0, 1.0, 0.0), 
    instanceList=('Membrane-1', 'cable-1', 'rib-1'), number=16, point=(0.0, 
    0.0, 0.0), totalAngle=360.0)

#step
mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None, 
    continueDampingFactors=False, initialInc=0.01, name='Step-1', nlgeom=ON, 
    previous='Initial', stabilizationMagnitude=DampingFactor, stabilizationMethod=
    DAMPING_FACTOR)
mdb.models['Model-1'].steps['Step-1'].setValues(minInc=1e-09)
mdb.models['Model-1'].steps['Step-1'].setValues(maxNumInc=1000)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'COORD'))


# #mesh

mdb.models['Model-1'].parts['Membrane'].setMeshControls(elemShape=QUAD, regions=
    mdb.models['Model-1'].parts['Membrane'].faces)
mdb.models['Model-1'].parts['Membrane'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=10.0)
mdb.models['Model-1'].parts['Membrane'].generateMesh()
mdb.models['Model-1'].parts['Membrane'].setElementType(elemTypes=(ElemType(
    elemCode=M3D4R, elemLibrary=STANDARD), ElemType(elemCode=M3D3, 
    elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=Region(faces=mem_faces), )

mdb.models['Model-1'].parts['cable'].seedEdgeByNumber(constraint=FINER, edges=
    mdb.models['Model-1'].parts['cable'].edges, number=5)
mdb.models['Model-1'].parts['cable'].generateMesh()
mdb.models['Model-1'].parts['cable'].setElementType(elemTypes=(ElemType(
    elemCode=T3D2, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].parts['cable'].edges, ))

mdb.models['Model-1'].parts['rib'].seedPart(deviationFactor=0.1, minSizeFactor=
    0.1, size=10.0)
mdb.models['Model-1'].parts['rib'].setElementType(elemTypes=(ElemType(
    elemCode=B31, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].parts['rib'].edges, ))
mdb.models['Model-1'].parts['rib'].generateMesh()

mdb.models['Model-1'].rootAssembly.makeIndependent(instances=(
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
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].nodes[39:40])

for i in range(2,m+1):
    mdb.models['Model-1'].rootAssembly.Set(name='fixed point of rib-'+ str(i), nodes=
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].nodes[39:40])



#INTERACTION

mdb.models['Model-1'].RigidBody(bodyRegion=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].sets['AllRib'], name=
    'RB-1', refPointRegion=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].sets['Set-RP'])        

for i in range(2,m+1):
    mdb.models['Model-1'].RigidBody(bodyRegion=
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].sets['AllRib'], name=
        'RB-'+ str(i), refPointRegion=
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].sets['Set-RP'])   


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

mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
    name='BC-1', region=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].sets['Set-RP'])

for i in range(2,m+1):
    mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
        name='BC-'+ str(i), region=
        mdb.models['Model-1'].rootAssembly.instances['rib-1-rad-'+ str(i)].sets['Set-RP'])

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

    mdb.models[itrModelName].rootAssembly.regenerate()

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


    mdb.models[itrModelName].rootAssembly.regenerate()
    
    ODB.close()


    # job 
    jobName = 'Iteration-' + str(ind+1) 
    mdb.Job(model=itrModelName, name=jobName)
    # mdb.saveAs(pathName=jobName)


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
