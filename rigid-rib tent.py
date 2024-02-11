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

# Parameters
Ec = 40000;
nu_c = 0.35
alpha_c = -3.0e-06;
Cable_Area = pi*(2.0/2.0)**2
Cable_T = 1;
delta_T_Cable = Cable_T/Cable_Area/Ec/alpha_c;

#
Em = 1.7;
nu_m = 0.3;
alpha_m = 1e-05;
sigma = 0.1;
delta_T_Membrane = sigma * (1-nu_m)/Em/alpha_m
Thickness_M = 0.08


Mdb()
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Membrane', type=
    DEFORMABLE_BODY)

mdb.models['Model-1'].parts['Membrane'].ReferencePoint(point=(0.0, 0.0, 0.0))


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
    edges=mdb.models['Model-1'].parts['Membrane'].edges.getByBoundingBox(xMin=x1[i],yMin=-9999,zMin=-9999,
                                                             xMax=x2[i+1],yMax=9999,zMax=9999) 
    mdb.models['Model-1'].parts['Membrane'].Set(edges=edges,name=setname)
    mdb.models['Model-1'].parts['Membrane'].CoverEdges(edgeList=edges,tryAnalytical=True)


#repalce faces
setname = 'mem_allfaces' 
faces=mdb.models['Model-1'].parts['Membrane'].faces.getByBoundingBox(xMin=x1[0],yMin=-9999,zMin=-9999,
                                                             xMax=x2[n-1],yMax=9999,zMax=9999) 
mdb.models['Model-1'].parts['Membrane'].Set(faces=faces,name=setname)
# mdb.models['Model-1'].parts['Membrane'].ReplaceFaces(faceList=faces, stitch=True)



# create rib
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='rib', type=
    DISCRETE_RIGID_SURFACE)

mdb.models['Model-1'].parts['rib'].ReferencePoint(point=(0.0, 0.0, 0.0))


# create rib point
name_point1 = []
name_point2 = []
for i in range(0,n):
    name_point1.append('P1_' + str(i+1))
    name_point2.append('P2_' + str(i+1))


for i in range(len(x)):
    name_point1[i] = mdb.models['Model-1'].parts['rib'].DatumPointByCoordinate(coords=(x1[i], y1[i], z1[i]))
    name_point2[i] = mdb.models['Model-1'].parts['rib'].DatumPointByCoordinate(coords=(x2[i], y2[i], z2[i]))

# create rib line
name_line1 = []
name_line2 = []

for i in range(0,n-1):
    name_line1.append('L1_' + str(i+1))
    name_line2.append('L2_' + str(i+1))
    

for i in range(len(x)-1):
    name_line1[i] = mdb.models['Model-1'].parts['rib'].WirePolyLine(mergeType=IMPRINT, meshable=
        ON, points=((mdb.models['Model-1'].parts['rib'].datums[name_point1[i].id], 
            mdb.models['Model-1'].parts['rib'].datums[name_point1[i+1].id])))
    name_line2[i] = mdb.models['Model-1'].parts['rib'].WirePolyLine(mergeType=IMPRINT, meshable=
        ON, points=((mdb.models['Model-1'].parts['rib'].datums[name_point2[i].id], 
            mdb.models['Model-1'].parts['rib'].datums[name_point2[i+1].id])))

#creat cable
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='cable', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['cable'].ReferencePoint(point=(0.0, 0.0, 0.0))


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
rib_edges=mdb.models['Model-1'].parts['rib'].edges.getByBoundingBox(xMin=-9999,yMin=-9999,zMin=-9999,
                                                             xMax=9999,yMax=9999,zMax=9999) 
mdb.models['Model-1'].parts['rib'].Set(edges=rib_edges,name=setname)

setname = 'AllMem'
mem_faces=mdb.models['Model-1'].parts['Membrane'].faces.getByBoundingBox(xMin=-9999,yMin=-9999,zMin=-9999,
                                                             xMax=9999,yMax=9999,zMax=9999) 
mdb.models['Model-1'].parts['Membrane'].Set(faces=mem_faces,name=setname)

setname = 'AllCable'
cable_edges=mdb.models['Model-1'].parts['cable'].edges.getByBoundingBox(xMin=-9999,yMin=-9999,zMin=-9999,
                                                             xMax=9999,yMax=9999,zMax=9999) 
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


#properties
# mdb.models['Model-1'].CircularProfile(name='Profile-1', r=Thickness_R)
# mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
#     DURING_ANALYSIS, material='Material-cables', name='Section-1', 
#     poissonRatio=0.0, profile='Profile-1', temperatureVar=LINEAR)
# mdb.models['Model-1'].parts['rib'].SectionAssignment(offset=0.0, offsetField=''
#     , offsetType=MIDDLE_SURFACE, region=Region(edges=rib_edges), sectionName='Section-1', 
#     thicknessAssignment=FROM_SECTION)
# mdb.models['Model-1'].parts['rib'].assignBeamSectionOrientation(method=
#     N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(edges=rib_edges))

mdb.models['Model-1'].MembraneSection(material='Material-membrane', name=
    'Section-2', poissonDefinition=DEFAULT, thickness=Thickness_M, thicknessField='', 
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
    mdb.models['Model-1'].parts['Membrane'].datums[201])
mdb.models['Model-1'].parts['Membrane'].MaterialOrientation(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, fieldName='', localCsys=
    mdb.models['Model-1'].parts['Membrane'].datums[700], orientationType=SYSTEM
    , region=Region(faces=mdb.models['Model-1'].parts['Membrane'].faces))

mdb.models['Model-1'].TrussSection(area=Cable_Area, material='Material-cables', 
    name='Section-3')
mdb.models['Model-1'].parts['cable'].SectionAssignment(offset=0.0, offsetField=
    '', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['cable'].sets['AllCable'], sectionName=
    'Section-3', thicknessAssignment=FROM_SECTION)


#assembly
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Membrane-1', 
    part=mdb.models['Model-1'].parts['Membrane'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='cable-1', part=
    mdb.models['Model-1'].parts['cable'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='rib-1', part=
    mdb.models['Model-1'].parts['rib'])
# mdb.models['Model-1'].rootAssembly.RadialInstancePattern(axis=(0.0, 1.0, 0.0), 
#     instanceList=('Membrane-1', 'cable-1', 'rib-1'), number=12, point=(0.0, 
#     0.0, 0.0), totalAngle=360.0)

#step
mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None, 
    continueDampingFactors=False, initialInc=0.1, name='Step-1', nlgeom=ON, 
    previous='Initial', stabilizationMagnitude=1e-10, stabilizationMethod=
    DAMPING_FACTOR)
mdb.models['Model-1'].steps['Step-1'].setValues(minInc=1e-09)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'COORD'))


# #mesh
mdb.models['Model-1'].parts['Membrane'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=50.0)
mdb.models['Model-1'].parts['Membrane'].setMeshControls(elemShape=QUAD, regions=
    mdb.models['Model-1'].parts['Membrane'].faces)
mdb.models['Model-1'].parts['Membrane'].generateMesh()
mdb.models['Model-1'].parts['Membrane'].setElementType(elemTypes=(ElemType(
    elemCode=M3D4R, elemLibrary=STANDARD), ElemType(elemCode=M3D3, 
    elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=Region(faces=mem_faces), )

mdb.models['Model-1'].parts['cable'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=50.0)
mdb.models['Model-1'].parts['cable'].generateMesh()
mdb.models['Model-1'].parts['cable'].setElementType(elemTypes=(ElemType(
    elemCode=T3D2, elemLibrary=STANDARD), ), regions=Region(edges=cable_edges), )

mdb.models['Model-1'].parts['rib'].seedPart(deviationFactor=0.1, minSizeFactor=
    0.1, size=100.0)
mdb.models['Model-1'].parts['rib'].generateMesh()


mdb.models['Model-1'].parts['Membrane'].Set(name='Membraneallnodes', nodes=
    mdb.models['Model-1'].parts['Membrane'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=-999999,
                                                             xMax=999999,yMax=999999,zMax=999999) )
mdb.models['Model-1'].parts['cable'].Set(name='cableallnodes', nodes=
    mdb.models['Model-1'].parts['cable'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=-999999,
                                                             xMax=999999,yMax=999999,zMax=999999) )
mdb.models['Model-1'].parts['rib'].Set(name='riballnodes', nodes=
    mdb.models['Model-1'].parts['rib'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=-999999,
                                                             xMax=999999,yMax=999999,zMax=999999) )

mdb.models['Model-1'].parts['cable'].Set(name='cable_twosides', 
    nodes=mdb.models['Model-1'].parts['cable'].nodes[0:1]+\
    mdb.models['Model-1'].parts['cable'].nodes[20:21])

mdb.models['Model-1'].parts['cable'].Set(name='cable_middle', 
    nodes=mdb.models['Model-1'].parts['cable'].nodes[1:20])

#
mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].sets['riballnodes'], 
    name='Constraint-1', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Membraneallnodes']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Constraint-1'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Membraneallnodes'], 
    name='Constraint-2', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.instances['cable-1'].sets['cableallnodes']
    , thickness=ON, tieRotations=ON)
mdb.models['Model-1'].constraints['Constraint-2'].setValues(positionTolerance=
    0.001, positionToleranceMethod=SPECIFIED)

#load
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=
    mdb.models['Model-1'].rootAssembly.instances['rib-1'].sets['riballnodes'], 
    u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(0.0, ), name='Predefined Field-1', region=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Membraneallnodes'])
mdb.models['Model-1'].predefinedFields['Predefined Field-1'].setValuesInStep(
    magnitudes=(-delta_T_Membrane, ), stepName='Step-1')

mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(0.0, ), name='Predefined Field-2', region=
    mdb.models['Model-1'].rootAssembly.instances['cable-1'].sets['cableallnodes'])
mdb.models['Model-1'].predefinedFields['Predefined Field-2'].setValuesInStep(
    magnitudes=(-delta_T_Cable, ), stepName='Step-1')

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

    itrModelName = 'Model-' + str(ind+1)

    mdb.Model(name=itrModelName, objectToCopy=mdb.models['Model-1'])
    
    ODB = openOdb(path = 'Iteration-' + str(ind) + '.odb')

    firstFrame = ODB.steps['Step-1'].frames[-1]
    coordinates = firstFrame.fieldOutputs['COORD']

 # Set new coordinates for membrane
    MembraneNodes=ODB.rootAssembly.instances['MEMBRANE-1'].nodeSets['MEMBRANEALLNODES']
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
    
    mdb.models[itrModelName].rootAssembly.regenerate()
    
    # Set new coordinates for rib
    RibNodes=ODB.rootAssembly.instances['RIB-1'].nodeSets['RIBALLNODES']
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


# Set new coordinates for cable
    CableNodes=ODB.rootAssembly.instances['CABLE-1'].nodeSets['CABLEALLNODES']
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


    ind = ind + 1

    
odb = openOdb(path = 'Iteration-' + str(n) + '.odb')
firstFrame = odb.steps['Step-1'].frames[-1]
coordinates = firstFrame.fieldOutputs['COORD']
itrModelName = 'Model-' + str(n)

RibNodes=odb.rootAssembly.instances['RIB-1'].nodeSets['ALLRIB']
RibCoordField = coordinates.getSubset(region=RibNodes)
RibCoord = RibCoordField.values
NewCoord_Rib = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['rib-1'].nodes), 3))

for i in range(len(mdb.models[itrModelName].rootAssembly.instances['rib-1'].nodes)):
    NewCoord_Rib[i][0]=np.float64(RibCoord[i].data[0])
    NewCoord_Rib[i][1]=np.float64(RibCoord[i].data[1])
    NewCoord_Rib[i][2]=np.float64(RibCoord[i].data[2])

    
    # odb = openOdb(path = 'Iteration-' + str(ind+1) + '.odb')
    # fm = -1
    # timeFrame = odb.steps['Step-1'].frames[fm]
    # readElement = odb.rootAssembly.instances['MEMBRANE-1'].elementSets['ALLMEM']
    # Stress = timeFrame.fieldOutputs['S']  # Remember to set field outputs manually
    # readElementStress = Stress.getSubset(position=CENTROID, region=readElement)

    # Membrane_Elements_S11 = np.zeros(len(readElementStress.values))
    # Membrane_Elements_S22 = np.zeros(len(readElementStress.values))

    # for i in range(len(readElementStress.values)):
    #    readElementFaceStress_S11 = readElementStress.values[i].data[0]  # data[0] - S11, data[1] - S22, data[2] - S33 , data[3] - S12
    #    readElementFaceStress_S22 = readElementStress.values[i].data[1]
    #    readElementFaceStress_EleLabel = readElementStress.values[i].elementLabel

    #    Membrane_Elements_S11[i] = readElementFaceStress_S11
    #    Membrane_Elements_S22[i] = readElementFaceStress_S22


    # prediction_stress = []

    # for i in range(len(readElementStress.values)):
    #     prediction_stress.append(sigma)


    # error = []
    # for i in range(len(readElementStress.values)):
    #     error.append(Membrane_Elements_S11[i] - prediction_stress[i])
     
     
    # print("Errors: ", error)
    # print(error)


    # squaredError = []
    # absError = []
    # for val in error:
    #     squaredError.append(val * val)#target-prediction之差平方 
    #     absError.append(abs(val))#误差绝对值


    # targetDeviation = []
    # targetMean = sum(Membrane_Elements_S11) / len(Membrane_Elements_S11)#target平均值
    # for val in Membrane_Elements_S11:
    #     targetDeviation.append((val - targetMean) * (val - targetMean))


    # results[ind]=(ind,sqrt(sum(squaredError) / len(squaredError)))

    # odb.close()    

    # ind = ind + 1        
     

    


csvFile = open("results3.csv","w")
writer = csv.writer(csvFile)
writer.writerow(["x","y"])

for value in NewCoord_Ribs:
    writer.writerow(value)
csvFile.close()