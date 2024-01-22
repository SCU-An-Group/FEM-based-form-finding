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

Mdb()

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# Parameters
Ec = 20000.0;  # MPa
nu_c = 0.3;
alpha_c = -2.0e-06; 
Cable_Area = pi*(1.0/2.0)**2;
Cable_T = 50.0; # N
delta_T_Cable = Cable_T/Cable_Area/Ec/alpha_c;
thickness_membrane = 0.0265;


Em = 2170.0;
nu_m = 0.34;
alpha_m = 29.0*1e-06;
sigma = 1.0; # MPa
delta_T_Membrane = sigma * (1-nu_m)/Em/alpha_m
 
L = 100.0; #size
element_number = 20; #mesh size

# Part membrane
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), point2=(L, L))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Membrane', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Membrane'].BaseShell(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

Plane_1 = mdb.models['Model-1'].parts['Membrane'].DatumPlaneByPrincipalPlane(offset=L/2.0, 
    principalPlane=YZPLANE)
Plane_2 = mdb.models['Model-1'].parts['Membrane'].DatumPlaneByPrincipalPlane(offset=L/2.0, 
    principalPlane=XZPLANE)
mdb.models['Model-1'].parts['Membrane'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Membrane'].datums[Plane_1.id], faces=
    mdb.models['Model-1'].parts['Membrane'].faces)
mdb.models['Model-1'].parts['Membrane'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Membrane'].datums[Plane_2.id], faces=
    mdb.models['Model-1'].parts['Membrane'].faces)


mdb.models['Model-1'].parts['Membrane'].Set(name='Set-MiddlePoint', vertices=
    mdb.models['Model-1'].parts['Membrane'].vertices.findAt(((0.0, L/2.0, 0.0), )))

# Part cables
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
    0.0, L))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, L), point2=(
    L, L))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(L, L), point2=(
    L, 0.0))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(L, 0.0), point2=(
    0.0, 0.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Cables', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Cables'].BaseWire(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

Plane_1 = mdb.models['Model-1'].parts['Cables'].DatumPlaneByPrincipalPlane(offset=L/2.0, 
    principalPlane=YZPLANE)
Plane_2 = mdb.models['Model-1'].parts['Cables'].DatumPlaneByPrincipalPlane(offset=L/2.0, 
    principalPlane=XZPLANE)
mdb.models['Model-1'].parts['Cables'].PartitionEdgeByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Cables'].datums[Plane_1.id], edges=
    mdb.models['Model-1'].parts['Cables'].edges)
mdb.models['Model-1'].parts['Cables'].PartitionEdgeByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Cables'].datums[Plane_2.id], edges=
    mdb.models['Model-1'].parts['Cables'].edges)

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

# Property
mdb.models['Model-1'].MembraneSection(material='Material-membrane', name=
    'Section-membrane', poissonDefinition=DEFAULT, thickness=thickness_membrane, 
    thicknessField='', thicknessType=UNIFORM)
mdb.models['Model-1'].TrussSection(area=Cable_Area, material='Material-cables', name=
    'Section-cables')

mdb.models['Model-1'].parts['Cables'].Set(edges=
    mdb.models['Model-1'].parts['Cables'].edges, 
    name='Set-Cables')
mdb.models['Model-1'].parts['Cables'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Cables'].sets['Set-Cables'], sectionName=
    'Section-cables', thicknessAssignment=FROM_SECTION)

mdb.models['Model-1'].parts['Membrane'].Set(faces=
    mdb.models['Model-1'].parts['Membrane'].faces, name='Set-Membrane')
mdb.models['Model-1'].parts['Membrane'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Membrane'].sets['Set-Membrane'], sectionName=
    'Section-membrane', thicknessAssignment=FROM_SECTION)

# Assembly
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Cables-1', 
    part=mdb.models['Model-1'].parts['Cables'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Membrane-1', 
    part=mdb.models['Model-1'].parts['Membrane'])

# Step
mdb.models['Model-1'].StaticStep(initialInc=0.01, name='Step-1', nlgeom=ON, 
    previous='Initial')
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'COORD'))

# Mesh
mdb.models['Model-1'].parts['Membrane'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=L/element_number)
mdb.models['Model-1'].parts['Membrane'].generateMesh()

mdb.models['Model-1'].parts['Membrane'].setElementType(elemTypes=(ElemType(
    elemCode=M3D4, elemLibrary=STANDARD, secondOrderAccuracy=OFF), ElemType(
    elemCode=M3D3, elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['Membrane'].faces, ))

mdb.models['Model-1'].parts['Membrane'].Set(name='Set-Membrane-Edgenodes_1', nodes=
    mdb.models['Model-1'].parts['Membrane'].nodes.getByBoundingBox(xMin=0,yMin=0,zMin=0,
                                                             xMax=0,yMax=L,zMax=0) )

mdb.models['Model-1'].parts['Membrane'].Set(name='Set-Membrane-Edgenodes_2', nodes=
    mdb.models['Model-1'].parts['Membrane'].nodes.getByBoundingBox(xMin=0,yMin=L,zMin=0,
                                                             xMax=L,yMax=L,zMax=0) )

mdb.models['Model-1'].parts['Membrane'].Set(name='Set-Membrane-Edgenodes_3', nodes=
    mdb.models['Model-1'].parts['Membrane'].nodes.getByBoundingBox(xMin=L,yMin=0,zMin=0,
                                                             xMax=L,yMax=L,zMax=0) )

mdb.models['Model-1'].parts['Membrane'].Set(name='Set-Membrane-Edgenodes_4', nodes=
    mdb.models['Model-1'].parts['Membrane'].nodes.getByBoundingBox(xMin=0,yMin=0,zMin=0,
                                                             xMax=L,yMax=0,zMax=0) )


mdb.models['Model-1'].parts['Cables'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=L/element_number)
mdb.models['Model-1'].parts['Cables'].generateMesh()
mdb.models['Model-1'].parts['Cables'].setElementType(elemTypes=(ElemType(
    elemCode=T3D2, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].parts['Cables'].edges, ))

# Set
mdb.models['Model-1'].parts['Cables'].Set(name='Set-Cables-Edgenodes_1', nodes=
    mdb.models['Model-1'].parts['Cables'].nodes.getByBoundingBox(xMin=0,yMin=0,zMin=0,
                                                             xMax=0,yMax=L,zMax=0) )

mdb.models['Model-1'].parts['Cables'].Set(name='Set-Cables-Edgenodes_2', nodes=
    mdb.models['Model-1'].parts['Cables'].nodes.getByBoundingBox(xMin=0,yMin=L,zMin=0,
                                                             xMax=L,yMax=L,zMax=0) )

mdb.models['Model-1'].parts['Cables'].Set(name='Set-Cables-Edgenodes_3', nodes=
    mdb.models['Model-1'].parts['Cables'].nodes.getByBoundingBox(xMin=L,yMin=0,zMin=0,
                                                             xMax=L,yMax=L,zMax=0) )
mdb.models['Model-1'].parts['Cables'].Set(name='Set-Cables-Edgenodes_4', nodes=
    mdb.models['Model-1'].parts['Cables'].nodes.getByBoundingBox(xMin=0,yMin=0,zMin=0,
                                                             xMax=L,yMax=0,zMax=0) )


# Interaction 
mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.instances['Cables-1'].sets['Set-Cables-Edgenodes_1']
    , name='Constraint-1', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Set-Membrane-Edgenodes_1']
    , thickness=ON, tieRotations=ON)

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.instances['Cables-1'].sets['Set-Cables-Edgenodes_2']
    , name='Constraint-2', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Set-Membrane-Edgenodes_2']
    , thickness=ON, tieRotations=ON)

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.instances['Cables-1'].sets['Set-Cables-Edgenodes_3']
    , name='Constraint-3', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Set-Membrane-Edgenodes_3']
    , thickness=ON, tieRotations=ON)

mdb.models['Model-1'].Tie(adjust=ON, master=
    mdb.models['Model-1'].rootAssembly.instances['Cables-1'].sets['Set-Cables-Edgenodes_4']
    , name='Constraint-4', positionToleranceMethod=COMPUTED, slave=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Set-Membrane-Edgenodes_4']
    , thickness=ON, tieRotations=ON)


# Load conditions
mdb.models['Model-1'].parts['Cables'].Set(name='Set-Vertice-nodes', vertices=
    mdb.models['Model-1'].parts['Cables'].vertices.findAt(((L, 0.0, 0.0), 
    ), ((L, L, 0.0), ), ((0.0, L, 0.0), ), ((0.0, 0.0, 0.0), ), 
    ))

mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
    name='BC-1', region=
    mdb.models['Model-1'].rootAssembly.instances['Cables-1'].sets['Set-Vertice-nodes'])

mdb.models['Model-1'].parts['Membrane'].Set(name='Set-Membrane-allnodes', 
    nodes=mdb.models['Model-1'].parts['Membrane'].nodes.getByBoundingBox(xMin=-9999,yMin=-9999,zMin=-9999,
                                                             xMax=9999,yMax=9999,zMax=9999) )

mdb.models['Model-1'].parts['Cables'].Set(name='Set-Cables-allnodes', 
    nodes=mdb.models['Model-1'].parts['Cables'].nodes.getByBoundingBox(xMin=-9999,yMin=-9999,zMin=-9999,
                                                             xMax=9999,yMax=9999,zMax=9999) )
mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(0.0, ), name='Predefined Field-1', region=
    mdb.models['Model-1'].rootAssembly.instances['Cables-1'].sets['Set-Cables-allnodes'])


mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(0.0, ), name='Predefined Field-2', region=
    mdb.models['Model-1'].rootAssembly.instances['Membrane-1'].sets['Set-Membrane-allnodes'])

#
mdb.models['Model-1'].predefinedFields['Predefined Field-1'].setValuesInStep(
    magnitudes=(-delta_T_Cable, ), stepName='Step-1')

mdb.models['Model-1'].predefinedFields['Predefined Field-2'].setValuesInStep(
    magnitudes=(-delta_T_Membrane, ), stepName='Step-1')

mdb.models['Model-1'].rootAssembly.regenerate()


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

n = 300 # iteration times
for ind in range(1, n):

    itrModelName = 'Iteration-' + str(ind+1)
    
    mdb.Model(name=itrModelName, objectToCopy=mdb.models['Model-1'])
    
    ODB = openOdb(path = 'Iteration-' + str(ind) + '.odb')


    firstFrame = ODB.steps['Step-1'].frames[-1]
    coordinates = firstFrame.fieldOutputs['COORD']


    
    # Set new coordinates for cables
    CablesNodes=ODB.rootAssembly.instances['CABLES-1'].nodeSets['SET-CABLES']
    CablesCoordField = coordinates.getSubset(region=CablesNodes)
    CablesCoord = CablesCoordField.values
    NewCoord_Cables = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['Cables-1'].nodes), 3))
    
    for i in range(len(mdb.models[itrModelName].rootAssembly.instances['Cables-1'].nodes)):
        NewCoord_Cables[i][0]=np.float64(CablesCoord[i].data[0])
        NewCoord_Cables[i][1]=np.float64(CablesCoord[i].data[1])
        NewCoord_Cables[i][2]=np.float64(CablesCoord[i].data[2])
    
    #SET THE NEW COORDINATES
    mdb.models[itrModelName].rootAssembly.editNode(
        nodes=mdb.models[itrModelName].rootAssembly.instances['Cables-1'].nodes,
        coordinates=NewCoord_Cables)
    
    
    # Set new coordinates for membrane
    MembraneNodes=ODB.rootAssembly.instances['MEMBRANE-1'].nodeSets['SET-MEMBRANE-ALLNODES']
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

