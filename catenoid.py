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
Em = 250.0;
nu_m = 0.34;
alpha_m = 1e-05;
sigma = 2.5;
delta_T_Membrane = sigma * (1-nu_m)/Em/alpha_m


# Part membrane
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(0.0, 40000.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='membrane', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['membrane'].BaseShellExtrude(depth=18339.5, draftAngle=-90.0+
    29.81742634, sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].parts['membrane'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=YZPLANE)
mdb.models['Model-1'].parts['membrane'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['membrane'].datums[2], faces=
    mdb.models['Model-1'].parts['membrane'].faces)

# Part cables
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(0.0, 40000.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-2', type=
    DISCRETE_RIGID_SURFACE)
mdb.models['Model-1'].parts['Part-2'].BaseShellExtrude(depth=1.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].parts['Part-2'].ReferencePoint(point=
    mdb.models['Model-1'].parts['Part-2'].InterestingPoint(
    mdb.models['Model-1'].parts['Part-2'].edges[0], CENTER))
del mdb.models['Model-1'].sketches['__profile__']


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(0.0, 8000.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-3', type=
    DISCRETE_RIGID_SURFACE)
mdb.models['Model-1'].parts['Part-3'].BaseShellExtrude(depth=1.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].parts['Part-3'].ReferencePoint(point=
    mdb.models['Model-1'].parts['Part-3'].InterestingPoint(
    mdb.models['Model-1'].parts['Part-3'].edges[0], CENTER))
del mdb.models['Model-1'].sketches['__profile__']




# Materials
mdb.models['Model-1'].Material(name='membrane')
mdb.models['Model-1'].materials['membrane'].Elastic(table=((Em, 
    nu_m), ))
mdb.models['Model-1'].materials['membrane'].Expansion(table=((alpha_m, 
    ), ))



# Property
mdb.models['Model-1'].MembraneSection(material='membrane', name='Section-1', 
    poissonDefinition=DEFAULT, thickness=1.0, thicknessField='', thicknessType=
    UNIFORM)
mdb.models['Model-1'].parts['membrane'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts['membrane'].faces), sectionName='Section-1', thicknessAssignment=FROM_SECTION)

mdb.models['Model-1'].parts['membrane'].DatumCsysByThreePoints(coordSysType=
    CYLINDRICAL, name='Datum csys-1', origin=(0.0, 0.0, 18339.5), point1=
    mdb.models['Model-1'].parts['membrane'].vertices.findAt((0.0, -7999.999998, 
    18339.5), ), point2=
    mdb.models['Model-1'].parts['membrane'].InterestingPoint(
    mdb.models['Model-1'].parts['membrane'].edges.findAt((5656.854248, 
    5656.854248, 18339.5), ), MIDDLE))
mdb.models['Model-1'].parts['membrane'].MaterialOrientation(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, fieldName='', localCsys=
    mdb.models['Model-1'].parts['membrane'].datums[5], orientationType=SYSTEM, 
    region=Region(faces=mdb.models['Model-1'].parts['membrane'].faces))


# Assembly
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
    part=mdb.models['Model-1'].parts['Part-2'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-3-1', 
    part=mdb.models['Model-1'].parts['Part-3'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='membrane-1', 
    part=mdb.models['Model-1'].parts['membrane'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-3-1', ), 
    vector=(0.0, 0.0, 18339.5))

# Step
mdb.models['Model-1'].StaticStep(initialInc=0.1, name='Step-1', nlgeom=ON, 
    previous='Initial')
mdb.models['Model-1'].steps['Step-1'].setValues(adaptiveDampingRatio=None, 
    continueDampingFactors=False, stabilizationMagnitude=0.00000002, 
    stabilizationMethod=DAMPING_FACTOR)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'COORD'))

# Mesh
mdb.models['Model-1'].parts['membrane'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=3000.0)

mdb.models['Model-1'].parts['membrane'].setMeshControls(elemShape=TRI, regions=
    mdb.models['Model-1'].parts['membrane'].faces, technique=STRUCTURED)
mdb.models['Model-1'].parts['membrane'].generateMesh()

mdb.models['Model-1'].parts['Part-2'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=1000.0)
mdb.models['Model-1'].parts['Part-2'].generateMesh()

mdb.models['Model-1'].parts['Part-3'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=1000.0)
mdb.models['Model-1'].parts['Part-3'].generateMesh()


mdb.models['Model-1'].parts['membrane'].setElementType(elemTypes=(ElemType(
    elemCode=M3D4, elemLibrary=STANDARD), ElemType(elemCode=M3D3, 
    elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=(
    mdb.models['Model-1'].parts['membrane'].faces, ))
mdb.models['Model-1'].parts['membrane'].generateMesh()


mdb.models['Model-1'].parts['Part-2'].Set(name='frame-1', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=-999999,
                                                             xMax=999999,yMax=999999,zMax=999999) )
mdb.models['Model-1'].parts['Part-3'].Set(name='frame-2', nodes=
    mdb.models['Model-1'].parts['Part-3'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=-999999,
                                                             xMax=999999,yMax=999999,zMax=999999) )

mdb.models['Model-1'].parts['membrane'].Set(name='membrane_all', nodes=
    mdb.models['Model-1'].parts['membrane'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=-999999,
                                                             xMax=999999,yMax=999999,zMax=999999) )
mdb.models['Model-1'].parts['membrane'].Set(name='membrane_bottomside', nodes=
    mdb.models['Model-1'].parts['membrane'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=0,
                                                             xMax=999999,yMax=999999,zMax=0) )
mdb.models['Model-1'].parts['membrane'].Set(name='membrane_topside', nodes=
    mdb.models['Model-1'].parts['membrane'].nodes.getByBoundingBox(xMin=-999999,yMin=-999999,zMin=18339.5,
                                                             xMax=999999,yMax=999999,zMax=18339.5) )


# Load conditions
mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
    name='BC-1', region=
    mdb.models['Model-1'].rootAssembly.instances['membrane-1'].sets['membrane_bottomside'])
mdb.models['Model-1'].EncastreBC(createStepName='Initial', localCsys=None, 
    name='BC-2', region=
    mdb.models['Model-1'].rootAssembly.instances['membrane-1'].sets['membrane_topside'])
mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(0.0, ), name='Predefined Field-1', region=
    mdb.models['Model-1'].rootAssembly.instances['membrane-1'].sets['membrane_all'])
mdb.models['Model-1'].predefinedFields['Predefined Field-1'].setValuesInStep(
    magnitudes=(-delta_T_Membrane, ), stepName='Step-1')


mdb.models['Model-1'].rootAssembly.regenerate()


#job 
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


# ---------------------------------------------------------------
# Iteration starts here
# ---------------------------------------------------------------

n = 300 # iteration times
results=np.zeros((n,2))
for ind in range(1, n):

    itrModelName = 'Iteration-' + str(ind+1)
    
    mdb.Model(name=itrModelName, objectToCopy=mdb.models['Model-1'])
    
    ODB = openOdb(path = 'Iteration-' + str(ind) + '.odb')


    firstFrame = ODB.steps['Step-1'].frames[-1]
    coordinates = firstFrame.fieldOutputs['COORD']

    
    
    # Set new coordinates for membrane
    MembraneNodes=ODB.rootAssembly.instances['MEMBRANE-1'].nodeSets['MEMBRANE_ALL']
    MembraneCoordField = coordinates.getSubset(region=MembraneNodes)
    MembraneCoord = MembraneCoordField.values
    NewCoord_Membrane = np.zeros((len(mdb.models[itrModelName].rootAssembly.instances['membrane-1'].nodes), 3))
    
    for i in range(len(mdb.models[itrModelName].rootAssembly.instances['membrane-1'].nodes)):
        NewCoord_Membrane[i][0]=np.float64(MembraneCoord[i].data[0])
        NewCoord_Membrane[i][1]=np.float64(MembraneCoord[i].data[1])
        NewCoord_Membrane[i][2]=np.float64(MembraneCoord[i].data[2])
    
    #SET THE NEW COORDINATES
    mdb.models[itrModelName].rootAssembly.editNode(
        nodes=mdb.models[itrModelName].rootAssembly.instances['membrane-1'].nodes,
        coordinates=NewCoord_Membrane)
    
    mdb.models[itrModelName].rootAssembly.regenerate()

    ODB.close()
    
    
    
    # job 
    jobName = 'Iteration-' + str(ind+1) 
    mdb.Job(model=itrModelName, name=jobName)
    # mdb.saveAs(pathName=jobName)#no cae


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