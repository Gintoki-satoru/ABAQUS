# Code for assigning material orientations for beam stiffeners (elements).
# This code has been developed by Dr.-Ing. Ahmad Alhaj Ahmad LSM, TU-Darmstadt, Germany 
from math import*
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
from numpy import* # Numeric module is working with ABAQUS 6.5. This version doesn't know numpy module. 
#from Numeric import *
from array import*
from abaqus import *
import regionToolset
import numpy as np
import math as math
import os

path_modules = 'U:\\Sachdeva\\MT_Nair\\ABAQUS\\MT'
os.chdir(path_modules)
# Input:

Lx= 600.0               # the length of the sequare panel.
W= 600.0                # The width of the panel.

SkinCorners=Lx/2.0     # coordinates of the corners of the panels to be used to draw the panel taking into account that the origin of coordinate system is located a the panel center.

NL=8	                 # number of layers for the skin
T=0.3225           # layer thickness
panelThickness=T*NL      # the total thickness of the skin

w_beamCrossSection= 6.192
h_beamCrossSection= 30.0-panelThickness/2.0
#--------------------------------------------------------------------------------------------------------

# Engineering constants 'MPa' for the panel:
E1=150.0e3  
E2=9080.0   
E3=9080.0   
nu12=0.32   
nu13=0.32   
nu23=0.458  
G12=5290    
G13=5290    
G23=4000    

LAM=[E1,E2,nu12,G12,G13,G23]

myModel=mdb.models['Model-1']
mat=myModel.Material(name='carbon-epoxy')
mat.Elastic(type=LAMINA, table=(tuple(LAM),))
myModel.materials['carbon-epoxy'].Density(table=((0.00135, ), ))

#Beam Material and Section defintions:

myModel.IProfile(b1=w_beamCrossSection, b2=w_beamCrossSection, h=h_beamCrossSection, l=-panelThickness/2.0, name='Profile-1', t1=h_beamCrossSection/2.0, t2=h_beamCrossSection/2.0, t3=w_beamCrossSection)
myModel.Material(name='Material-1')
myModel.materials['Material-1'].Elastic(table=((E1,E2,E3, nu12, nu13, nu23,G12, G13, G23), ), type=ENGINEERING_CONSTANTS)
myModel.materials['Material-1'].Density(table=((0.0016, ), ))
myModel.BeamSection(consistentMassMatrix=False, integration=DURING_ANALYSIS, material='Material-1', name='beam', profile='Profile-1')

##########################################################################################################################################################################################
# Draw the flat skin:
myModel.ConstrainedSketch(name='SkinSketch', sheetSize=4.0*Lx)
myModel.sketches['SkinSketch'].rectangle(point1=(Lx/2.0, W/2.0), point2=(-Lx/2.0, -W/2.0))
myModel.Part(dimensionality=THREE_D, name='Skin', type=  DEFORMABLE_BODY)
myModel.parts['Skin'].BaseShell(sketch= myModel.sketches['SkinSketch'])

skinPart=myModel.parts['Skin'] 
###############################################################################################################################################

mySpline=myModel.ConstrainedSketch(name='Besier', sheetSize=4.0*Lx)
execfile('Besier_curve.py')
pointsxy = []                          
for j in range(len(XP)): 
   pointsxy.append([XP[j],YP[j]])

mySpline.Spline(points=pointsxy)
###################################################################################################################################################

#Partition for family of stringer-1:
edges=skinPart.getFeatureEdges('Shell planar-1')
edgeAtPoint=edges.getClosest(coordinates=((Lx/2.0,0.25*W,0.0),))
RightEdge=edgeAtPoint[0][0]
facesToBePartitioned=skinPart.getFeatureFaces('Shell planar-1')
faceAsSketchPlane=skinPart.getFeatureFaces('Shell planar-1')[0]

mdb.models['Model-1'].ConstrainedSketch(gridSpacing=70.71, name='__profile__',sheetSize=2828.42, transform=mdb.models['Model-1'].parts['Skin'].MakeSketchTransform(
    sketchPlane=faceAsSketchPlane,sketchPlaneSide=SIDE1,sketchUpEdge=RightEdge, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
skinPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])

mdb.models['Model-1'].sketches['__profile__'].retrieveSketch(sketch=mdb.models['Model-1'].sketches['Besier'])

skinPart.PartitionFaceBySketch(faces=facesToBePartitioned, sketch=mdb.models['Model-1'].sketches['__profile__'], sketchUpEdge=RightEdge)
#

# Create the Stringer family-1:
edgesFamily1=skinPart.getFeatureEdges('Partition face-1')
skinPart.Stringer(name='Stringer-1', edges=edgesFamily1)

#Section Assignment:
skinPart.SectionAssignment(offset=0.0,offsetField='', offsetType=MIDDLE_SURFACE, region=Region(stringerEdges=((
    'Stringer-1', edgesFamily1), )), sectionName='beam', thicknessAssignment=FROM_SECTION)
#Beam orientation of the stringers:
skinPart.assignBeamSectionOrientation(method=N1_COSINES, n1=(0.0, 1.0, 0.0), region=Region(stringerEdges=(('Stringer-1', edgesFamily1), )))
#####################################################################################################################################################
#Meshing:
allFaces=skinPart.getFeatureFaces('Shell planar-1')
skinPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=Lx/60.0)

skinPart.setMeshControls(elemShape=TRI, regions= allFaces)
s=skinPart.Set(name='panel', faces=skinPart.faces)
skinPart.setElementType(elemTypes=(ElemType( elemCode=S4R, elemLibrary=STANDARD), ElemType(elemCode=S3,elemLibrary=STANDARD, secondOrderAccuracy=OFF)), regions=s)
skinPart.generateMesh()
####################################################################################################################################################
# define the material orientations:
setFaces=skinPart.Set(name='Skin_set', faces=skinPart.faces)
skinPart.DatumCsysByThreePoints(coordSysType=CARTESIAN,name='CSYS-orientation', origin=(-Lx/2.0, -W/2.0, 0.0), point1=(Lx/2.0,-W/2.0, 0.0), point2=(0.0, 0.0, 0.0))
k=skinPart.datums.keys()
dcsys=skinPart.datums[k[-1]]
#mat.UserOutputVariables(1)
skinPart.MaterialOrientation(localCsys=dcsys, axis=AXIS_3, angle=0., region=setFaces)

############################################################################################################################################# 
# The beam elements associated with the stringers:
elementsOfStringer1=skinPart.stringers['Stringer-1'].elements

k=0
for element in skinPart.elements:
    k=k+1
    elemLabel=element.label
    ELtype=element.type
    nodalConnect=element.connectivity
    # convert the tuple to a list:
    nNodes=len(nodalConnect)
    x = [0.0]*nNodes
    y = [0.0]*nNodes
    z = [0.0]*nNodes
    j=-1
    for i in nodalConnect:
        j=j+1
        nodeCoord=skinPart.nodes[i].coordinates
        # be carefull and make sure that y and z ( or 1 an 2 indices) are compatible with your model coord sys.
        x[j]=nodeCoord[0]
        y[j]=nodeCoord[1]
        z[j]=nodeCoord[2]
    if ELtype==B31:
        Cx=(x[0]+x[1])/2.0
        Cy=(y[0]+y[1])/2.0
        Cz=(z[0]+z[1])/2.0   
        n='ELstr1_%02d' %(k)
        s=skinPart.SetFromElementLabels(name=n,elementLabels=(elemLabel,))
        localAngle=math.atan((y[1]-y[0])/(x[1]-x[0]))*180/np.pi
        skinPart.MaterialOrientation(additionalRotationField='', additionalRotationType=ROTATION_ANGLE, angle=
		localAngle, axis=AXIS_3, fieldName='', localCsys=dcsys, orientationType=SYSTEM,region=s, stackDirection=STACK_3)
 
 
#######################################################################################################################################