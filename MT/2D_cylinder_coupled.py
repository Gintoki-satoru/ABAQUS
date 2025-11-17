import __main__
import os
import time
import datetime as dt
import operator

from abaqus import *
from abaqusConstants import *
#session.setValues(kernelMemoryLimit=16384)
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgot
import displayGroupMdbToolset as dgmt
import connectorBehavior
import numpy as np
import math

path_modules = 'U:\\Sachdeva\\MT_Nair\\FE'

os.chdir(path_modules)

# Further packages:
# import coordinateTransformation as ct  # REPLACED with working functions

# =============================================================================
# WORKING COORDINATE TRANSFORMATION FUNCTIONS (FIXED)
# =============================================================================

def pol2cart_x(r, theta):
    """
    Convert polar coordinates to Cartesian X-coordinate.
    
    Args:
        r (float): Radius
        theta (float): Angle in radians
        
    Returns:
        float: X-coordinate (r * cos(theta))
    """
    return r * math.cos(theta)

def pol2cart_y(r, theta):
    """
    Convert polar coordinates to Cartesian Y-coordinate.
    
    Args:
        r (float): Radius  
        theta (float): Angle in radians
        
    Returns:
        float: Y-coordinate (r * sin(theta))
    """
    return r * math.sin(theta)

def cart2pol_radius(x, y):
    """
    Convert Cartesian coordinates to polar radius.
    
    Args:
        x (float or array): X-coordinate(s)
        y (float or array): Y-coordinate(s)
        
    Returns:
        float or array: Radius (sqrt(x^2 + y^2))
    """
    # Handle both scalar and array inputs
    try:
        # Try numpy first for array compatibility
        return np.sqrt(x*x + y*y)
    except:
        # Fallback to math for scalar values
        return math.sqrt(x*x + y*y)

def cart2pol_theta(x, y):
    """
    Convert Cartesian coordinates to polar angle (theta) in radians.
    
    Args:
        x (float or array): X-coordinate(s)
        y (float or array): Y-coordinate(s)
        
    Returns:
        float or array: Angle in radians (atan2(y, x))
    """
    # Handle both scalar and array inputs
    try:
        # Try numpy first for array compatibility
        return np.arctan2(y, x)
    except:
        # Fallback to math for scalar values
        return math.atan2(y, x)

# Create coordinate transformation module class for compatibility
class CoordinateTransformation:
    """
    Coordinate transformation module class providing the working functions.
    This class can be used as a drop-in replacement for the broken coordinateTransformation module.
    """
    
    @staticmethod
    def pol2cart_x(r, theta):
        """Static method version of pol2cart_x for module compatibility."""
        return r * math.cos(theta)
    
    @staticmethod  
    def pol2cart_y(r, theta):
        """Static method version of pol2cart_y for module compatibility."""
        return r * math.sin(theta)
    
    @staticmethod
    def cart2pol_radius(x, y):
        """Convert Cartesian coordinates to polar radius."""
        # Handle both scalar and array inputs
        try:
            # Try numpy first for array compatibility
            return np.sqrt(x*x + y*y)
        except:
            # Fallback to math for scalar values
            return math.sqrt(x*x + y*y)
    
    @staticmethod
    def cart2pol_theta(x, y):
        """Convert Cartesian coordinates to polar angle (theta) in radians."""
        # Handle both scalar and array inputs
        try:
            # Try numpy first for array compatibility
            return np.arctan2(y, x)
        except:
            # Fallback to math for scalar values
            return math.atan2(y, x)

# Create ct module instance for compatibility with existing code
ct = CoordinateTransformation()

analysis_Path = 'U:\\Sachdeva\\MT_Nair\\FE'

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

############################################################################
# PARAMETER ZUM ANPASSEN

# Geometrie-Endwinkel (für L-Bracket: π/2 = 90°)
angleEnd = np.pi / 2
# Auswertungswinkel für Post-Processing (für element classification)
angleEval = 0
#Einstellung, ob bei 0° ausgewertet werden soll
Null_Auswaertung = True

# Meshparameter für den gekrümmten Bereich
# Anzahl Elemente radial, tangential
# NOTE: mR must be >= 3 for post-processing to work correctly
# The post-processing algorithm extrapolates stresses to nodes and needs sufficient elements
mR, mT = 3,250#5,1000
mFlanken = mT*2
# Ratio für gekrümmten Bereich
mRRatio, mTRatio = 1, 6#25

#Meshparameter für Flanken
#mFlanken = 500  # 50
#Verteilungsparameter
mBias = mTRatio*2   
rBias = 7  

#Probekörpergeometrie
#Innenradius
ri = 152.5
#Länge der Tangente
ltan=300
#3D tiefe
length=2
heatcoeff = 9.9934e-06

compositeMaterialName = 'cfk'  # 'cfk', 'cfkDuro', 'cfk'
AngleName = 'bestSymm'  # 'bestSymm', 'bestUnsymm', 'Vergleichslaminat', 'Test'

if compositeMaterialName == 'cfk':
	if AngleName == 'Test':
		ric_plyAngle = '45 -45 30 -30'.split()#45 -45 -45 45 90 30 -30 -30 30 90
	elif AngleName == 'Thick':
		ric_plyAngle = '45 -45 -45 45 30 -30 -30 30 60 -60 -60 60 75 -75 -75 75 15 -15 -15 15 90 90 90 90'.split()
	elif AngleName == 'bestSymm':
		ric_plyAngle = '0'  #00010110111010111011011101011101101000
	elif AngleName == 'onlyB':
		ric_plyAngle = '00011011011011011011011011011011011000'  
	elif AngleName == 'bestUnsymm':
		ric_plyAngle = '00111011101110100101101110111011001000' 
	elif AngleName == 'Vergleichslaminat':
		ric_plyAngle = '00110010011001011100101001101101001110'  
	elif AngleName == 'Foto':
		ric_plyAngle = '0010110'  
	 
elif compositeMaterialName == 'gfk':
	if AngleName == 'Test':
		ric_plyAngle = '01010101010101010101010101'  
	elif AngleName == 'bestSymm':
		ric_plyAngle = '00101101101011010110110100'  
	elif AngleName == 'bestUnsymm':
		ric_plyAngle = '00110110110110110110100100' 
	elif AngleName == 'Vergleichslaminat':
		ric_plyAngle = '00110010010011001001001100'  


############################################################################
# Bestimmung des Zeitpunktes der FE-Analyse:
analysis_currentDateTime = dt.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
modelName = 'internal_pressure_' + '_mR' + str(mR) +'_mT'+ str(mT) + '_'+ analysis_currentDateTime
#modelName = '0_90_2_R4h_L4h_T1_AO90_AE45_CFK_PS'
# Erstellen eines neuen Ordners in dem angebeben Pfad und Definition dessen als Arbeitsverzeichnis:
analysis_newPath = analysis_Path + '\\' + modelName
os.makedirs(analysis_newPath)
os.chdir(analysis_newPath)
#------------------------------------------------

#------------------------------------------------
# Transformieren in  0 90 Schreibweise
#plyAngle = [90 if char == '1' else 0 for char in ric_plyAngle]
plyAngle = [float(a) for a in ric_plyAngle if -180 <= float(a) <= 180]
N = len(plyAngle)
# Auswertungsinterfaces (Interfaces zweier aufeinanderfolgender Schichten mit unterschiedlichen Faserorientierungswinkeln):
iInterfaceEval = list(range(1,N))
# Startwinkel des Modells, Winkelbereich des Sub-Modells, Oeffnungswinkel,
angleStart,angleOpening = 0.0,np.pi/2

# Materialparameter der 0-Schicht (1 - Radialkoordinate, 2 - Tangentialkoordinate, 3 - Axialkoordinate):
if compositeMaterialName == 'cfk':
	E1, E2, E3 = 7460.0, 118148.0, 7460.0
	Nu12, Nu13, Nu23 = 0.021, 0.37, 0.34
	G12, G13, G23 = 4800.0, 2701.0, 4800.0
	alpha11,alpha22,alpha33=2.6e-5,-1.0e-6,2.6e-5
	# dsingle = 1.25
	dsingle = 2
elif compositeMaterialName == 'gfk':
	E1, E2, E3 = 9552.6, 39296.0, 9552.6
	Nu21, Nu13, Nu23 = 0.29, 0.38, 0.29
	Nu12=E1/E2*Nu21
	G12, G13, G23 = 3080.5, 3449.0, 3080.5
	alpha11,alpha22,alpha33=2.6e-5,8.6e-6,2.6e-5
	dsingle = 0.190
elif compositeMaterialName == 'cfkDuro':
	E1, E2, E3 = 6895.0, 172375.0, 6895.0
	Nu21, Nu13, Nu23 = 0.25, 0.25, 0.25
	Nu12=E1/E2*Nu21
	G12,G13,G23=3448.0,1379.0,3448.0
	alpha11,alpha22,alpha33=2.6e-5,-1.0e-6,2.6e-5
	dsingle = 1.25
else:
	pass

#--------------------------------------------------------------
#Berechnung verschiedener Punkte

#Gesamtdicke des Pruefkörpers
d=dsingle*N

#Außenradius
ra=ri+d

#Mittelpunkt der Kurvenradien
centerpoint=(ra,ra)
rc = np.sqrt(2*ra**2)
#Äußerer Halbkreis Punkt oben x y
outer_upper_arc = [0,ra]

#Äußerer Halbkreis Punkt unten
outer_lower_arc = [ra,0]

#Innerer Halbkreis punkt oben
inner_upper_arc = [d,ra]

#innerer Halbkreis Punkt unten
inner_lower_arc = [ra,d]

#Radien der Mittelflächen in der Kruemmung
rk =[ri]
for ii in range (N):
	rk.append (rk[ii]+dsingle)

print("rk array:", rk)
print("len(rk):", len(rk))
print("N:", N)

rN = rk[-1]

# Betrag des angreifenden Momentenflusses:
bendingMoment = 0.0
if bendingMoment != 0.0:
	boolBendingMoment = True
else:
	boolBendingMoment = False

# Betrag der angreifenden Radiallast:
radialForce = 0.0
if radialForce != 0.0:
	boolRadialForce = True
else:
	boolRadialForce = False

# Betrag des angreifenden Normalkraftflusses:
circForce = 0.0
if circForce != 0.0:
	boolCircForce = True
else:
	boolCircForce = False

# Betrag des angreifenden Innen-/Aussendruckes:
OuterPressure, InnerPressure = 0.0, 1.0
if OuterPressure != 0.0 or InnerPressure != 0.0:
	boolPressure = True
else:
	boolPressure = False

# Betrag des angreifenden zylindrischen Innen-/Aussendruckes:
cylindricalOuterPressure, cylindricalInnerPressure = 0.0, 0.0
if cylindricalOuterPressure != 0.0 or cylindricalInnerPressure != 0.0:
	boolCylindricalPressure = True
else:
	boolCylindricalPressure = False

# Temperaturdifferenz:
tempDif = 0.0
if tempDif != 0.0:
	boolTempDif = True
else:
	boolTempDif = False

# Feuchtekonzentration:
moistDif = 0.0
if moistDif != 0.0:
	boolMoistDif = True
else:
	boolMoistDif = False

if (any([boolBendingMoment, boolRadialForce, boolCircForce]) and any([boolPressure, boolCylindricalPressure])) or (any([boolBendingMoment, boolRadialForce, boolCircForce]) and any([boolTempDif, boolMoistDif])) or (any([boolPressure, boolCylindricalPressure]) and any([boolTempDif, boolMoistDif])):
	raise ValueError('Please check the modelled structural situation!')
else:
	pass

#--------------------------------------------------------------
#Parameter
def radialGeometryParameters():
	"""
	Berechnet und gibt verschiedene radiale Geometrieparameter für Laminatschichten zurück.

	Returns:
		tuple: Ein Tupel, das Folgendes enthält:
			- rm (list): Radien der Mittelflächen der Laminat-Einzelschichten.
			- rEdgesEval (list): Radien aller zu evaluierenden Interfaces.
			- rEdgesAll (list): Radien aller Partitionierungsinterfaces.
			- rFacesAll (list): Mittlere Radien aller partitionierten physikalischen Schichten.
	"""
	# Radien der Mittelflächen der Laminat-Einzelschichten:
	rm = [(rk[ii] + rk[ii + 1]) / 2 for ii in range(N)]
	 
	# Radien aller zu evaluierenden Interfaces:
	print("Creating rEdgesEval...")
	print("iInterfaceEval:", iInterfaceEval)
	print("rk length:", len(rk))
	print("rk:", rk)
	
	try:
		rEdgesEval = [rk[iInterEval] for iInterEval in iInterfaceEval]
		print("rEdgesEval created successfully:", rEdgesEval)
	except IndexError as e:
		print("ERROR creating rEdgesEval:", e)
		print("Trying to access indices:", iInterfaceEval, "in rk of length", len(rk))
		rEdgesEval = []
	 
	# Radien aller Interfaces:
	try:
		rInterfaceEval = [rk[iInterEval] for iInterEval in iInterfaceEval]
		print("rInterfaceEval created successfully:", rInterfaceEval)
	except IndexError as e:
		print("ERROR creating rInterfaceEval:", e)
		print("Trying to access indices:", iInterfaceEval, "in rk of length", len(rk))
		rInterfaceEval = []
	 
	# Radien aller Partitionierungsinterfaces:
	rEdgesAll = []
	for ii in range(N + 1):
		rEdgesAll.append(rk[ii])
	 
	# Mittlere Radien aller partitionierten physikalischen Schichten:
	rFacesAll = [(rEdgesAll[ii] + rEdgesAll[ii + 1]) / 2 for ii in range(N)]
	 
	return (rm, rEdgesEval, rEdgesAll, rInterfaceEval, rFacesAll)

def circGeometryParameters():
	"""
	Berechnet und gibt verschiedene Umfangswinkel-Geometrieparameter zurück.

	Returns:
		tuple: Ein Tupel, das Folgendes enthält:
			- thetaPartition (list): Umfangswinkel-Partitionierungen.
			- thetaEdges (list): Umfangswinkel-Kanten.
			- thetaFaces (list): Umfangswinkel-Mittenflächen.
			- thetaPartitionFaces (list): Mittlere Umfangswinkel-Partitionierungen.
	"""
	# Umfangswinkel-Kanten:
	thetaEdges = [angleStart, angleEnd]#, angleOpening]
	 
	# Umfangswinkel-Partitionierungen:
	thetaPartition = [angleEnd]
	 
	# Umfangswinkel-Mittenflächen:
	refinementInterval = np.pi * 10 / 180
	thetaPartiotionsAll = [
		angleStart,
		angleEnd - refinementInterval,
		angleEnd]#,
#        angleEnd + refinementInterval,
#        angleOpening
#    ]
	thetaPartitionFaces = [(thetaPartiotionsAll[ii] + thetaPartiotionsAll[ii + 1]) / 2 for ii in range(len(thetaPartiotionsAll) - 1)]
	thetaFaces = [(thetaEdges[ii] + thetaEdges[ii + 1]) / 2 for ii in range(len(thetaEdges) - 1)]
	 
	return (thetaPartition, thetaEdges, thetaFaces, thetaPartitionFaces)

rm,rEdgesEval,rEdgesAll,rInterfaceEval,rFacesAll = radialGeometryParameters()
thetaPartition,thetaEdges,thetaFaces,thetaPartitionFaces = circGeometryParameters()

#--------------------------------------------------------------
#Dateien für Auswärtung
Ratio = (ri + d) / ri 

#--------------------------------------------------------------
#SKIZZE ERSTELLEN
Mdb()
# Name ändern
mdb.models.changeKey(fromName='Model-1', toName= modelName)
mdb.models[modelName].rootAssembly.clearGeometryCache()
mdb.models[modelName].rootAssembly.regenerate()
LModel = mdb.models[modelName]

def create_closed_L_shape_sketch_internal_pressure(LModel, rk, ltan, thetaEdges, sketch_name='L-Shape-Closed-Internal-Pressure'):
    """
    Creates a closed L-shaped sketch for internal pressure model using proper circular arcs and connected path.
    
    This function creates the L-shaped geometry as ONE CONTINUOUS CLOSED PATH,
    which solves the "Shell extrude feature failed" error by ensuring BaseShell
    can recognize a proper closed region.
    
    This version is adapted for the internal pressure geometry which uses TWO_D_PLANAR
    parts and has a specific L-bracket pattern (curved quarter-circle + downward flank only).
    
    Args:
        LModel: Abaqus model object
        rk (list): Interface radii [inner_radius, ..., outer_radius]
        ltan (float): Length of the tangent (flank height)
        thetaEdges (list): Theta angles for the arcs [0, pi/2]
        sketch_name (str): Name for the sketch
        
    Returns:
        tuple: (LModelSketch, OuterCirc_id, LowerOuterTan_id)
    """
    
    print("Creating closed L-shaped sketch for internal pressure model...")
    
    # Get geometry parameters
    r_inner = rk[0]   # Inner radius
    r_outer = rk[-1]  # Outer radius  
    height = ltan     # Flank height
    
    print("Geometry parameters: r_inner={}, r_outer={}, height={}".format(
        r_inner, r_outer, height))
    
    # Create the sketch - using sheet size from original
    LModelSketch = LModel.ConstrainedSketch(name=sketch_name, sheetSize=(2 * (ltan + ra)))
    
    # Add construction line for axis of symmetry (Y-axis for axisymmetric analysis)
    print("Adding construction line for Y-axis (axis of symmetry)...")
    LModelSketch.ConstructionLine(point1=(0.0, -height-10.0), point2=(0.0, r_outer+10.0))
    print("SUCCESS: Construction line added for axis of symmetry")
    
    # Create the CLOSED L-shaped path for internal pressure geometry
    print("Creating complete closed L-shaped profile for internal pressure geometry...")
    
    # This internal pressure model creates a simpler L-shape:
    # - Quarter-circle curved section from 0° to 90° (0 to pi/2)
    # - Downward vertical flank from curved section
    # - No upper horizontal extension (unlike the angleply model)
    
    # THE CRITICAL INSIGHT: Create ONE CONTINUOUS CLOSED PATH
    # Start from bottom-left and go clockwise to ensure proper closure
    
    # 1. Bottom horizontal line (left to right)
    print("1. Creating bottom horizontal line")
    LModelSketch.Line(point1=(r_inner, -height), point2=(r_outer, -height))
    
    # 2. Right vertical line (bottom to curved section start)  
    print("2. Creating right vertical line")
    LowerOuterTan = LModelSketch.Line(point1=(r_outer, -height), point2=(r_outer, 0.0))
    LowerOuterTan_id = LowerOuterTan.id
    
    # 3. Outer arc: quarter-circle from (r_outer, 0) to (0, r_outer)
    outer_start = (r_outer, 0.0)
    outer_end = (ct.pol2cart_x(r_outer, thetaEdges[1]), ct.pol2cart_y(r_outer, thetaEdges[1]))
    
    print("3. Creating outer arc: start={}, end={}".format(outer_start, outer_end))
    OuterCirc = LModelSketch.ArcByCenterEnds(
        center=(0.0, 0.0), 
        point1=outer_start, 
        point2=outer_end, 
        direction=COUNTERCLOCKWISE)
    OuterCirc_id = OuterCirc.id
    
    # 4. Top radial line (connects arc ends at Y-axis, outer to inner)
    print("4. Creating top radial line")
    LModelSketch.Line(
        point1=(ct.pol2cart_x(r_outer, thetaEdges[1]), ct.pol2cart_y(r_outer, thetaEdges[1])), 
        point2=(ct.pol2cart_x(r_inner, thetaEdges[1]), ct.pol2cart_y(r_inner, thetaEdges[1])))
    
    # 5. Inner arc: quarter-circle from (0, r_inner) to (r_inner, 0) - CLOCKWISE to close properly
    inner_start = (ct.pol2cart_x(r_inner, thetaEdges[1]), ct.pol2cart_y(r_inner, thetaEdges[1]))
    inner_end = (r_inner, 0.0)
    
    print("5. Creating inner arc: start={}, end={}".format(inner_start, inner_end))
    LModelSketch.ArcByCenterEnds(
        center=(0.0, 0.0), 
        point1=inner_start, 
        point2=inner_end, 
        direction=CLOCKWISE)
    
    # 6. Left vertical line (middle to bottom) - closes the L-shape
    print("6. Creating left vertical line (closes the shape)")
    LModelSketch.Line(point1=(r_inner, 0.0), point2=(r_inner, -height))
    
    print("SUCCESS: Created CLOSED L-shape for internal pressure model")
    print("  - 2 quarter-circle arcs (mathematically perfect)")
    print("  - 4 straight connecting lines")
    print("  - Total: 6 geometric entities forming ONE CLOSED REGION")
    print("  - Path: bottom→right→outer_arc→top→inner_arc→left→close")
    
    return LModelSketch, OuterCirc_id, LowerOuterTan_id

def sketchLModel():
	"""
	Erstellt die Skizze des L-Profils und die IDs der äußeren Elemente zurück.
	
	UPDATED: Now uses the working closed-path approach to avoid "Shell extrude feature failed" error.

	Returns:
		tuple:
			- LModelPart: Die erstellte Part-Instanz des L-Profils.
			- OuterCirc_id: Die ID des äußeren Halbkreises.
			- LowerOuterTan_id: Die ID der unteren äußeren Tangente.
	"""
	
	# Use the working closed-path sketch approach
	LModelSketch, OuterCirc_id, LowerOuterTan_id = create_closed_L_shape_sketch_internal_pressure(
		LModel, rk, ltan, thetaEdges, 'L-Profil-Sketch-Working')
	 
	# Create the part using TWO_D_PLANAR (correct for both plane strain and axisymmetric)
	LModelPart = LModel.Part(name='L-Profil-Part', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
	
	try:
		LModelPart.BaseShell(sketch=LModelSketch)
		print("SUCCESS: BaseShell operation completed successfully with closed-path approach!")
	except Exception as e:
		print("BaseShell failed even with closed-path approach: {}".format(str(e)))
		print("This suggests a more fundamental issue with the geometry or Abaqus setup.")
		raise e
	 
	return (LModelPart, OuterCirc_id, LowerOuterTan_id)#, UpperOuterTan_id)

#LModelPart,OuterCirc_id,LowerOuterTan_id,UpperOuterTan_id=sketchLModel()
LModelPart,OuterCirc_id,LowerOuterTan_id=sketchLModel()

#--------------------------------------------------------------
#Abkürzungen als Variablen speichern
LModelAssembly = LModel.rootAssembly
LModelPart = mdb.models[modelName].parts['L-Profil-Part']
LModelFaces = LModelPart.faces

def instanceLModel():
	"""
	Erstellt Instanzen und positioniert sie entsprechend.

	Returns:
		Instance: Die letzte erstellte Instanz des starren Körpers.
	"""
	LModelinstance = LModelAssembly.Instance(name='L_Profil_Instance', part=LModelPart, dependent=ON)
	
	return (LModelinstance)



LModelInstance = instanceLModel()
#--------------------------------------------------------------
#Partionierung

def partitionComposite():
    """
    Creates partitions in the L-profile based on ply angles and layer thickness.
    """
    LModelPartition = mdb.models[modelName].parts['L-Profil-Part']
    partitionfaces = LModelPartition.faces
    # Create sketch transformation for partitioning
    LModelTransform = LModelPartition.MakeSketchTransform(
        sketchPlane=partitionfaces[0],
        sketchPlaneSide=SIDE1,
        origin=(0, 0, 0))
    # Create sketch for partitioning
    LModelPartitionSketch = mdb.models[modelName].ConstrainedSketch(
        name='L-Profil-Partition-Sketch',
        sheetSize=(2 * (ltan + ra)),
        transform=LModelTransform)
    g = LModelPartitionSketch.geometry
    # Project references onto the sketch
    LModelPartition.projectReferencesOntoSketch(sketch=LModelPartitionSketch, filter=COPLANAR_EDGES)
    # Create offsets and partitions for each layer
    for ii in range(N - 1):
        LModelPartitionSketch.offset(
            distance=((ii + 1) * dsingle),
            objectList=(
                g.findAt((ct.pol2cart_x(rk[0], thetaFaces[0]), ct.pol2cart_y(rk[0], thetaFaces[0]))),
                g.findAt((ct.pol2cart_x(rk[0], thetaEdges[0]), 0 - ltan / 2))),
            side=RIGHT)
        LModelPickedFaces = partitionfaces.findAt((ct.pol2cart_x(rm[ii], thetaFaces[0]), ct.pol2cart_y(rm[ii], thetaFaces[0]), 0),)
        LModelPartition.PartitionFaceBySketch(faces=LModelPickedFaces, sketch=LModelPartitionSketch)
    
    # Add simple theta partitions for edge set limits
    angleLimit = 30.0 * np.pi / 180.0  # 30 degrees limit for curved section
    if Null_Auswaertung:
        angleEval_partition = 0
    else:
        angleEval_partition = angleEval
    
    # Create theta partitions at the limits for curved edges
    thetaPartitionStart = max(0.0, angleEval_partition - angleLimit)
    thetaPartitionEnd = min(angleEnd, angleEval_partition + angleLimit)
    
    # Add these angles to the existing thetaPartition list for partitioning
    try:
        # Create datum points for theta partitioning
        for theta in [thetaPartitionStart, thetaPartitionEnd]:
            if theta != 0.0 and theta != angleEnd:  # Don't create partitions at existing boundaries
                LModelPartition.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[0], theta), ct.pol2cart_y(rk[0], theta), 0.0))
                LModelPartition.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], theta), ct.pol2cart_y(rk[-1], theta), 0.0))
                LModelPartition.DatumPointByCoordinate(coords=(ct.pol2cart_x(rk[-1], theta), ct.pol2cart_y(rk[-1], theta), length))
                
                # Get the last three datum points created
                datumKeys = list(LModelPartition.datums.keys())
                datumKeys.sort()
                
                # Create datum plane through these three points for partitioning
                LModelPartition.DatumPlaneByThreePoints(
                    point1=LModelPartition.datums[datumKeys[-3]],
                    point2=LModelPartition.datums[datumKeys[-2]], 
                    point3=LModelPartition.datums[datumKeys[-1]])
                
                # Get the datum plane ID and partition faces
                planeKeys = [key for key in LModelPartition.features.keys() if key.startswith('Datum plane-')]
                if planeKeys:
                    planeKeys.sort()
                    latestPlaneKey = planeKeys[-1]
                    planeId = LModelPartition.features[latestPlaneKey].id
                    LModelPartition.PartitionFaceByDatumPlane(faces=partitionfaces, datumPlane=LModelPartition.datums[planeId])
    except:
        print("Warning: Could not create theta partitions for edge limits")

partitionComposite()
#--------------------------------------------------------------
#materialRigid()

#--------------------------------------------------------------
#Koordinatensystem festlegen
# For axisymmetric model: Flank section uses CYLINDRICAL coordinate system
# r in (1,0,0), theta in (0,0,1), z is axial
LModelPart.DatumCsysByThreePoints(name='Cartesian-KOS', coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
cartesianCompositeCSCylId = LModelPart.features['Cartesian-KOS'].id

# For axisymmetric model: Curved section uses SPHERICAL coordinate system
# r in (1,0,0), phi (azimuthal) in (0,0,1), theta (polar) defined by geometry
#LModelPart.DatumCsysByThreePoints(origin=(ra,ra,0), point1=LModelPart.vertices.findAt(coordinates=(0, ra, 0.0)), point2=LModelPart.vertices.findAt(coordinates=(ra, 0.0, 0.0)), name='Cylindrical-KOS', coordSysType=SPHERICAL)
#LModelPart.DatumCsysByThreePoints(origin=(0,0,0), point1=LModelPart.vertices.findAt(coordinates=(0, ra, 0.0)), point2=LModelPart.vertices.findAt(coordinates=(ra, 0.0, 0.0)), name='Cylindrical-KOS', coordSysType=SPHERICAL)
LModelPart.DatumCsysByThreePoints(origin=(rk[-1],rk[-1],0), point1=(rk[-1], -ltan, 0.0), point2=(-ltan, rk[-1], 0.0), name='Cylindrical-KOS', coordSysType=SPHERICAL)
#LModelPart.DatumCsysByThreePoints(name='Cylindrical-KOS', coordSysType=CARTESIAN, origin=(ra, ra, 0.0), point1=(0.0, 0.0, 0.0), line2=(-ltan, ltan , 0.0))
#LModelPart.DatumCsysByThreePoints(name='Cylindrical-KOS', coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
curvedCompositeCSCylId = LModelPart.features['Cylindrical-KOS'].id
curvedCompositeCylCoordSys = LModel.rootAssembly.instances['L_Profil_Instance'].datums[curvedCompositeCSCylId]

# Spherical coordinate system at origin for orientation
# r in (1,0,0), phi in (0,0,1)
LModelPart.DatumCsysByThreePoints(origin=(0,0,0), point1=(1, 0 , 0.0), point2=(0, 1.0, 0.0), name='Cylindrical-KOS-Orientation', coordSysType=SPHERICAL)
orientationKOSid = LModelPart.features['Cylindrical-KOS-Orientation'].id


#KOS Für Boundary Contditions
LModelAssembly.DatumCsysByThreePoints(name='BoundaryCD-Datum', coordSysType=CARTESIAN, origin=(ra, ra, 0.0), 
									  point1=(0.0, 0.0, 0.0), line2=(-ltan, ltan , 0.0))
BoundaryCDDatum_ID = LModelAssembly.features['BoundaryCD-Datum'].id
BoundaryCDDatum = LModelAssembly.datums[BoundaryCDDatum_ID]

#--------------------------------------------------------------
#Partitionierung der Flanken für besseres Mesh
def partitionflanken():
	"""
	Partitioniert die Flanken des L-Profils für eine bessere Vernetzung.
	"""
	try:
		LModelPart = mdb.models[modelName].parts['L-Profil-Part']
		LModelfaces = LModelPart.faces
		
		# Create sketch transformation for partitioning
		LModelTransform = LModelPart.MakeSketchTransform(
			sketchPlane=LModelfaces.findAt((ct.pol2cart_x(rm[0], thetaFaces[0]), ct.pol2cart_y(rm[0], thetaFaces[0]), 0.0), normal=(0.0, 0.0, 1.0)), 
			sketchPlaneSide=SIDE1, 
			origin=(0, 0, 0.0))
		
		LModelSketch = mdb.models[modelName].ConstrainedSketch(name='FlankenPartition', sheetSize=200, transform=LModelTransform)
		
		# Project references onto the sketch first
		LModelPart.projectReferencesOntoSketch(sketch=LModelSketch, filter=COPLANAR_EDGES)
		
		# Partitionierung der Flanken - simplified approach
		# Create radial lines for partitioning at theta edges
		for jj in range(len(thetaEdges)):
			if jj < len(thetaEdges):  # Safety check
				try:
					# Create radial line from inner to outer radius at each theta edge
					LModelSketch.Line(
						point1=(ct.pol2cart_x(rk[0], thetaEdges[jj]), ct.pol2cart_y(rk[0], thetaEdges[jj])), 
						point2=(ct.pol2cart_x(rk[-1], thetaEdges[jj]), ct.pol2cart_y(rk[-1], thetaEdges[jj])))
				except:
					print("Warning: Could not create partition line for theta edge", jj)
					continue
		
		# Add horizontal partition line for right edge limit
		# Create horizontal line at certain percentage of -ltan (20% for consistency with setsCurvedCompositeCirc)
		flankLimitPercentage = 0.15  # 20% of flank length from corner
		yPartitionLimit = -ltan * flankLimitPercentage
		
		try:
			# Create horizontal line from inner radius to outer radius at the partition y-coordinate
			LModelSketch.Line(
				point1=(rk[0], yPartitionLimit), 
				point2=(rk[-1], yPartitionLimit))
			print("Created horizontal partition line at y =", yPartitionLimit)
		except:
			print("Warning: Could not create horizontal partition line for right edge limit")
		
		# Apply partitioning to all faces once
		try:
			LModelPart.PartitionFaceBySketch(faces=LModelfaces, sketch=LModelSketch)
		except:
			print("Warning: Could not partition faces with flank sketch")
			
	except Exception as e:
		print("Warning: partitionflanken failed:", str(e))
		print("Continuing without flank partitioning...")

partitionflanken()

#--------------------------------------------------------------

#--------------------------------------------------------------
#Meshgenerierung
flankenFaces = []
for ii in range(1,6,2):
	flankenFaces.append(-float(ltan)*ii/6)

def meshLPart(mR, mT, mFlanken, mBias, mRRatio, mTRatio):
	"""
	Erstellt ein Netz für das L-Profil unter Verwendung der angegebenen Parameter.
	 
	Args:
		mR (int): Anzahl der Elemente in radialer Richtung.
		mT (int): Anzahl der Elemente in tangentialer Richtung.
		mFlanken (int): Anzahl der Elemente entlang der Flanken.
		mBias (float): Bias-Wert für das Netz.
		mRRatio (float): Verhältnis für das Bias im radialen Netz.
		mTRatio (float): Verhältnis für das Bias im tangentialen Netz.
	 
	Raises:
		ValueError: Wenn eine Kante keine Richtung hat.
	"""
	for ii in range(len(rFacesAll)):
		for jj in range(len(thetaEdges)):
			curvedCompositeMeshEdge = LModelPart.edges.findAt(
				((ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj]),
				  ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj]), 0.0),)
			)
			LModelPart.seedEdgeByNumber(
				edges=curvedCompositeMeshEdge,
				number=mR,
				constraint=FINER
			)
	 
	for ii in range(len(rk)):
		  curvedCompositeMeshEdge = LModelPart.edges.findAt(
				((ct.pol2cart_x(rk[ii], thetaFaces[-1]),
				  ct.pol2cart_y(rk[ii], thetaFaces[-1]), 0.0),)
			)
		  LModelPart.seedEdgeByNumber(
			edges=curvedCompositeMeshEdge,
			number=mT*2,
			constraint=FINER
			) 
	 
	for jj in range(len(flankenFaces)):
		for ii in range(len(rk)):
			CurvedCompositeMeshEdge = LModelPart.edges.findAt(
				((rk[ii], flankenFaces[jj], 0.0),)
			)
			LModelPart.seedEdgeByNumber(
				edges=CurvedCompositeMeshEdge,
				number=mFlanken,
				constraint=FINER
			)
	 
	# for jj in range(len(flankenEdges)):
	# 	for ii in range(len(rFacesAll)):
	# 		curvedCompositeMeshEdge = LModelPart.edges.findAt(
	# 			((flankenEdges[jj],
	# 			  ct.pol2cart_y(rm[ii], thetaEdges[-1]), 0.0),)
	# 		)
	# 		LModelPart.seedEdgeByBias(
	# 			biasMethod=DOUBLE,
	# 			endEdges=curvedCompositeMeshEdge,
	# 			ratio=mRRatio,
	# 			number=mR,
	# 			constraint=FINER
	# 		)
	compositeElementType = mesh.ElemType(elemCode=CAX8RT, elemLibrary=STANDARD)
	LModelPart.setMeshControls(
		regions=LModelPart.faces,
		elemShape=QUAD,
		technique=STRUCTURED
	)
	LModelPart.setElementType(
		regions=(LModelPart.faces[:],),
		elemTypes=(compositeElementType,)
	)
	LModelPart.seedPart(size=0.4, deviationFactor=0.1, constraint=FINER)
	LModelPart.generateMesh()

meshLPart(mR, mT, mFlanken, mBias, mRRatio, mTRatio)


#--------------------------------------------------------------
#Sets zur Auswärtung

if Null_Auswaertung == True:
	angleEval_FEInnner = 0
else:
	angleEval_FEInnner =angleEval

def setsCurvedCompositeCirc(rFacesAll, angleEval_FEInnner, iInterfaceEval, thetaPartitionFaces, angleEnd, rk):
	"""
	Erstellt Sets zur Verifizierung der Inner-Solution und zur Bestimmung des Abklingverhaltens der singulären Spannungen.
	Nutzt die in partitionComposite() erstellten Theta-Partitionen für begrenzte Edge-Sets.

	Args:
		rFacesAll (list): Liste aller Radien der Flächen.
		angleEval_FEInnner (float): Winkel zur Auswertung der inneren Lösung.
		iInterfaceEval (list): Liste der auszuwertenden Interfaces.
		thetaPartitionFaces (list): Liste der Winkelpartitionen.
		angleEnd (float): Geometrie-Endwinkel.
		rk (list): Liste der Radien.
	"""
	# Calculate the same limits as used in partitionComposite
	# Sets zur Verifizierung der Inner-Solution:
	for ii in range(len(rFacesAll)):
		try:
			compositeSetEvalEdge = LModelInstance.edges.findAt(((ct.pol2cart_x(rFacesAll[ii], angleEval_FEInnner), ct.pol2cart_y(rFacesAll[ii], angleEval_FEInnner), 0.0),))
			LModelAssembly.Set(edges=compositeSetEvalEdge, name='FEInner' + str(ii + 1))
		except:
			print("Warning: Could not create FEInner set for layer", ii + 1)
	
	# Sets für Interface-Analyse mit begrenzten Theta-Bereichen:
	# Right side sets - straight flank (no partitioning needed, just use a point along the flank)
	for kk in range(len(iInterfaceEval)):
		ii = iInterfaceEval[kk]
		try:
			# Right side: straight flank at constant radius, limited to upper part near corner
			curvedCompositeEdgesRight = LModelInstance.edges.findAt(((rk[ii], -ltan * 0.1, 0.0),))
			LModelAssembly.Set(edges=curvedCompositeEdgesRight, name='FEInterfaceCircRight' + str(iInterfaceEval[kk]))
		except:
			print("Warning: Could not find right flank edge for interface", iInterfaceEval[kk])
	
	# Left side sets - curved section using partitioned edges
	for kk in range(len(iInterfaceEval)):
		ii = iInterfaceEval[kk]
		try:
			# Left side: curved section, use partitioned edge in the middle of the limited angular range
			midAngle = angleEval_FEInnner  # Start from evaluation angle
			curvedCompositeEdgesLeft = LModelInstance.edges.findAt(((ct.pol2cart_x(rk[ii], midAngle + 0.00001), ct.pol2cart_y(rk[ii], midAngle + 0.00001), 0.0),))
			LModelAssembly.Set(edges=curvedCompositeEdgesLeft, name='FEInterfaceCircLeft' + str(iInterfaceEval[kk]))
		except:
			print("Warning: Could not find partitioned left edge for interface", iInterfaceEval[kk])

setsCurvedCompositeCirc(rFacesAll, angleEval_FEInnner, iInterfaceEval, thetaPartitionFaces, angleEnd, rk)

#### Material properties ####
mdb.models[modelName].Material(name='Material-1')
mdb.models[modelName].materials['Material-1'].Elastic(
    temperatureDependency=ON, 
    table=((85700.0, 0.318, 0.0),
           (84500.0, 0.32, 100.0),
           (81200.0, 0.325, 200.0),
           (77400.0, 0.33, 300.0))
)
mdb.models[modelName].materials['Material-1'].Conductivity(
    temperatureDependency=ON, 
    table=((0.0306, 20.0),
           (0.056, 73.0),
           (0.077, 123.0),
           (0.107, 223.0),
           (0.123, 300.0))
)
mdb.models[modelName].materials['Material-1'].Expansion(
    table=((1.44e-05, 20.0),
           (1.72e-05, 73.0),
           (1.94e-05, 123.0),
           (2.13e-05, 223.0),
           (2.28e-05, 300.0)),
    zero=293.0,
    temperatureDependency=ON
)

mdb.models[modelName].HomogeneousSolidSection(
    name='Section-1', 
    material='Material-1', 
    thickness=None
)

# --- SECTION ASSIGNMENT USING BOUNDING BOX ---
p = mdb.models[modelName].parts['L-Profil-Part']

# Select the entire geometry using bounding box
faces = p.faces.getByBoundingBox(
    xMin=-1e6, yMin=-1e6, zMin=-1e6,
    xMax= 1e6, yMax= 1e6, zMax= 1e6
)

region = regionToolset.Region(faces=faces)

p.SectionAssignment(
    region=region,
    sectionName='Section-1',
    offset=0.0,
    offsetType=MIDDLE_SURFACE,
    offsetField='',
    thicknessAssignment=FROM_SECTION
)


#--------------------------------------------------------------
#def stepInitialise():
	#"""
	#Erstellt einen quasi-statischen Schritt.

	#Returns:
	#    str: Der Name des erstellten Schritts.
	#"""
	#step = 'BendingMoment'
	#LModel.StaticStep(name=step, previous='Initial', description='Introduce bending moment', nlgeom=ON, initialInc=0.1)
	#return step

#step = stepInitialise()



def stepMechanicalLoads():
	stepPrevious,step,stepDescription = 'Initial','Mechanical_Loading','Introduce mechanical loadings'
	LModel.StaticStep(name=step, previous=stepPrevious, description=stepDescription, nlgeom=OFF, initialInc=1.0)
	LModel.CoupledTempDisplacementStep(name=step, previous=stepPrevious, maintainAttributes=True, 
    response=STEADY_STATE, deltmx=None, cetol=None, creepIntegration=None, 
    amplitude=RAMP)
	return (stepPrevious,step)

#bendingMoment = -1.0
stepPrevious,step = stepMechanicalLoads()

#---------------------------------------------------------------
# Randbedingungen:
def boundaryConditionBendingMoment():
	pass

def boundaryConditionEdgeLoad():
	# Sperren der radialen Bewegung:
	curvedCompositeFixRadialEdge = LModelInstance.vertices.findAt(((ct.pol2cart_x(ri,  thetaEdges[-1]), ct.pol2cart_y(ri, thetaEdges[-1]), 0.0),))
	curvedCompositeFixRadialSet = LModelAssembly.Set(vertices=curvedCompositeFixRadialEdge, name='Set_Fix_U_Interface_'+ str(0))
	LModel.DisplacementBC(name='Fix_U_Interface_'+ str(0), createStepName=stepPrevious, region=curvedCompositeFixRadialSet, u1=0.0)
	# Sperren der tangentialen Bewegung
	for ii in range(len(rFacesAll)):
		curvedCompositeStringFixCircFaces = ''
		curvedCompositeStringFixCircFaces = curvedCompositeStringFixCircFaces + str(((ct.pol2cart_x(rFacesAll[ii], thetaEdges[-1]), ct.pol2cart_y(rFacesAll[ii], thetaEdges[-1]), 0.0),),) + ','
		curvedCompositeFixCircFacesExec = 'curvedCompositeFixCircFaces = LModelInstance.edges.findAt(' + curvedCompositeStringFixCircFaces + ')'
		exec(curvedCompositeFixCircFacesExec)
		
		curvedCompositeFixCircSet = LModelAssembly.Set(edges=curvedCompositeFixCircFaces, name='Set_Fix_V_Layer_'+ str(ii+1))
		LModel.DisplacementBC(name='Fix_V_Layer_'+ str(ii+1), createStepName=stepPrevious, region=curvedCompositeFixCircSet, u2 = 0.0)

def boundaryConditionSurfaceLoad():
	for jj in range(1,3):#enumerate([thetaEdges[0],thetaEdges[-1]]):
		# Build list of edge coordinates
		edge_coords = []
		for ii in range(len(rFacesAll)):
			coord = (ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj-1]), 
			         ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj-1]) + (jj-2)*ltan, 
			         0.0)
			edge_coords.append((coord,))
		
		# Find edges at these coordinates
		curvedCompositeFixCircFaces = LModelInstance.edges.findAt(*edge_coords)
		curvedCompositeFixCircSet = LModelAssembly.Set(edges=curvedCompositeFixCircFaces, name='Set_Fix_U_SS_'+ str(jj))
		
		if jj == 2:
			LModel.DisplacementBC(name='Fix_U_SS_'+ str(jj), createStepName=stepPrevious, region=curvedCompositeFixCircSet, u1=0.0)
		else:
			LModel.DisplacementBC(name='Fix_U_SS_'+ str(jj), createStepName=stepPrevious, region=curvedCompositeFixCircSet, u2=0.0)#, u1=0.0)

if boolBendingMoment and not any([boolRadialForce,boolCircForce]):
	boundaryConditionBendingMoment()
elif any([boolBendingMoment,boolRadialForce,boolCircForce]):
	boundaryConditionEdgeLoad()
elif any([boolPressure,boolCylindricalPressure]):
	boundaryConditionSurfaceLoad()
else:
	pass

def loadsCreateKinCoup(jj):
	# Erstellen eines Referenzpunktes und Fixierung dessen:
	if jj == 1:
		LModelAssembly.ReferencePoint(point=((
			ct.pol2cart_x(rk[N//2], thetaEdges[jj-1]) + (jj-1)*ltan, 
			ct.pol2cart_y(rk[N//2], thetaEdges[jj-1]) - (jj)*ltan, 
			0)))
		
		curvedCompositeLoadRefPoint = (
			LModelAssembly.referencePoints.findAt((
				ct.pol2cart_x(rk[N//2], thetaEdges[jj-1]) + (jj-1)*ltan, 
				ct.pol2cart_y(rk[N//2], thetaEdges[jj-1]) - (jj)*ltan, 
				0),),)
	else:
		LModelAssembly.ReferencePoint(point=((
			ct.pol2cart_x(rk[N//2], thetaEdges[jj-1]), 
			ct.pol2cart_y(rk[N//2], thetaEdges[jj-1]), 
			0)))
		
		curvedCompositeLoadRefPoint = (
			LModelAssembly.referencePoints.findAt((
				ct.pol2cart_x(rk[N//2], thetaEdges[jj-1]), 
				ct.pol2cart_y(rk[N//2], thetaEdges[jj-1]), 
				0),),)
	
	curvedCompositeLoadRefPointRegion = regionToolset.Region(referencePoints=curvedCompositeLoadRefPoint)
	
	curvedCompositeFixRefPointSet = LModelAssembly.Set(
		region=curvedCompositeLoadRefPointRegion, 
		name='Set_Fix_RefPoint_Pos_' + str(jj*int(180*angleOpening/np.pi)))
	
	# LModel.DisplacementBC(name='Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeFixRefPointSet, localCsys=curvedCompositeCylCoordSys, u3=0.0, ur1=0.0, ur2=0.0)
	LModel.DisplacementBC(name='Fix_RefPoint_Pos_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeFixRefPointSet, u3=0.0, ur1=0.0, ur2=0.0)
	
	# LModel.DisplacementBC(
		# name='Fix_RefPoint_Pos_' + str(jj*int(180*angleOpening/np.pi)), 
		# createStepName=step, 
		# region=curvedCompositeFixRefPointSet, 
		# localCsys=LModelPart.datums[cartesianCompositeCSCylId], 
		# u3=0.0, ur1=0.0, ur2=0.0)
	
	# Kinematic Coupling:
	curvedCompositeStringKinCoupFaces = ''
	for ii in range(len(rFacesAll)):
		curvedCompositeStringKinCoupFaces = curvedCompositeStringKinCoupFaces + str(((ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj-1]), ct.pol2cart_y(rFacesAll[ii], thetaEdges[jj-1]) + (jj-2)*ltan, 0.0),),) + ','#ct.pol2cart_x(rFacesAll[ii], thetaEdges[jj-1]) + (jj-1)*ltan
	curvedCompositeKinCoupFacesExec = 'curvedCompositeKinCoupFaces = LModelInstance.edges.findAt(' + curvedCompositeStringKinCoupFaces + ')'
	exec(curvedCompositeKinCoupFacesExec)
	curvedCompositeKinCoupRegion = regionToolset.Region(edges=curvedCompositeKinCoupFaces)
	# LModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion,influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=OFF, u2=ON, ur3=ON)
	if jj == 2:
		LModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=ON, u2=ON, ur3=ON)#couplingType=KINEMATIC, u1=ON, u2=OFF, ur3=ON)
	else:
		LModel.Coupling(name='Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi)), controlPoint=curvedCompositeLoadRefPointRegion, surface=curvedCompositeKinCoupRegion, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=ON, u2=ON, ur3=ON)#u1=OFF, u2=ON, ur3=ON)
	# LModel.constraints['Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi))].setValues(localCsys=curvedCompositeCylCoordSys)
	LModel.constraints['Kinematic_Coupling_'+ str(jj*int(180*angleOpening/np.pi))].setValues()
	return curvedCompositeLoadRefPointRegion

def loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion):
	# Aeussere Last - Moment:
	#LModel.Moment(name='Load_BendingMoment_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeLoadRefPointRegion, cm3=((-1)**(jj))*bendingMoment, localCsys=curvedCompositeCylCoordSys, distributionType=UNIFORM)
	LModel.Moment(name='Load_BendingMoment_'+ str(jj*int(180*angleOpening/np.pi)), createStepName=step, region=curvedCompositeLoadRefPointRegion, cm3=((-1)**(jj))*bendingMoment, distributionType=UNIFORM)

def loadsRadialForce():
	# Aeussere Last - Querkraft:
	LModel.ConcentratedForce(name='Load_RadialForce', createStepName=step, region=curvedCompositeLoadRefPointRegion, cf1=radialForce, distributionType=UNIFORM)

def loadsSurfaceLoad():
	if OuterPressure != 0.0 or cylindricalOuterPressure != 0.0:
		# Aeussere Last - Aussendruck:
		curvedCompositeStringFaces = ''
		for jj in range(len(thetaFaces)):
			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[-1], thetaFaces[jj]), ct.pol2cart_y(rk[-1], thetaFaces[jj]), 0.0),),) + ','
		
		curvedCompositeStringFacesExec = 'curvedCompositeOuterPressureFaces = LModelInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
		exec(curvedCompositeStringFacesExec)
		
		if boolCylindricalPressure:
			LModel.ExpressionField(name='Load_OuterPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			LModel.Pressure(name='Load_OuterPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeOuterPressureFaces), magnitude=cylindricalOuterPressure, distributionType=FIELD, field='Load_OuterPressure_CylindricalBending', amplitude=UNSET)
		else:
			LModel.Pressure(name='Load_OuterPressure', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeOuterPressureFaces), magnitude=OuterPressure, distributionType=UNIFORM)
	if InnerPressure != 0.0 or cylindricalInnerPressure != 0.0:
		# Aeussere Last - Innendruck:
		curvedCompositeStringFaces = ''
#		for jj in range(len(thetaFaces)):
#			curvedCompositeStringFaces = curvedCompositeStringFaces + str(((ct.pol2cart_x(rk[0], thetaFaces[jj]), ct.pol2cart_y(rk[0], thetaFaces[jj]), 0.0),),) + ','
#		
#		curvedCompositeStringFacesExec = 'curvedCompositeInnerPressureFaces = LModelInstance.edges.findAt(' + curvedCompositeStringFaces + ')'
#		exec(curvedCompositeStringFacesExec)
#		for jj in range(len(thetaPartitionFaces) / 2):
#			curvedCompositeStringEdgesRight = curvedCompositeStringEdgesRight + str(((ct.pol2cart_x(rk[0], angleEval_FEInnner - 0.00001), ct.pol2cart_y(rk[0], angleEval_FEInnner - 0.00001), 0.0),),) + ','
#			curvedCompositeStringEdgesLeft = curvedCompositeStringEdgesLeft + str(((ct.pol2cart_x(rk[0], angleEval - 0.00001), ct.pol2cart_y(rk[0], angleEval - 0.00001), 0.0),),) + ','
#		curvedCompositeStringFacesRightExec = 'curvedCompositeInnerPressureFaces = LModelInstance.edges.findAt(' + curvedCompositeStringEdgesRight + ')'
#		curvedCompositeStringFacesLeftExec = 'curvedCompositeInnerPressureFaces = LModelInstance.edges.findAt(' + curvedCompositeStringEdgesLeft + ')'
#		exec(curvedCompositeStringFacesLeftExec)
#		exec(curvedCompositeStringFacesRightExec)
		# 1) build two Python lists of “1-tuples”:
		angleLimit = 30.0 * np.pi / 180.0
		flankLimitPercentage = 0.1
		yPartitionLimit = -ltan * flankLimitPercentage
		
		inner_edge_coords = []
		
		# Curved section: edges at angleLimit/2 (below) and above angleLimit
		# Edge below angleLimit (in the small angle region)
		inner_edge_coords.append(((ct.pol2cart_x(rk[0], angleLimit/2), 
		                           ct.pol2cart_y(rk[0], angleLimit/2), 
		                           0.0),))
		
		# Edges above angleLimit (in the main curved region)
		for jj in range(len(thetaFaces)):
			if thetaFaces[jj] > angleLimit:
				inner_edge_coords.append(((ct.pol2cart_x(rk[0], thetaFaces[jj]),
				                           ct.pol2cart_y(rk[0], thetaFaces[jj]),
				                           0.0),))
		
		# Straight flank: edges at yPartitionLimit and below
		# Edge above yPartitionLimit (near corner)
		inner_edge_coords.append(((rk[0], yPartitionLimit, 0.0),))
		
		# Edges below yPartitionLimit (main flank region)
		for jj in range(len(flankenFaces)):
			if flankenFaces[jj] < yPartitionLimit:
				inner_edge_coords.append(((rk[0], flankenFaces[jj], 0.0),))
		
		curvedCompositeInnerPressureFaces = LModelInstance.edges.findAt(*inner_edge_coords)
		if boolCylindricalPressure:
			LModel.ExpressionField(name='Load_InnerPressure_CylindricalBending', localCsys=curvedCompositeCylCoordSys, description='',expression='sin(pi*Th/angleOpening)')
			LModel.Pressure(name='Load_InnerPressure_CylindricalBending', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeInnerPressureFaces), magnitude=cylindricalInnerPressure, distributionType=FIELD, field='Load_InnerPressure_CylindricalBending', amplitude=UNSET)
		else:
			LModel.Pressure(name='Load_InnerPressure', createStepName=step, region=regionToolset.Region(side1Edges=curvedCompositeInnerPressureFaces), magnitude=InnerPressure, distributionType=UNIFORM)
	region = regionToolset.Region(edges=curvedCompositeInnerPressureFaces)
	LModel.TemperatureBC(name='Temp', createStepName='Mechanical_Loading', region=region, fixed=OFF, 
		distributionType=UNIFORM, fieldName='', magnitude=20.0, amplitude=UNSET)

#### Interaction ####
angleLimit = 30.0 * np.pi / 180.0
flankLimitPercentage = 0.1
yPartitionLimit = -ltan * flankLimitPercentage

inner_edge_coords = []

# Curved section: edges at angleLimit/2 (below) and above angleLimit
# Edge below angleLimit (in the small angle region)
inner_edge_coords.append(((ct.pol2cart_x(rk[-1], angleLimit/2), 
						ct.pol2cart_y(rk[-1], angleLimit/2), 
						0.0),))

# Edges above angleLimit (in the main curved region)
for jj in range(len(thetaFaces)):
	if thetaFaces[jj] > angleLimit:
		inner_edge_coords.append(((ct.pol2cart_x(rk[-1], thetaFaces[jj]),
								ct.pol2cart_y(rk[-1], thetaFaces[jj]),
								0.0),))

# Straight flank: edges at yPartitionLimit and below
# Edge above yPartitionLimit (near corner)
inner_edge_coords.append(((rk[-1], yPartitionLimit, 0.0),))

# Edges below yPartitionLimit (main flank region)
for jj in range(len(flankenFaces)):
	if flankenFaces[jj] < yPartitionLimit:
		inner_edge_coords.append(((rk[-1], flankenFaces[jj], 0.0),))

curvedCompositeOuter = LModelInstance.edges.findAt(*inner_edge_coords)
region=regionToolset.Region(side1Edges=curvedCompositeOuter)
mdb.models[modelName].FilmCondition(
    name='Int-1', createStepName='Mechanical_Loading', surface=region, 
    definition=EMBEDDED_COEFF, filmCoeff=heatcoeff, filmCoeffAmplitude='', 
    sinkTemperature=300.0, sinkAmplitude='', sinkDistributionType=UNIFORM, 
    sinkFieldName='')

# ---------------------------------------------------------------
if boolBendingMoment and not any([boolRadialForce,boolCircForce]):
	noKinCoup = 1
	for jj in range(noKinCoup,3):
		curvedCompositeLoadRefPointRegion = loadsCreateKinCoup(jj)
		loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion)
elif any([boolBendingMoment,boolRadialForce,boolCircForce]):
	noKinCoup = 1
	for jj in range(noKinCoup):
		curvedCompositeLoadRefPointRegion = loadsCreateKinCoup(jj)
		if boolBendingMoment:
			loadsBendingMoment(jj,curvedCompositeLoadRefPointRegion)
	if boolRadialForce:
		loadsRadialForce()
elif any([boolPressure,boolCylindricalPressure]):
	loadsSurfaceLoad()
else:
	raise ValueError('Please check the loading situation.')

# ---------------------------------------------------------------
# Field Output:
#LModel.FieldOutputRequest(name='F-Output-1', createStepName='Initial', variables=('S', 'U', 'COORD'))

# LModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'U', 'COORD'))

#-------------------------------------------------------------------
# Speichern der CAE-Datei:
'''mdb.saveAs(pathName=analysis_newPath + '\\' + modelName + '_CAE_Model')

#--------------------------------------------------------------
#Job
jobName = modelName + '_Job'

def LmodelJob():
	"""
	Erstellt und führt einen Job für das Modell aus.

	Returns:
		float: Simulationszeit in Sekunden.
	"""
	global startTime
	mdb.Job(
		name=jobName, 
		model=modelName, 
		description='', 
		type=ANALYSIS, 
		atTime=None, 
		waitMinutes=0, 
		waitHours=0, 
		queue=None, 
		memory=99, 
		memoryUnits=PERCENTAGE, 
		getMemoryFromAnalysis=True, 
		explicitPrecision=SINGLE, 
		nodalOutputPrecision=SINGLE, 
		echoPrint=OFF, 
		modelPrint=OFF, 
		contactPrint=OFF, 
		historyPrint=OFF, 
		userSubroutine='', 
		scratch='', 
		resultsFormat=ODB, 
		multiprocessingMode=DEFAULT, 
		numCpus=6, 
		numDomains=6
	)
	startTime = time.time()
	localTime = time.strftime("%H:%M", time.localtime(startTime))
	print('Start:' + localTime + ' Uhr')
	mdb.jobs[jobName].submit(consistencyChecking=OFF)
	mdb.jobs[jobName].waitForCompletion()
	simTime = time.time() - startTime
	return(simTime)

simTime = LmodelJob()'''
