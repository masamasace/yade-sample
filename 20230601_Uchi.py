from __future__ import print_function
from yade import pack,plot,polyhedra_utils,geom
from yade import export,qt
from yade.gridpfacet import *
import math
from yade import plot
import numpy as np
import os

frictangle = np.radians(35)#45do
density = 3700.0
young = 3e8
gravity = (0.0, 0.0, 0)

	
mat_sp = CohFrictMat(alphaKr=0.5,alphaKtw=0.5,young = young, poisson = 0.33,frictionAngle=(frictangle), density=density)
O.materials.append(mat_sp)
	
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (6, 6, 2), rMean=0.10, rRelFuzz=0, periodic=True) 
	# insert the packing into the simulation
sp.toSimulation(color=(0, 0, 1))  # pure blue

O.materials[0].frictionAngle = 0

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb(),Bo1_Box_Aabb()]),
        InteractionLoop(
                # interaction loop
                [Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()],
                [Law2_ScGeom6D_CohFrictPhys_CohesionMoment(always_use_moment_law=True),Law2_ScGeom_FrictPhys_CundallStrack()]
        ),
        NewtonIntegrator(gravity=gravity,damping=0.2),
        PyRunner(command='checkStress1()', realPeriod=1, label='checker'),
        PyRunner(command='checkStress2()', realPeriod=1, label='checker'),
        PyRunner(iterPeriod=1,command="addPlotData()")
        ]
O.dt = .1 * PWaveTimeStep()

plot.plots = {"strainXZ":"stressXZ"}              
def checkStress1():
	stress1 = getStress()[0,0]
	stress2 = getStress()[1,1]
	stress3 = getStress()[2,2]
	print(stress1,stress2,stress3)
	
	
O.cell.velGrad = Matrix3(-10,0,0, 0,-10,0, 0,0,-10)
def checkStress2():
	if abs(getStress()[0,0]) > 1e5:
		checker.command = 'checkDistorsion1()'
	
def checkDistorsion1():
	# positive shear
	O.cell.velGrad =  Matrix3(0,0,+10, 0,0,0, 0,0,0)
	if (O.cell.trsf[0, 2]) > 0.40:
		checker.command = 'checkDistorsion2()'
		
def checkDistorsion2():
	# negative shear
	O.cell.velGrad =  Matrix3(0,0,-10, 0,0,0, 0,0,0)
	if (O.cell.trsf[0, 2]) < -0.40:
		checker.command = 'checkDistorsion3()'

def checkDistorsion3():
	# positive shear
	O.cell.velGrad =  Matrix3(0,0,+10, 0,0,0, 0,0,0)
	if (O.cell.trsf[0, 2]) > 0.40:
		checker.command = 'checkDistorsion4()'
		
def checkDistorsion4():
	# negative shear
	O.cell.velGrad =  Matrix3(0,0,-10, 0,0,0, 0,0,0)
	if (O.cell.trsf[0, 2]) < -0.40:
		O.pause()
		#checker.command = 'checkDistorsion5()'	
			
def addPlotData():
    stress = getStress()
    plot.addData(
        strainXZ = O.cell.trsf[0,2],
        stressXZ = stress[0,2],
    )
plot.plots = {"strainXZ":"stressXZ"}

plot.plot()

