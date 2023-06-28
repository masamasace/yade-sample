

from __future__ import print_function
from yade import pack, plot
from yade.gridpfacet import *
import numpy as np
import datetime
from pathlib import Path
import os
import sys

# Dividing Lines to make the output more comprehensive
print("")
print("-------------------------------------")

############## Constant and Variable ##############
pile_height = 0.5
state_index = 0

############## Contact Model ##############

############## Boundary Box ##############
O.periodic = True
O.cell.hSize = Matrix3(1, 0, 0,
                       0, 3, 0,
                       0, 0, 1)

lowBox = box(center=(0.5, 0.005, 0.5), extents=(5, 0.005, 5),
             fixed=True, wire=False)
O.bodies.append(lowBox)

############## Sphere ##############
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), rMean=.04, rRelFuzz=0, seed=1)
sp.toSimulation()

############## Sphere ##############
cylIds = []
nodesIds = []
pile = cylinder((0.5, 1.25, 0.5),(0.5, 1.25+pile_height, 0.5),
                cylIds=cylIds, nodesIds=nodesIds, radius=0.2, fixed=True)
print(cylIds)

############## Engine ##############
# "allowBiggerThanPeriod=True" in InsertionSortCollider is needed to avoid spheres from falling down out of the box
O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Box_Aabb()], allowBiggerThanPeriod=True),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(gravity=(0, -9.81, 0), damping=0.4),
    PyRunner(iterPeriod=100, command="ShowValue()"),
    PyRunner(iterPeriod=1, command="ChangeState()")
]

O.dt = .5 * PWaveTimeStep()
O.trackEnergy = True

# print(O.bodies[cylIds[0]].shape.node1.shape.dict())
print(nodesIds)

def ShowValue():
    # energy_dict = dict(O.energy.items())
    
    output_str = "UnF: " + '{:> 5.2f}'.format(utils.unbalancedForce()) + \
                 " Pos: " + '{:> 6.3f}'.format(O.bodies[cylIds[0]].state.pos[1]) + \
                 " Fx_cy: " + '{:> 8.2f}'.format(O.forces.f(cylIds[0])[0]) + \
                 " Fy_cy: " + '{:> 8.2f}'.format(O.forces.f(cylIds[0])[1]) + \
                 " Fz_cy: " + '{:> 8.2f}'.format(O.forces.f(cylIds[0])[2]) + \
                 " Fx_no: " + '{:> 8.2f}'.format(O.forces.f(nodesIds[0])[0]) + \
                 " Fy_no: " + '{:> 8.2f}'.format(O.forces.f(nodesIds[0])[1]) + \
                 " Fz_no: " + '{:> 8.2f}'.format(O.forces.f(nodesIds[0])[2])  
    
    print(output_str)
    # print(utils.unbalancedForce(), O.bodies[cylIds[0]].state.pos,O.bodies[cylIds[0]].dynamic)


def ChangeState():
    global state_index
    # print(O.forces.f(cylIds[0]), O.forces.t(cylIds[0]))
    # print(utils.bodyStressTensors()[cylIds[0]])
    # print(O.bodies[cylIds[0]].state.dict())
    # print(utils.unbalancedForce())
    if state_index == 0:
        if utils.unbalancedForce() < 0.1 and O.iter > 10000:
            O.bodies[cylIds[0]].state.blockedDOFs = "xyzXYZ"
            O.bodies[cylIds[0]].state.vel = Vector3(0, -0.2, 0)
            O.bodies[nodesIds[0]].state.blockedDOFs = "xyzXYZ"
            O.bodies[nodesIds[0]].state.vel = Vector3(0, -0.2, 0)
            O.bodies[nodesIds[1]].state.blockedDOFs = "xyzXYZ"
            O.bodies[nodesIds[1]].state.vel = Vector3(0, -0.2, 0)
            state_index = 1
    elif  state_index == 1:
        if O.bodies[cylIds[0]].state.pos[1] <= pile_height / 2:
            O.pause()
        
