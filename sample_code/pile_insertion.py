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

############## Contact Model ##############

############## Boundary Box ##############
# TODO: 周期境界にしたい
O.periodic = True
O.cell.hSize = Matrix3(1, 0, 0,
                       0, 1, 0,
                       0, 0, 1)

lowBox = box(center=(0.5, 0.005, 0.5), extents=(5, 0.005, 5),
             fixed=True, wire=False)
O.bodies.append(lowBox)


############## Sphere ##############
sp = pack.SpherePack()
sp.makeCloud((0.05, 0.5, 0.05), (0.95, 0.95, 0.95), rMean=.02, rRelFuzz=0)
sp.toSimulation()

############## Engine ##############
# "allowBiggerThanPeriod=True" is needed to avoid particle falling down out of the box
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
    PyRunner(iterPeriod=10, command="checkUnbalncedForce()"),
]

O.dt = .5 * PWaveTimeStep()
O.trackEnergy = True


def checkUnbalncedForce():
    print(utils.unbalancedForce())
