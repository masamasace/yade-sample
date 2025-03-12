############## 基本的なモジュールのインポート ##############
from __future__ import print_function
from yade import pack, plot, O, Matrix3, utils
from yade.wrapper import FrictMat, ForceResetter, \
    InsertionSortCollider, InteractionLoop, NewtonIntegrator, PyRunner
import numpy as np
import datetime
from pathlib import Path
import os
import sys
import json

############## シミュレーションの設定 ##############
# given constants
SEED = 1
TARGET_VOID_RATIO = 0.6
CONS_STRAIN_RATE = 1e-3

## simulation cell size
CELL_SIZE = 1.0

## sphere parameters
SPHERE_RADIUS = 0.05

## model parameters
CONTACT_YOUNG_MODULUS = 1e7
CONTACT_POISSON_RATIO = 0.3
CONTACT_FRICTION_ANGLE = np.radians(30)
SPHERE_DENSITY = 2650

# derived constants
TARGET_POROSITY = TARGET_VOID_RATIO / (1 + TARGET_VOID_RATIO)
CONS_VELGRAD = Matrix3(-CONS_STRAIN_RATE, 0, 0,
                        0, -CONS_STRAIN_RATE, 0,
                        0, 0, -CONS_STRAIN_RATE)

# set model
mat_sp = FrictMat(young=CONTACT_YOUNG_MODULUS, 
                  poisson=CONTACT_POISSON_RATIO, 
                  frictionAngle=CONTACT_FRICTION_ANGLE, 
                  density=SPHERE_DENSITY)
O.materials.append(mat_sp)

# set sphere
pack_sp = pack.SpherePack()
pack_sp.makeCloud(minCorner=(0, 0, 0), 
                  maxCorner=(CELL_SIZE, CELL_SIZE, CELL_SIZE), 
                  rMean=SPHERE_RADIUS, 
                  periodic=True,
                  seed=SEED)
pack_sp.toSimulation(material=mat_sp)

# set PWaveTimeStepper
O.dt = 0.1 * utils.PWaveTimeStep()

# set engines
O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(damping=0.2),
    PyRunner(iterPeriod=1, command="checkState()")
]

# change settings
## change velocity gradient
O.cell.velGrad = CONS_VELGRAD

## change track energy
O.trackEnergy = True

# define functions

## check state
def checkState():

    # get stress and strain data
    [s00, s01, s02], [s10, s11, s12], [s20, s21, s22] = utils.getStress()
    [e00, e01, e02], [e10, e11, e12], [e20, e21, e22] = O.cell.trsf

    #     
############## シミュレーションの実行 ##############

