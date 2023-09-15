from __future__ import print_function
from yade import pack, utils

mat = FrictMat()
mat_id = O.materials.append(mat)

base_facet_1 = utils.facet([Vector3(-0.3, 0.03, -0.3), 
                             Vector3(0.2, 0.03, -0.3),
                             Vector3(-0.3, 0.03, 0.2)],
                            fixed=True,
                            wire=False,
                            material=mat_id)

base_facet_id = O.bodies.append([base_facet_1])


pack_sp = pack.SpherePack()
pack_sp.makeCloud((0, 0.02, 0), 
                  (0.1, 0.98, 0.1), 
                  rMean=0.02)
sphere_id = pack_sp.toSimulation(material=mat_id)

print(O.bodies[0].dict())
print(O.bodies[1].dict())


O.periodic = True
O.cell.hSize = Matrix3(0.1, 0, 0,
                       0, 1, 0,
                       0, 0, 0.1)

O.dt = PWaveTimeStep() / 20

O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()], allowBiggerThanPeriod=True),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]),
    NewtonIntegrator(damping=0.2, gravity=(0, -9.81, 0))
    ]

O.run()