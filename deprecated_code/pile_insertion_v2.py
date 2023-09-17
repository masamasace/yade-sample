from yade import pack, utils
import numpy as np

state_index = 0

pack_sp = pack.SpherePack()
pack_sp.makeCloud((0.02, 0.02, 0.02), 
                  (0.08, 0.98, 0.08), 
                  rMean=0.01)
sphere_id = pack_sp.toSimulation()

O.periodic = True

base_facet_1 = utils.facet([Vector3(-0.02, 0, -0.02), 
                            Vector3(0.12,  0, -0.02),
                            Vector3(-0.02, 0, 0.12)],
                            fixed=True)
base_facet_2 = utils.facet([Vector3(0.12, 0,  0.12), 
                            Vector3(0.12,  0, -0.02),
                            Vector3(-0.02, 0, 0.12)],
                            fixed=True)

base_facet_id = O.bodies.append([base_facet_1, base_facet_2])

O.cell.hSize = Matrix3(0.1, 0, 0,
                       0, 1, 0,
                       0, 0, 0.1)
O.periodic = True

O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()], allowBiggerThanPeriod=True),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]),
    NewtonIntegrator(damping=0.2, gravity=(0, -9.81, 0)),
    PyRunner(iterPeriod=100, command="checkState()")
    ]

O.dt = PWaveTimeStep() * 0.5

def checkState():
    global sphere_id, state_index
    
    if state_index == 0:
        
        temp_sphere_coord_num = [len(O.interactions.withBody(i)) for i in sphere_id]
        temp_sphere_Y = [O.bodies[i].state.pos[1] for i in sphere_id]
        temp_sphere_data = np.array([sphere_id, temp_sphere_coord_num, temp_sphere_Y])
        
        if temp_sphere_data[1, :].mean() >=3:
            
            temp_unused_sphere_id = temp_sphere_data[0, temp_sphere_data[2, :] > 0.1]
            
            for temp_unused_sphere_id_each in temp_unused_sphere_id:
                O.bodies.erase(int(temp_unused_sphere_id_each))
                sphere_id.remove(int(temp_unused_sphere_id_each))
            
            print("some particles are deleted at", O.iter)
            
            state_index = 1
    
    elif state_index == 1:
        pass