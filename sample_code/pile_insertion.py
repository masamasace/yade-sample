

from __future__ import print_function
from yade import pack, plot
from yade.gridpfacet import *
import numpy as np
import datetime
from pathlib import Path
import os
import sys
import json
import math

# Dividing Lines to make the output more comprehensive
print("")
print("-------------------------------------")

dt_start = datetime.datetime.now()
print("Start simulation started at " +
      dt_start.strftime('%Y/%m/%d %H:%M:%S.%f'))


############## Constants ##############
initial_parameters = {
    "pile_radius": 0.125,
    "pile_height": 1,
    "pile_insertion_velocity": -0.02,
    "sphere_diameter_mean": 0.03,
    "sphere_diameter_std_dev": 0,
    "sphere_pack_initial_height": 5,
    "base_box_height_ratio_to_mean_diameter": 5,
    "simulation_box_width": 1,
    "flag_import_existing_pack_file" : False,
    "manual_contact_model": True,
    "flag_output_VTK" : True,
    "export_data_iter_interval" : 100,
    "local_voxel_of_interest": []
}


############## Temporal Variables ##############
state_index = 0

temp_base_box_height = initial_parameters["sphere_diameter_mean"] * initial_parameters["base_box_height_ratio_to_mean_diameter"]

temp_lower_bound = initial_parameters["sphere_diameter_mean"]
temp_lower_bound_y = temp_lower_bound + temp_base_box_height/2
temp_upper_bound = initial_parameters["simulation_box_width"] - initial_parameters["sphere_diameter_mean"]
temp_upper_bound_y = temp_lower_bound_y + initial_parameters["sphere_pack_initial_height"] - initial_parameters["sphere_diameter_mean"]

temp_pile_initial_position_x = initial_parameters["simulation_box_width"] / 2
temp_pile_initial_position_y = temp_upper_bound_y + initial_parameters["sphere_diameter_mean"] + initial_parameters["pile_radius"]
temp_pile_initial_position_z = temp_pile_initial_position_x

temp_hsize_y = math.ceil(temp_pile_initial_position_y + 
                          initial_parameters["pile_height"] + 
                          initial_parameters["pile_radius"])

temp_pause_pile_position_y = temp_base_box_height / 2 + initial_parameters["pile_radius"] + initial_parameters["sphere_diameter_mean"]

temp_initial_pile_disp = 0
# print temporal variables
print("Box width: ", '{:> 4.2f}'.format(initial_parameters["simulation_box_width"]))
print("Base box height: ", '{:> 4.2f}'.format(temp_base_box_height))
print("Initial Y of top of sphere pack: ", '{:> 4.2f}'.format(temp_upper_bound_y))
print("Initial Y of bottom of pile: ", '{:> 4.2f}'.format(temp_pile_initial_position_y))
print("Initial Y of top of pile: ", '{:> 4.2f}'.format(temp_pile_initial_position_y + initial_parameters["pile_height"] + initial_parameters["pile_radius"]))
print("Box height: ", '{:> 4.2f}'.format(temp_hsize_y))
print("Final Y of bottom of pile", '{:> 4.2f}'.format(temp_pause_pile_position_y))


############## Make some directories##############
# Result directory
output_folder_path = Path(os.path.abspath(
    os.path.dirname(sys.argv[0]))).parent / "result" / dt_start.strftime('%Y-%m-%d_%H-%M-%S')
output_folder_path.mkdir(exist_ok=True)
output_file_path = output_folder_path / \
    (dt_start.strftime('%Y-%m-%d_%H-%M-%S') + "_output.csv")
    
# VTK directory
output_VTK_folder_path = output_folder_path / "VTK"
output_VTK_folder_path.mkdir(exist_ok=True)

# Sphere alignment file path
temp_folder_path = Path(os.path.abspath(
    os.path.dirname(sys.argv[0]))).parent / "temp"
temp_folder_path.mkdir(exist_ok=True)
temp_sp_file_path = temp_folder_path / "sphere_pack_pile_insertion.txt"

# Initial parameter json file path
json_file_path = output_folder_path / "initial_parameters.json"
json_file = open(json_file_path, mode="w")
json.dump(initial_parameters, json_file, indent=4)
json_file.close()


############## Contact Model ##############
if initial_parameters["manual_contact_model"]:
    pass


############## Base Boundary Box ##############

base_box_width = initial_parameters["simulation_box_width"]
lowBox = box(center=(base_box_width/2, 0, base_box_width/2), 
             extents=(base_box_width*2, temp_base_box_height/2, base_box_width*2),
             fixed=True, wire=False)
O.bodies.append(lowBox)


############## Sphere ##############
pack_sp = pack.SpherePack()

if initial_parameters["flag_import_existing_pack_file"]:
    pack_sp.load(str(temp_sp_file_path))
    O.periodic = True
else:    
    pack_sp.makeCloud((temp_lower_bound, temp_lower_bound_y, temp_lower_bound), 
                      (temp_upper_bound, temp_upper_bound_y, temp_upper_bound), 
                      rMean=initial_parameters["sphere_diameter_mean"], 
                      rRelFuzz=initial_parameters["sphere_diameter_std_dev"],
                      seed=1,
                      periodic=True)
    pack_sp.save(str(temp_sp_file_path))

pack_sp.toSimulation()
O.periodic = True
O.cell.hSize = Matrix3(initial_parameters["simulation_box_width"], 0, 0,
                       0, temp_hsize_y, 0,
                       0, 0, initial_parameters["simulation_box_width"])

print(len(O.bodies), " spheres are genereted")

############## Pile ##############
cylIds = []
nodesIds = []

pile = cylinder((temp_pile_initial_position_x,
                 temp_pile_initial_position_y,
                 temp_pile_initial_position_z),
                (temp_pile_initial_position_x,
                 temp_pile_initial_position_y+initial_parameters["pile_height"],
                 temp_pile_initial_position_z),
                cylIds=cylIds, nodesIds=nodesIds,
                radius=initial_parameters["pile_radius"], fixed=True)


############## Engine ##############
O.dt = .5 * PWaveTimeStep()

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
    NewtonIntegrator(gravity=(0, -9.81, 0), damping=0.2),
    PyRunner(iterPeriod=initial_parameters["export_data_iter_interval"], command="ExportData()"),
    PyRunner(iterPeriod=1, command="CheckState()")
]


############## Microcodes for exporting dataset ##############
# Enable tracking energy
O.trackEnergy = True

# make an instance for VTK dataset
vtk_recorder = VTKRecorder(fileName=str(output_VTK_folder_path)+'/vtk-', recorders=['spheres', 'intr', 'coordNumber', 'stress', 'force', 'bstresses', 'velocity'])


def ExportData():
    global state_index
    
    i = O.iter
    
    temp_pile_position_y = O.bodies[cylIds[0]].state.pos[1]
    temp_pile_force_x = O.forces.f(cylIds[0])[0]
    temp_pile_force_y = O.forces.f(cylIds[0])[1]
    temp_pile_force_z = O.forces.f(cylIds[0])[2]
    temp_unbalanced_force = utils.unbalancedForce()
    temp_friction_angle = O.materials[0].frictionAngle
    temp_coord_num = utils.avgNumInteractions()
        
    output_str = "UnF: " + '{:> 5.2f}'.format(temp_unbalanced_force) + \
                 " Pos_cy: " + '{:> 6.3f}'.format(temp_pile_position_y) + \
                 " Pos_nd0: " + '{:> 6.3f}'.format(O.bodies[nodesIds[0]].state.pos[1]) + \
                 " Pos_nd1: " + '{:> 6.3f}'.format(O.bodies[nodesIds[1]].state.pos[1]) + \
                 " Fx_cy: " + '{:> 8.2f}'.format(temp_pile_force_x) + \
                 " Fy_cy: " + '{:> 8.2f}'.format(temp_pile_force_y) + \
                 " Fz_cy: " + '{:> 8.2f}'.format(temp_pile_force_z)
    
    print(output_str)
    
    # export some global values to csv file
    output_values = [i, state_index, temp_unbalanced_force, 
                     temp_friction_angle, temp_coord_num,
                     temp_pile_position_y,
                     temp_pile_force_x, 
                     temp_pile_force_y,
                     temp_pile_force_z]

    output_values_str = ""

    for temp in output_values:
        output_values_str += str(temp) + ","
        
    output_values_str = output_values_str[:-1] + "\n"
    with open(output_file_path, 'a') as f:
        f.write(output_values_str)
        
    # export some values to matplotlib figure
    if state_index == 0:
        temp_pile_position_y = None
    plot.addData(i=i,
                 UnF=temp_unbalanced_force,
                 Pos=temp_pile_position_y,
                 Fy_cy=temp_pile_force_y)
    
    # export VTK file
    if initial_parameters["flag_output_VTK"] and state_index != 0:
        vtk_recorder()



def CheckState():
    global state_index, temp_initial_pile_disp

    if state_index == 0:
        if utils.unbalancedForce() < 0.01 and O.iter > 10000:
            
            temp_pile_body_id = cylIds + nodesIds
            temp_body_max_y = max(O.bodies[i].state.pos[1] for i in range(len(O.bodies)) if i not in temp_pile_body_id)
            print("Highest Y position: ", '{:> 4.2f}'.format(temp_body_max_y))
                        
            temp_body_max_y += initial_parameters["sphere_diameter_mean"] + initial_parameters["pile_radius"]
            
            O.bodies[cylIds[0]].state.blockedDOFs = "xyzXYZ"
            O.bodies[cylIds[0]].state.vel = Vector3(0, initial_parameters["pile_insertion_velocity"], 0)
            O.bodies[cylIds[0]].state.pos = Vector3(temp_pile_initial_position_x, temp_body_max_y, temp_pile_initial_position_z)
            
            temp_initial_pile_disp = O.bodies[cylIds[0]].state.displ()[1]
            
            O.bodies[nodesIds[0]].state.blockedDOFs = "xyzXYZ"
            O.bodies[nodesIds[0]].state.vel = Vector3(0, initial_parameters["pile_insertion_velocity"], 0)
            O.bodies[nodesIds[0]].state.pos = Vector3(temp_pile_initial_position_x, temp_body_max_y, temp_pile_initial_position_z)
            
            O.bodies[nodesIds[1]].state.blockedDOFs = "xyzXYZ"
            O.bodies[nodesIds[1]].state.vel = Vector3(0, initial_parameters["pile_insertion_velocity"], 0)
            O.bodies[nodesIds[1]].state.pos = Vector3(temp_pile_initial_position_x, temp_body_max_y+initial_parameters["pile_height"], temp_pile_initial_position_z)
            
            state_index = 1

    elif  state_index == 1:
        
        flag_bottom_reached = O.bodies[cylIds[0]].state.pos[1] <= initial_parameters["sphere_diameter_mean"] + initial_parameters["pile_radius"]
        flag_pile_length = abs(O.bodies[cylIds[0]].state.displ()[1] - temp_initial_pile_disp) >= initial_parameters["pile_height"]
        
        if flag_bottom_reached or flag_pile_length:
            O.pause()
        
plot.plots = {"i": ("UnF"),
              "Fy_cy ": ("Pos")}

plot.plot()