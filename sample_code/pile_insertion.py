

from __future__ import print_function
from yade import pack, plot, ymport
from yade.gridpfacet import *
import numpy as np
import datetime
from pathlib import Path
import os
import sys
import json
import math
import collections

# Dividing line to make the output more comprehensive
print("")
print("-------------------------------------")

dt_start = datetime.datetime.now()
print("Start simulation started at " +
      dt_start.strftime('%Y/%m/%d %H:%M:%S.%f'))


############## Constants ##############
initial_parameters = {
    "pile_height": 0.3,            # Temporary value to be updated after loading stl file
    "pile_insertion_velocity": -0.02,
    "flag_uniform_sphere_diameter": False,
    "sphere_diameter_mean": 0.003,
    "sphere_diameter_std_dev": 0,
    "sphere_size_distribution_size": [0.000050, 0.000075, 0.000106, 0.000250, 0.000425], # Ref: 平成25年度地盤材料試験の技能試験報告書
    "sphere_size_distribution_size_ratio": 20,
    "sphere_size_distribution_cumm": [0.0     , 0.04    , 0.12    , 0.88    , 1.0     ], # Ref: 平成25年度地盤材料試験の技能試験報告書
    "sphere_density": 2650.0,
    "manual_contact_model": True,
    "sphere_contact_normal_stiffness": 8.0e7,
    "sphere_contact_stiffness_ratio": 0.25,
    "sphere_contact_frictional_coeffcient": 0.5,
    "pile_contact_normal_stiffness": 6.0e12,
    "pile_contact_stiffness_ratio": 0.25,
    "pile_contact_frictional_coeffcient": 0.5,
    "sphere_fall_height": 0.25,
    "sphere_pack_initial_height": 0.5,
    "base_box_height_ratio_to_mean_diameter": 5,
    "simulation_box_width": 0.04,
    "flag_import_existing_pack_file" : False,
    "flag_import_heavy_stl_model": True,
    "flag_output_VTK" : True,
    "check_state_iter_interval" : 50,
    "export_data_iter_interval" : 100,
    "local_voxel_of_interest": []
}


############## Temporary Variables ##############
state_index = 0
temp_sphere_size_distribution_size = [ssd*initial_parameters["sphere_size_distribution_size_ratio"] for ssd in initial_parameters["sphere_size_distribution_size"]]

if initial_parameters["flag_uniform_sphere_diameter"]:
    temp_max_sphere_size = initial_parameters["sphere_diameter_mean"]
    temp_min_sphere_size = initial_parameters["sphere_diameter_mean"]
else:
    temp_max_sphere_size = max(temp_sphere_size_distribution_size)
    temp_min_sphere_size = min(temp_sphere_size_distribution_size)

temp_base_box_height = temp_max_sphere_size * initial_parameters["base_box_height_ratio_to_mean_diameter"]


temp_lower_bound = 0
# temp_lower_bound = temp_max_sphere_size
temp_lower_bound_y = temp_max_sphere_size / 2
# temp_lower_bound_y = temp_lower_bound * 2 + temp_base_box_height / 2

temp_upper_bound = initial_parameters["simulation_box_width"]
# temp_upper_bound = initial_parameters["simulation_box_width"] - temp_max_sphere_size
temp_upper_bound_y = temp_lower_bound_y + initial_parameters["sphere_pack_initial_height"]
# temp_upper_bound_y = temp_lower_bound_y + initial_parameters["sphere_pack_initial_height"] - temp_max_sphere_size

temp_pile_initial_position_x = initial_parameters["simulation_box_width"] / 2
temp_pile_initial_position_y = temp_upper_bound_y + temp_max_sphere_size
temp_pile_initial_position_z = temp_pile_initial_position_x

temp_hsize_y = math.ceil((temp_pile_initial_position_y + initial_parameters["pile_height"]) * 10) / 10

temp_pause_pile_position_y = temp_base_box_height / 2 + temp_max_sphere_size

temp_initial_pile_disp = 0


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
    frict_mat_sphere = FrictMat(young=initial_parameters["sphere_contact_normal_stiffness"], 
                               poisson=initial_parameters["sphere_contact_stiffness_ratio"],
                               frictionAngle=initial_parameters["sphere_contact_frictional_coeffcient"], 
                               density=initial_parameters["sphere_density"])
    material_sphere_Id = O.materials.append(frict_mat_sphere)
    
    frict_mat_facet = FrictMat(young=initial_parameters["pile_contact_normal_stiffness"], 
                              poisson=initial_parameters["pile_contact_stiffness_ratio"],
                              frictionAngle=initial_parameters["pile_contact_frictional_coeffcient"], 
                              density=1000)
    material_facet_Id = O.materials.append(frict_mat_facet)


############## Base Boundary Box ##############
# utils.wallを使う方法もある

base_facet_1 = utils.facet([Vector3(0, 0, 0), 
                             Vector3(initial_parameters["simulation_box_width"], 0, 0),
                             Vector3(0, 0, initial_parameters["simulation_box_width"])],
                            fixed=True,
                            wire=False,
                            material=material_facet_Id)
base_facet_2 = utils.facet([Vector3(initial_parameters["simulation_box_width"], 0, initial_parameters["simulation_box_width"]), 
                             Vector3(initial_parameters["simulation_box_width"], 0, 0),
                             Vector3(0, 0, initial_parameters["simulation_box_width"])],
                            fixed=True,
                            wire=False,
                            material=material_facet_Id)
top_facet_1 = utils.facet([Vector3(0, temp_pile_initial_position_y, 0), 
                             Vector3(initial_parameters["simulation_box_width"]*2, temp_pile_initial_position_y, 0),
                             Vector3(0, temp_pile_initial_position_y, initial_parameters["simulation_box_width"]*2)],
                            fixed=True,
                            wire=False,
                            material=material_facet_Id)
base_top_facet_id = O.bodies.append([base_facet_1, base_facet_2, top_facet_1])


############## Sphere ##############
pack_sp = pack.SpherePack()

if initial_parameters["flag_import_existing_pack_file"]:
    pack_sp.load(str(temp_sp_file_path))
    O.periodic = True
else: 
    if initial_parameters["flag_uniform_sphere_diameter"]:   
        pack_sp.makeCloud((temp_lower_bound, temp_lower_bound_y, temp_lower_bound), 
                        (temp_upper_bound, temp_upper_bound_y, temp_upper_bound), 
                        rMean=initial_parameters["sphere_diameter_mean"], 
                        rRelFuzz=initial_parameters["sphere_diameter_std_dev"],
                        seed=1,
                        periodic=True)
        pack_sp.save(str(temp_sp_file_path))
    else:
        
        pack_sp.makeCloud((temp_lower_bound, temp_lower_bound_y, temp_lower_bound), 
                        (temp_upper_bound, temp_upper_bound_y, temp_upper_bound), 
                        psdSizes=temp_sphere_size_distribution_size,
                        psdCumm=initial_parameters["sphere_size_distribution_cumm"],
                        seed=1,num=3000,
                        periodic=True,)
        
        pack_sp.save(str(temp_sp_file_path))


sphere_id = pack_sp.toSimulation(material=material_sphere_Id)

O.periodic = True
O.cell.hSize = Matrix3(initial_parameters["simulation_box_width"], 0, 0,
                       0, temp_hsize_y, 0,
                       0, 0, initial_parameters["simulation_box_width"])

print(len(sphere_id), "spheres are genereted")


############## Pile ##############
if not initial_parameters["flag_import_heavy_stl_model"]:
    stl_file_path = r"./temp/pile_v2_light.stl"
    
else:
    stl_file_path = r"./temp/pile_v1_heavy.stl"

pile_facets = ymport.stl(stl_file_path, fixed=True, 
                        wire=True, noBound=False, 
                        material=material_facet_Id, 
                        scale=1/1000, 
                        shift=Vector3(temp_pile_initial_position_x, temp_pile_initial_position_y, temp_pile_initial_position_z))
print(len(pile_facets), "facets of pile model are imported")

temp_pile_facets_id = O.bodies.append(pile_facets)
temp_pile_facets_pos_Y = [each_pile_facets.state.pos[1] for each_pile_facets in pile_facets]
pile_facet_data = np.array([temp_pile_facets_id, temp_pile_facets_pos_Y]).T

pile_facets_pos_Y_top = pile_facet_data[:, 1].max()
pile_facets_pos_Y_bottom = pile_facet_data[:, 1].min()
initial_parameters["pile_height"] = pile_facets_pos_Y_top - pile_facets_pos_Y_bottom

print("Pile is initially located exists across coordinates", 
      '{:4.2f}'.format(pile_facets_pos_Y_bottom),
      "to",
      '{:4.2f}'.format(pile_facets_pos_Y_top))

temp_flag_pile_facets_on_top = pile_facet_data[:, 1] == pile_facets_pos_Y_top
temp_flag_pile_facets_at_bottom = pile_facet_data[:, 1] == pile_facets_pos_Y_bottom
temp_flag_pile_facets_on_lateral = ~(temp_flag_pile_facets_on_top | temp_flag_pile_facets_at_bottom)

pile_facets_id_top = np.array(pile_facet_data[temp_flag_pile_facets_on_top, 0], dtype=np.int32)
pile_facets_id_bottom = np.array(pile_facet_data[temp_flag_pile_facets_at_bottom, 0], dtype=np.int32)
pile_facets_id_lateral = np.array(pile_facet_data[temp_flag_pile_facets_on_lateral, 0], dtype=np.int32)

# update simulation box size
temp_hsize_y = math.ceil(pile_facets_pos_Y_top * 10) / 10
O.cell.hSize = Matrix3(initial_parameters["simulation_box_width"], 0, 0,
                       0, temp_hsize_y, 0,
                       0, 0, initial_parameters["simulation_box_width"])


############## Parameter Check ##############
print("Simulation box width:", '{:> 4.2f}'.format(initial_parameters["simulation_box_width"]))
print("Base box height:", '{:> 4.2f}'.format(temp_base_box_height))
print("Initial Y of top of sphere pack:", '{:> 4.2f}'.format(temp_upper_bound_y))
print("Initial Y of bottom of pile:", '{:> 4.2f}'.format(temp_pile_initial_position_y))
print("Initial Y of top of pile:", '{:> 4.2f}'.format(pile_facets_pos_Y_top))
print("Simulation box height:", '{:> 4.2f}'.format(temp_hsize_y))
print("Final Y of bottom of pile:", '{:> 4.2f}'.format(temp_pause_pile_position_y))


############## Engine ##############
O.dt = .5 * PWaveTimeStep()

# "allowBiggerThanPeriod=True" in InsertionSortCollider is needed to avoid spheres from falling down out of the box
# Though this is periodical simulation, Bo1_Box_Aabb is still needed to prevent spheres from falling out of the bottom base
O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()], allowBiggerThanPeriod=True),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(damping=0.2, gravity=(0, -9.81, 0)),
    PyRunner(iterPeriod=initial_parameters["export_data_iter_interval"], command="exportData()"),
    PyRunner(iterPeriod=1, command="checkState()")
]


############## Microcodes for exporting dataset ##############
# Enable tracking energy
O.trackEnergy = True

# make an instance for VTK dataset
vtk_recorder = VTKRecorder(fileName=str(output_VTK_folder_path)+'/vtk-', recorders=['spheres', 'facets', 'intr', 'coordNumber', 'stress', 'force', 'bstresses', 'velocity'])


############## Subroutine ##############
def exportData():
    global state_index
    
    iter = O.iter
    temp_pile_position_y, temp_pile_force_x, temp_pile_force_y, temp_pile_force_z = calcuratePileMechanicalValues()
    temp_unbalanced_force = utils.unbalancedForce()
    temp_friction_angle = O.materials[0].frictionAngle
    temp_coord_num = utils.avgNumInteractions()
        
    output_str = "Iter: " + str(iter) + \
                 " Stage: " + str(state_index) + \
                 " UnF: " + '{:> 5.2f}'.format(temp_unbalanced_force) + \
                 " Cn: " + '{:> 6.3f}'.format(temp_coord_num) + \
                 " Pos_cy: " + '{:> 6.3f}'.format(temp_pile_position_y) + \
                 " Fx_cy: " + '{:> 8.2f}'.format(temp_pile_force_x) + \
                 " Fy_cy: " + '{:> 8.2f}'.format(temp_pile_force_y) + \
                 " Fz_cy: " + '{:> 8.2f}'.format(temp_pile_force_z)
    
    print(output_str)
    
    # export some global values to csv file
    output_values = [iter, state_index, temp_unbalanced_force, 
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
    if state_index != 2:
        temp_pile_position_y = None
        
    plot.addData(i=iter,
                 UnF=temp_unbalanced_force,
                 Pos=temp_pile_position_y,
                 Fy_cy=temp_pile_force_y)
    
    # export VTK file
    if initial_parameters["flag_output_VTK"] and state_index != 0:
        vtk_recorder()


def calcuratePileMechanicalValues(debug=False):
    
    pile_force = Vector3(0, 0, 0)
    
    # maybe faster when using utils.sumForces((list)ids, (Vector3)direction?)
    for pile_facet_id in pile_facet_data[:, 0]:
        pile_force += O.forces.f(int(pile_facet_id))
        if debug:
            print(O.forces.f(int(pile_facet_id)))
    
    temp_pile_position_y = O.bodies[int(pile_facets_id_bottom[0])].state.pos[1]
    return (temp_pile_position_y, pile_force[0], pile_force[1], pile_force[2])


def checkState():
    global state_index, temp_initial_pile_disp

    # initial particle diameter expansion without gravity force
    # expanding radius 
    if state_index == 0:
        temp_body_max_y = max(O.bodies[i].state.pos[1] for i in sphere_id)
        temp_coord_num = utils.avgNumInteractions()
        
        if O.iter % 100 == 0:
            print(O.forces.f(base_top_facet_id[2]))
        
        if temp_coord_num > 1.5:
            for i in sphere_id:
                O.bodies[i].state.vel = Vector3(0, 0, 0)
            state_index = 1
        else:
            utils.growParticles(1.0001)

            
    # gravity depostion stage (air pluviation stage)
    elif state_index == 1:
        if O.iter % 100 == 0:
            print(O.forces.f(base_top_facet_id[2]))

        if utils.unbalancedForce() < 0.01:
            O.pause()
        
    
    # wait until the assembly is stabilized
    elif state_index == 2:
        
        if utils.unbalancedForce() < 0.05 and O.iter > 10000:
            
            temp_sphere_velocity_max = max(O.bodies[i].state.vel.norm() for i in range(sphere_id_max))
            
            if temp_sphere_velocity_max < 0.1:
                temp_body_max_y = max(O.bodies[i].state.pos[1] for i in range(sphere_id_max))
                print("Highest Y position: ", '{:> 4.2f}'.format(temp_body_max_y))
                            
                temp_body_max_y += temp_max_sphere_size 
                temp_initial_pile_disp = O.bodies[int(pile_facets_id_bottom[0])].state.pos[1]
                
                temp_offset = temp_initial_pile_disp - temp_body_max_y
                
                for pile_facet_id in pile_facet_data[:, 0]:
                    O.bodies[int(pile_facet_id)].state.blockedDOFs = "xyzXYZ"
                    O.bodies[int(pile_facet_id)].state.vel = Vector3(0, initial_parameters["pile_insertion_velocity"], 0)
                    
                    temp_pile_facet_pos_each = O.bodies[int(pile_facet_id)].state.pos
                    temp_pile_facet_pos_each[1] -= temp_offset 
                    O.bodies[int(pile_facet_id)].state.pos = temp_pile_facet_pos_each
                
                temp_initial_pile_disp = temp_body_max_y
                
                state_index = 2

    elif  state_index == 3:
        
        flag_bottom_reached = O.bodies[int(pile_facets_id_bottom[0])].state.pos[1] <= temp_max_sphere_size
        flag_pile_length = abs(temp_initial_pile_disp - O.bodies[int(pile_facets_id_bottom[0])].state.pos[1]) >= initial_parameters["pile_height"]
        
        if flag_bottom_reached or flag_pile_length:
            O.pause()


plot.plots = {"i": ("UnF"),
              "Fy_cy ": ("Pos")}

plot.plot()