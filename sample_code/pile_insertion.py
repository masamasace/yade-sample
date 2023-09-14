

from __future__ import print_function
from yade import pack, plot, ymport, export
from yade.gridpfacet import *
import numpy as np
import datetime
from pathlib import Path
import os
import sys
import json
import math


# Dividing line to make the output more comprehensive
print("")
print("-------------------------------------")

dt_start = datetime.datetime.now()
print("Start simulation started at " +
      dt_start.strftime('%Y/%m/%d %H:%M:%S.%f'))


############## Pre-difined Subroutines ##############

def calcurateVoxelRegion(center_line_x_of_interest_local_voxel, size_of_interest_local_voxel, sphere_pack_target_Y, simulation_box_width, total_num_of_local_voxel=-1, debug=False):

    temp_num_voxel_along_line = int(sphere_pack_target_Y // size_of_interest_local_voxel)
    temp_voxel_y_min = [y*size_of_interest_local_voxel for y in range(temp_num_voxel_along_line)]
    
    if total_num_of_local_voxel != -1:
        temp_voxel_y_min_id = list(np.linspace(0, len(temp_voxel_y_min), total_num_of_local_voxel + 2).astype('int32')[1:-1])
    else:
        temp_voxel_y_min_id = list(range(len(temp_voxel_y_min)))
        
    target_region = []
    
    for i in range(len(center_line_x_of_interest_local_voxel)):
        for temp_voxel_y_min_id_each in temp_voxel_y_min_id:
            target_region_each_x_min = center_line_x_of_interest_local_voxel[i] - size_of_interest_local_voxel / 2
            target_region_each_x_max = center_line_x_of_interest_local_voxel[i] + size_of_interest_local_voxel / 2
            target_region_each_y_min = temp_voxel_y_min[temp_voxel_y_min_id_each]
            target_region_each_y_max = temp_voxel_y_min[temp_voxel_y_min_id_each] + size_of_interest_local_voxel
            target_region_each_z_min = simulation_box_width / 2 - size_of_interest_local_voxel / 2
            target_region_each_z_max = simulation_box_width / 2 + size_of_interest_local_voxel / 2
            target_region_each = [[target_region_each_x_min, target_region_each_y_min, target_region_each_z_min],
                                  [target_region_each_x_max, target_region_each_y_max, target_region_each_z_max]]
            target_region.append(target_region_each)
    
    print(len(target_region), "regional voxels are created:")
    for target_region_each in target_region:
        print(target_region_each)
    
        
    return target_region


############## Constants ##############
initial_parameters = {
    "pile_height": 0.3,            # Temporary value to be updated after loading stl file
    "pile_insertion_velocity": -0.2,
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
    "sphere_pack_target_Y": 0.5,
    "base_facet_Y_ratio_to_mean_diameter": 5,
    "simulation_box_width": 0.04,
    "flag_import_existing_pack_file" : False,
    "flag_import_heavy_stl_model": True,
    "flag_output_VTK" : True,
    "check_state_iter_interval" : 50,
    "export_data_iter_interval" : 1000,
    "center_line_x_of_interest_local_voxel": [0.01],
    "size_of_interest_local_voxel": 0.01,
    "total_num_of_local_voxel": 4
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

temp_base_facet_Y = temp_max_sphere_size / 2

temp_sphere_generation_bound_lower_xz = 0
# temp_sphere_generation_bound_lower_xz = temp_max_sphere_size
temp_sphere_generation_bound_lower_y = temp_base_facet_Y + temp_max_sphere_size / 2
# temp_sphere_generation_bound_lower_y = temp_sphere_generation_bound_lower * 2 + temp_base_facet_height / 2

temp_upper_bound = initial_parameters["simulation_box_width"]
# temp_upper_bound = initial_parameters["simulation_box_width"] - temp_max_sphere_size
temp_upper_bound_y = temp_sphere_generation_bound_lower_y + initial_parameters["sphere_pack_target_Y"]
# temp_upper_bound_y = temp_sphere_generation_bound_lower_y + initial_parameters["sphere_pack_target_Y"] - temp_max_sphere_size

temp_pile_initial_position_x = initial_parameters["simulation_box_width"] / 2
temp_pile_initial_position_y = temp_upper_bound_y + temp_max_sphere_size
temp_pile_initial_position_z = temp_pile_initial_position_x

temp_pile_top_height = temp_pile_initial_position_y + initial_parameters["pile_height"]
temp_pile_top_height_ceiled = math.ceil(temp_pile_top_height * 10) / 10
temp_hsize_y = temp_pile_top_height_ceiled

temp_voxel_region = calcurateVoxelRegion(initial_parameters["center_line_x_of_interest_local_voxel"], 
                                         initial_parameters["size_of_interest_local_voxel"],
                                         initial_parameters["sphere_pack_target_Y"],
                                         initial_parameters["simulation_box_width"],
                                         total_num_of_local_voxel=initial_parameters["total_num_of_local_voxel"],
                                         debug=True)

sphere_id = []

temp_prev_stage_iter = 0


############## Initilize periodic cell ##############

O.periodic = True
O.cell.hSize = Matrix3(initial_parameters["simulation_box_width"], 0, 0,
                       0, temp_hsize_y, 0,
                       0, 0, initial_parameters["simulation_box_width"])

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
    
    frict_mat_facet_pile = FrictMat(young=initial_parameters["pile_contact_normal_stiffness"], 
                              poisson=initial_parameters["pile_contact_stiffness_ratio"],
                              frictionAngle=initial_parameters["pile_contact_frictional_coeffcient"], 
                              density=1000)
    material_facet_pile_Id = O.materials.append(frict_mat_facet_pile)
    
    frict_mat_facet_base = FrictMat(young=initial_parameters["pile_contact_normal_stiffness"] * 10, 
                              poisson=initial_parameters["pile_contact_stiffness_ratio"] * 10,
                              frictionAngle=initial_parameters["pile_contact_frictional_coeffcient"] * 10, 
                              density=1000)
    material_facet_base_Id = O.materials.append(frict_mat_facet_base)


############## Base Boundary Box ##############
temp_base_facet_x_min = -temp_max_sphere_size / 2
temp_base_facet_z_min = temp_base_facet_x_min
temp_base_facet_x_max = temp_max_sphere_size / 2 * math.sqrt(2) + 2 * initial_parameters["simulation_box_width"]
temp_base_facet_z_max = temp_base_facet_x_max

base_facet_1 = utils.facet([Vector3(temp_base_facet_x_min, temp_base_facet_Y, temp_base_facet_z_min), 
                             Vector3(temp_base_facet_x_max, temp_base_facet_Y, temp_base_facet_z_min),
                             Vector3(temp_base_facet_x_min, temp_base_facet_Y, temp_base_facet_z_max)],
                            fixed=True,
                            wire=False,
                            material=material_facet_base_Id)
""
"""base_facet_2 = utils.facet([Vector3(initial_parameters["simulation_box_width"], temp_base_facet_Y, initial_parameters["simulation_box_width"]), 
                             Vector3(initial_parameters["simulation_box_width"], temp_base_facet_Y, 0),
                             Vector3(0, temp_base_facet_Y, initial_parameters["simulation_box_width"])],
                            fixed=True,
                            wire=False,
                            material=material_facet_base_Id)"""

base_top_facet_id = O.bodies.append([base_facet_1])


############## Sphere ##############

pack_sp = pack.SpherePack()
pack_sp.makeCloud((temp_sphere_generation_bound_lower_xz, temp_sphere_generation_bound_lower_y, temp_sphere_generation_bound_lower_xz), 
                  (initial_parameters["simulation_box_width"], initial_parameters["sphere_pack_target_Y"]+temp_max_sphere_size, initial_parameters["simulation_box_width"]), 
                  psdSizes=temp_sphere_size_distribution_size,
                  psdCumm=initial_parameters["sphere_size_distribution_cumm"],
                  seed=1,
                  periodic=True)
sphere_id = pack_sp.toSimulation(material=material_sphere_Id)

O.periodic = True
O.cell.hSize = Matrix3(initial_parameters["simulation_box_width"], 0, 0,
                       0, temp_hsize_y, 0,
                       0, 0, initial_parameters["simulation_box_width"])

print(len(sphere_id), "spheres are genereted")


############## Pile ##############
if not initial_parameters["flag_import_heavy_stl_model"]:
    stl_file_path = r"../stl/pile_v2_light.stl"
    
else:
    stl_file_path = r"../stl/pile_v1_heavy.stl"

stl_file_path = Path(stl_file_path).resolve()
print("Pile stl is found in", str(stl_file_path))
stl_file_path = str(stl_file_path)

pile_facets = ymport.stl(stl_file_path, fixed=True, 
                        wire=True, noBound=True, 
                        material=material_facet_pile_Id, 
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
print("Base Facet Y:", '{:> 4.2f}'.format(temp_base_facet_Y))
print("Initial Y of top of sphere pack:", '{:> 4.2f}'.format(temp_upper_bound_y))
print("Initial Y of bottom of pile:", '{:> 4.2f}'.format(temp_pile_initial_position_y))
print("Initial Y of top of pile:", '{:> 4.2f}'.format(pile_facets_pos_Y_top))
print("Simulation box height:", '{:> 4.2f}'.format(temp_hsize_y))


############## Engine ##############
O.dt = PWaveTimeStep() / 2

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
    
    temp_sphere_vel = [O.bodies[i].state.vel[1] for i in sphere_id]
    temp_sphere_Y = [O.bodies[i].state.pos[1] for i in sphere_id]
    temp_sphere_coord_num = [len(O.interactions.withBody(i)) for i in sphere_id]
    temp_sphere_data = np.array([sphere_id, temp_sphere_vel, temp_sphere_Y, temp_sphere_coord_num])
    
    # export some values to matplotlib figure
    if state_index <= 2:
        temp_pile_position_y = 0

    temp_voxel_porosity = [utils.voxelPorosity(start=temp_voxel_region_each[0], end=temp_voxel_region_each[1]) for temp_voxel_region_each in temp_voxel_region]

    if len(temp_sphere_data[2, temp_sphere_data[3, :]>=3]) != 0:
        temp_maxY_with_3cn = temp_sphere_data[2, temp_sphere_data[3, :]>=3].max()
    else:
        temp_maxY_with_3cn = 0
        
    output_str = "Iter: " + str(iter) + \
                 " Stage: " + str(state_index) + \
                 " UnF: " + '{:> 5.2f}'.format(temp_unbalanced_force) + \
                 " Cn: " + '{:> 7.4f}'.format(temp_sphere_data[3, :].mean()) + \
                 " Y_pi " + '{:> 7.4f}'.format(temp_pile_position_y) + \
                 " MaxY_sp: " + '{:> 7.4f}'.format(temp_sphere_data[2, :].max()) + \
                 " MaxY_sp_cn3: " + '{:> 7.4f}'.format(temp_maxY_with_3cn) + \
                 " Fx_pi: " + '{:> 8.2f}'.format(temp_pile_force_x) + \
                 " Fy_pi: " + '{:> 8.2f}'.format(temp_pile_force_y) + \
                 " Fz_pi: " + '{:> 8.2f}'.format(temp_pile_force_z)
    
    print(output_str)
    
    # export some global values to csv file
    output_values = [iter, state_index, temp_unbalanced_force, 
                     temp_friction_angle, temp_sphere_data[3, :].mean(),
                     temp_pile_position_y, temp_sphere_data[2, :].max(),
                     temp_maxY_with_3cn,
                     temp_pile_force_x, 
                     temp_pile_force_y,
                     temp_pile_force_z] + temp_voxel_porosity

    output_values_str = ""

    for temp in output_values:
        output_values_str += str(temp) + ","
        
    output_values_str = output_values_str[:-1] + "\n"
    with open(output_file_path, 'a') as f:
        f.write(output_values_str)
        
    plot.addData(i=iter,
                 Cn=temp_sphere_data[3, :].mean(),
                 maxY_sp=temp_sphere_data[2, :].max(),
                 maxY_sp_3cn=temp_maxY_with_3cn,
                 Pos=temp_pile_position_y,
                 Fy_pi=temp_pile_force_y,
                 VoxPo_1=temp_voxel_porosity[0],
                 VoxPo_2=temp_voxel_porosity[int(len(temp_voxel_porosity)/3)],
                 VoxPo_3=temp_voxel_porosity[int(len(temp_voxel_porosity)/3*2)])
    
    # export VTK file
    if initial_parameters["flag_output_VTK"] and state_index >= 2:
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
    global state_index, temp_initial_pile_disp, sphere_id, temp_prev_stage_iter

    # initial gravity deposition
    if state_index == 0:
        if O.iter % 100 == 0:
            temp_sphere_vel = [O.bodies[i].state.vel[1] for i in sphere_id]
            temp_sphere_Y = [O.bodies[i].state.pos[1] for i in sphere_id]
            temp_sphere_coord_num = [len(O.interactions.withBody(i)) for i in sphere_id]
            temp_sphere_data = np.array([sphere_id, temp_sphere_vel, temp_sphere_Y, temp_sphere_coord_num])
            
            if len(temp_sphere_data[2, temp_sphere_data[3, :]>=3]) != 0:
                temp_maxY_with_3cn = temp_sphere_data[2, temp_sphere_data[3, :]>=3].max()
            else:
                temp_maxY_with_3cn = 0
            
            temp_sphere_id_positive_vel = temp_sphere_data[0, temp_sphere_data[1, :] > 0]
            
            for temp_sphere_id_positive_vel_each in temp_sphere_id_positive_vel:
                O.bodies[int(temp_sphere_id_positive_vel_each)].state.vel[1] = 0
                        
            if temp_sphere_data[2, :].max() == temp_maxY_with_3cn and temp_sphere_data[3, :].mean() > 3:
                export.text(str(temp_sp_file_path))
                
                temp_num_layers = int(math.ceil(initial_parameters["sphere_pack_target_Y"] / temp_sphere_data[2, :].max()))
                for i in range(temp_num_layers):
                    temp_additional_spheres = ymport.text(str(temp_sp_file_path), shift=Vector3(0, temp_sphere_data[2, :].max()*(i+1), 0), material=material_sphere_Id)
                    temp_additional_spheres_id = O.bodies.append(temp_additional_spheres)
                    sphere_id.extend(temp_additional_spheres_id)
                    
                print("Some spheres are added. now", len(sphere_id), "spheres exist")
                    
                
                temp_prev_stage_iter = O.iter
                state_index = 1
            """
            print(O.iter, state_index,
                  '{:> 7.5f}'.format(temp_sphere_data[2, :].max()), 
                  '{:> 7.5f}'.format(temp_sphere_data[3, :].mean()),
                  '{:> 7.5f}'.format(temp_sphere_data[1, :].mean()),
                  '{:> 7.5f}'.format(temp_maxY_with_3cn))
            """
            

    # copy particle assemblies to make the layer height equal to the target height
    elif state_index == 1:
        if O.iter % 100 == 0 and O.iter - temp_prev_stage_iter > 10000:
            temp_sphere_vel = [O.bodies[i].state.vel[1] for i in sphere_id]
            temp_sphere_Y = [O.bodies[i].state.pos[1] for i in sphere_id]
            temp_sphere_coord_num = [len(O.interactions.withBody(i)) for i in sphere_id]
            temp_sphere_data = np.array([sphere_id, temp_sphere_vel, temp_sphere_Y, temp_sphere_coord_num])
            
            temp_sphere_id_positive_vel = temp_sphere_data[0, temp_sphere_data[1, :] > 0]
            
            for temp_sphere_id_positive_vel_each in temp_sphere_id_positive_vel:
                O.bodies[int(temp_sphere_id_positive_vel_each)].state.vel[1] = 0
            
            if len(temp_sphere_data[2, temp_sphere_data[3, :]>=3]) != 0:
                temp_maxY_with_3cn = temp_sphere_data[2, temp_sphere_data[3, :]>=3].max()
            else:
                temp_maxY_with_3cn = 0
                
            if temp_maxY_with_3cn > initial_parameters["sphere_pack_target_Y"] and temp_sphere_data[3, :].mean() > 3:
                temp_unused_sphere_id = temp_sphere_data[0, temp_sphere_data[2, :] > initial_parameters["sphere_pack_target_Y"]]
                for temp_unused_sphere_id_each in temp_unused_sphere_id:
                     O.bodies.erase(int(temp_unused_sphere_id_each))
                     sphere_id.remove(int(temp_unused_sphere_id_each))
                
                print(len(temp_unused_sphere_id), " spheres are deleted")
                print("Now", len(sphere_id), "spheres exist")
                
                temp_prev_stage_iter = O.iter
                state_index = 2
            
            """
            print(O.iter, state_index,
                  '{:> 7.5f}'.format(temp_sphere_data[2, :].max()), 
                  '{:> 7.5f}'.format(temp_sphere_data[3, :].mean()),
                  '{:> 7.5f}'.format(temp_sphere_data[1, :].mean()),
                  '{:> 7.5f}'.format(temp_maxY_with_3cn))
            """
        
    
    elif state_index == 2:
        temp_sphere_Y = np.array([O.bodies[i].state.pos[1] for i in sphere_id])
        
        
        for pile_facet_id in pile_facet_data[:, 0]:
            O.bodies[int(pile_facet_id)].bounded = True
            O.bodies[int(pile_facet_id)].state.blockedDOFs = "xyzXYZ"
            O.bodies[int(pile_facet_id)].state.vel = Vector3(0, initial_parameters["pile_insertion_velocity"], 0)
                                    
        state_index = 3

    elif state_index == 3:
        
        if O.iter % 100 == 0:
            temp_sphere_Y = np.array([O.bodies[i].state.pos[1] for i in sphere_id])
            
            flag_bottom_reached = O.bodies[int(pile_facets_id_bottom[0])].state.pos[1] <= temp_max_sphere_size
            flag_pile_length = (temp_sphere_Y.max() - O.bodies[int(pile_facets_id_top[0])].state.pos[1]) >= 0
        
            if flag_bottom_reached or flag_pile_length:
                O.pause()


plot.plots = {"i": ("maxY_sp", "maxY_sp_3cn"),
              " i": ("VoxPo_1","VoxPo_2","VoxPo_3"),
              "i ": ("Cn"),
              "Fy_pi": ("Pos")}

plot.plot()