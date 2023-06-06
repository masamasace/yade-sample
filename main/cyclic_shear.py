"""
以下は下のプログラムコードを作るにあたって書いた大雑把なメモ書きです。
わからないところがあれば答えられる範囲で答えるので、なんでも聞いてください。
https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Peri3dController

# 基本的な内容
Yadeには境界を制御するためのクラスが用意されているようです。
このクラスにはいくつか種類があるのですが、今回はPeri3dControllerというものを使いたいと思います。

# TODO あるいは要確認事項
- 基本的な物性値が結果に与える影響、論文にある？
    - 内部摩擦角、粒子単体のヤング率、接触モデルの各パラメータ
    - 豊浦砂や珪砂といった実際の粒径加積曲線を用いる (yade._packSpheres.SpherePack.makeCloudの引数としてある)
    - 粒子数
- 等方圧密をしているのにせん断ひずみが出る
- 重力堆積の方向が結果に与える影響
- 重力堆積の際に毎回生成される際の粒子にばらつきを与える(Seed値をいじる)
"""

from __future__ import print_function
from yade import pack, plot
from yade.gridpfacet import *
import numpy as np
import datetime
from pathlib import Path
import os
import sys

# ターミナルの出力が詰まっていて見づらいための区切り線
print("")
print("-------------------------------------")


############## 定数・変数の定義 (単位は全てSI単位系(長さはm, 時間はsec, 質量はkg)) ##############

# 粒子と粒子間の接触モデルに関する定数・変数
frictangle = 35 / 180 * np.pi           # (要検討)内部摩擦角、35度をラジアンに変換しています。
density = 2650.0                        # 粒子単体の密度 (kg/m3)です。
young = 3e8                             # (要検討)粒子単体のヤング率 (kPa)です。
maxCorner_x = 6                         # 直方体の境界の一番上のx座標
maxCorner_y = 6                         # 直方体の境界の一番上のy座標
maxCorner_z = 2                         # 直方体の境界の一番上のz座標      
rMean = 0.1                             # 粒子の平均粒径
rRelFuzz = 0                            # 粒子の粒径の分散値
particle_num = 50                     # 粒子数
voidRatio = 0.5                         # 重力堆積終了後の目標間隙比(圧密終了後ではないことに注意)
porosity = voidRatio / (1 + voidRatio)  # 空隙率(間隙比から計算)。
print('\nTarget Porosity: {:.3f}'.format(porosity))

# 載荷プロセス全体に関する定数・変数
state_index = 0                         # 0:重力堆積、1:圧密、2:繰り返しせん断(Forward)、3:繰り返しせん断(Backward)の状態管理をするための変数です。
gravity_vector = [0, -9.81, 0]          # 重力ベクトル(従来はz軸負の方向に作用させるが、粒子の異方性に与える影響が無視できない(要確認)ため、今回はy軸負の方向に作用させる)

# 圧密に関する定数・変数
consolidation_stress = 100e3            # 等方圧密時の目標鉛直応力です。：100kPa=100×10^3Pa
consolidation_nSteps = 2000             # 圧密時での最低限必要なステップ数です。
consolidation_max_strain_rate = 5       # 圧密時での最大ひずみ速度です。

# 繰り返しせん断に関する定数・変数
flag_constant_pressure = False           # 繰り返しせん断時に定圧条件で行うか定体積条件で行うかを決めます。Trueが定圧条件です。
flag_stress_threshold = True           # 繰り返しせん断時の反転を応力で定義するか、ひずみで定義するかを決めます。Trueが応力定義です。
cyclic_loading_max_strain_rate = 0.005   # 繰り返しせん断時での最大ひずみ速度です。(反転条件がひずみの場合には後で更新されます。)
cyclic_shear_stress_amplitude = 10e3    # 繰り返しせん断時の目標最大せん断応力振幅です。：20kPa
cyclic_shear_strain_amplitude = 0.0005    # 繰り返しせん断時の目標最大せん断ひずみ振幅です。：0.05→5%
target_num_cycle = 2                    # 目標の繰り返し回数です。
current_num_cycle = 0                   # 現時点での繰り返し回数です。
cyclic_loading_nSteps = 10000            # 繰り返しせん断時での最低限必要なステップ数です。

# 出力に関する定数・変数
flag_header_done = False                # 出力ファイルにヘッダーが生成されたかどうかのフラグ

# 結果を保存するためのファイルとフォルダを作成します。
dt_start = datetime.datetime.now()
print("Start simulation started at " + dt_start.strftime('%Y/%m/%d %H:%M:%S.%f'))

output_folder_path = Path(os.path.abspath(os.path.dirname(sys.argv[0]))).parent / "result"
output_folder_path.mkdir(exist_ok=True)
output_file_path = output_folder_path / (dt_start.strftime('%Y-%m-%d_%H-%M-%S') + "_output.csv")


############## 粒子と接触モデルの定義 ##############
# 2つの粒子の間に働く力関係を計算するためのモデルを定義します。 
# CohFrictMatクラスは粘着力を考慮できる材料のクラスです。
mat_sp = CohFrictMat(alphaKr=0.5,
                     alphaKtw=0.5,
                     young=young,
                     poisson=0.33,
                     frictionAngle=frictangle,
                     density=density)
O.materials.append(mat_sp)

# 境界と粒子の生成
O.periodic = True
O.cell.setBox((maxCorner_x, maxCorner_y, maxCorner_z))

# 重力堆積の際に粒子を生成する領域
# seed値は呼び出されるごとに変更するようにしたほうがいい
spheres_set = pack.SpherePack().makeCloud(Vector3(0, maxCorner_y*0.8, 0),
                                          Vector3(maxCorner_x, maxCorner_y, maxCorner_z),
                                          rMean = rMean,
                                          rRelFuzz = rRelFuzz,
                                          num = particle_num,
                                          periodic = True,
                                          seed = 1)
O.bodies.append(spheres_set)

# 解析の際の時間ステップを決めるコードです。
# MEMO: 本来であればすべての2粒子間接触から求まる固有周期の最小値を時間ステップとして採用するべきなんですが、初期状態では接触点の数が0の場合もあるので、近似的に粒子のP波速度の0.1倍を使っています。(Class Referenceと計算式が違うがほとんどのサンプルコードがこれぐらいの値を採用している...)
O.dt = .1 * PWaveTimeStep()

# 繰り返しせん断のときの最大ひずみ速度、どうやら定圧でひずみ反転条件の場合には、O.dt、cyclic_loading_nSteps、cyclic_shear_strain_amplitudeの3つから計算される値を用いたほうがよいらしい
if not flag_stress_threshold and flag_constant_pressure:
    print("Strain Rate is updated! :" + '{:.4g}'.format(cyclic_loading_max_strain_rate), end="")
    cyclic_loading_max_strain_rate = cyclic_shear_strain_amplitude / (O.dt * cyclic_loading_nSteps)
    print(" -> " + '{:.4g}'.format(cyclic_loading_max_strain_rate))
    

############## 境界制御モデルの定義 ##############
# DEMのシミュレーションではプログラムで作った仮想的な板(境界)を動かすことで、シミュレーションを実行します。
# このときに境界をどのように制御するかで用いるクラスが異なってきます。
# 今回使用するのは(調べたところ)最も細かい制御ができるであろうPeri3dControllerを用います
# 今回のシミュレーションでは、圧密→繰り返しせん断という2つのステップが存在します。
# また繰り返しせん断にはせん断応力が増加する場合と、減少する場合の2種類があるため、
# 合計3つの境界制御用のインスタンスを作成します

# 等方圧密用の境界制御のためのクラスです。
Peri3D_iso = Peri3dController(
    # Goalでは応力テンソル、あるいはひずみテンソルの対角成分を含む上半分の目標数値を決めます。
    # 順番は00,11,22,12,02,01の順です。
    # まずは最初等方的に圧密をしたいので、(consolidation_stress, consolidation_stress, consolidation_stress, 0, 0, 0)とします。
    
    goal = (-consolidation_stress,
            -consolidation_stress, 
            -consolidation_stress,
            0, 0, 0),
    
    # このstressMaskは説明が少し難しいんですが、順番にします。
    # まず0bというのはpythonでは0bのあとに続く数字は2進数だよという目印です。
    # 問題は111111の部分なんですが6桁の数字のそれぞれが上のテンソルの位置に対応します。
    # 注意点として"最下位ビット"から00, 11, 22, 12, 02, 01の順で並んでいるので以下のような対応関係になります。
    # stressMask = 0b    1    1    1    1    1    1
    #                   01   02   12   22   11   00
    # そして、それぞれの数字は0か1の値を取ります。
    # 0の場合はそのテンソルの目標値がひずみであること、1の場合は応力であることを表します。
    # なので上の例だと、すべての成分の目標値は応力であることを示します。
    stressMask = 0b111111,
    
    # このステップ数が終了するとdoneHookに指定した関数が呼ばれます。
    nSteps = consolidation_nSteps,
    doneHook = "checkState()",
    
    # 最大のひずみ速度の値を代入します。
    maxStrainRate = consolidation_max_strain_rate
)

# 繰り返しせん断の反転条件が応力の場合と、ひずみの場合でこの後のインスタンス生成時に代入するパラメータが変わってきます。
# なのでここでif文を使って条件整理をします。
if flag_constant_pressure:
    if flag_stress_threshold:
        goal_forward_temp = (-consolidation_stress,
                            -consolidation_stress, 
                            -consolidation_stress,
                            0, cyclic_shear_stress_amplitude, 0)
        goal_backward_temp = (-consolidation_stress,
                            -consolidation_stress, 
                            -consolidation_stress,
                            0, -cyclic_shear_stress_amplitude, 0)
        stressMask_temp = 0b111111
    else:
        goal_forward_temp = (-consolidation_stress,
                            -consolidation_stress, 
                            -consolidation_stress,
                            0, cyclic_shear_strain_amplitude, 0)
        goal_backward_temp = (-consolidation_stress,
                            -consolidation_stress, 
                            -consolidation_stress,
                            0, -cyclic_shear_strain_amplitude, 0)
        
        # ここでstressMaskの説明をもう一回します。今回の場合の対応関係は以下のとおりです。
        # goalの成分の順番と、逆になっている点に注意しましょう。
        # stressMask = 0b    1    0    1    1    1    1
        #                   01   02   12   22   11   00
        # 0の場合はそのテンソルの目標値がひずみであること、1の場合は応力であることを表すので、
        # 今回の場合は02成分(xz成分)のみが目標値がひずみであり、他の5成分は全て応力であることを示します。
        stressMask_temp = 0b101111
else:
    # 定体積条件では圧密後のひずみ値を記録する必要があるので、goalの5番目の成分以外は暫定値です。
    if flag_stress_threshold:
        goal_forward_temp = (0, 0, 0, 0, cyclic_shear_stress_amplitude, 0)
        goal_backward_temp = (0, 0, 0, 0, -cyclic_shear_stress_amplitude, 0)
        stressMask_temp = 0b010000
    else:
        goal_forward_temp = (0, 0, 0, 0, cyclic_shear_strain_amplitude, 0)
        goal_backward_temp = (0, 0, 0, 0, -cyclic_shear_strain_amplitude, 0)
        stressMask_temp = 0b000000

# 繰り返し載荷(せん断力を増やす際)の際の境界制御のクラスです。
# xxPath = ((0, 1), (1, 1))の部分について。これがないと載荷初期段階の応力やひずみが0スタートになってしまい、挙動が不自然になります。
Peri3D_cyclic_forward = Peri3dController(
    goal = goal_forward_temp,
    stressMask = stressMask_temp,
    nSteps = cyclic_loading_nSteps,
    doneHook = "checkState()",
    maxStrainRate = cyclic_loading_max_strain_rate,
    xxPath = ((0, 1), (1, 1)),
    yyPath = ((0, 1), (1, 1)),
    zzPath = ((0, 1), (1, 1)),
)

# 繰り返し載荷(せん断力を減らす際際)の際の境界制御のクラスです。
Peri3D_cyclic_backward = Peri3dController(
    goal = goal_backward_temp,
    stressMask = stressMask_temp,
    nSteps = cyclic_loading_nSteps,
    doneHook = "checkState()",
    maxStrainRate = cyclic_loading_max_strain_rate,
    xxPath = ((0, 1), (1, 1)),
    yyPath = ((0, 1), (1, 1)),
    zzPath = ((0, 1), (1, 1)),
)


############## エンジンの定義 ##############
# エンジンと呼ばれる部分です。
# ちょっと現段階では自分も完全に理解できていないので、説明は割愛させてください。
O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Facet_Aabb(), Bo1_Box_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(),
         Ig2_Box_Sphere_ScGeom()],
        [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()],
        [Law2_ScGeom6D_CohFrictPhys_CohesionMoment(
            always_use_moment_law=True), Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(damping=0.2, gravity=gravity_vector),
    PyRunner(iterPeriod=1, command="checkState()"),
    # PyRunner(realPeriod=0.5, command="addPlotData()")
]

# エネルギーを追跡します
O.trackEnergy = True

# 状態遷移用の関数です
def checkState():
    global state_index, porosity
    global consolidation_stress, cyclic_shear_stress_amplitude
    global current_num_cycle, target_num_cycle
    
    # 重力堆積の状態確認をします。
    if state_index == 0:
        print(unbalancedForce(), len(O.bodies), utils.porosity())
        if (unbalancedForce() < 0.05 or unbalancedForce() == None) and utils.porosity() > porosity:
            O.bodies.append(spheres_set)
        O.pause()
        
    
    elif state_index == 1:
        e00, e11, e22, e12, e02, e01 = O.engines[3].strain
        s00, s11, s22, s12, s02, s01 = O.engines[3].stress
    
        print("Consolidation has finished! Now proceeding to cyclic loading")
        state_index += 1
        print(len(O.bodies), utils.porosity())
        
        O.engines = O.engines[0:3] + [Peri3D_cyclic_forward] + O.engines[4:]
        O.engines[3].strain = (e00, e11, e22, e12, e02, e01)
        O.engines[3].stress = (s00, s11, s22, s12, s02, s01)
        O.engines[3].stressIdeal = (s00, s11, s22, s12, s02, s01)
        O.engines[3].stressRate = (0, 0, 0, 0, 0, 0)
        
        if not flag_constant_pressure:
            if flag_stress_threshold:
                O.engines[3].goal = (e00, e11, e22, e12, cyclic_shear_stress_amplitude, e01)
            else:
                O.engines[3].goal = (e00, e11, e22, e12, cyclic_shear_strain_amplitude, e01)
        O.engines[3].progress = 0

    else:
        e00, e11, e22, e12, e02, e01 = O.engines[3].strain
        s00, s11, s22, s12, s02, s01 = O.engines[3].stress
        
        if (current_num_cycle > target_num_cycle):
            print("Cyclic loading has just finished!")
            finish_simulation()
            
        elif (current_num_cycle <= target_num_cycle) and (state_index == 2):
            
            if flag_stress_threshold:
                flag_under_torrelance = s02 > cyclic_shear_stress_amplitude
            else:
                flag_under_torrelance = e02 > cyclic_shear_strain_amplitude
                
            if flag_under_torrelance:
                print("Current Cycle: " + str(current_num_cycle), end="") 
                
                state_index = 3
                
                O.engines = O.engines[0:3] + [Peri3D_cyclic_backward] + O.engines[4:]
                O.engines[3].strain = (e00, e11, e22, e12, e02, e01)
                O.engines[3].stress = (s00, s11, s22, s12, s02, s01)
                O.engines[3].stressIdeal = (s00, s11, s22, s12, s02, s01)
                O.engines[3].stressRate = (0, 0, 0, 0, 0, 0)
                if not flag_constant_pressure:
                    if flag_stress_threshold:
                        O.engines[3].goal = (e00, e11, e22, e12, -cyclic_shear_stress_amplitude, e01)   
                    else:
                        O.engines[3].goal = (e00, e11, e22, e12, -cyclic_shear_strain_amplitude, e01)   
                O.engines[3].progress = 0
        
                current_num_cycle += 0.5
                print(" Next Cycle: " + str(current_num_cycle)) 

        elif (current_num_cycle <= target_num_cycle) and (state_index == 3):
                        
            if flag_stress_threshold:
                flag_under_torrelance = s02 < -cyclic_shear_stress_amplitude
            else:
                flag_under_torrelance = e02 < -cyclic_shear_strain_amplitude
                
            if flag_under_torrelance: 
                print("Current Cycle: " + str(current_num_cycle), end="") 

                state_index = 2

                O.engines = O.engines[0:3] + [Peri3D_cyclic_forward] + O.engines[4:]
                O.engines[3].strain = (e00, e11, e22, e12, e02, e01)
                O.engines[3].stress = (s00, s11, s22, s12, s02, s01)
                O.engines[3].stressIdeal = (s00, s11, s22, s12, s02, s01)
                O.engines[3].stressRate = (0, 0, 0, 0, 0, 0)
                if not flag_constant_pressure:
                    if flag_stress_threshold:
                        O.engines[3].goal = (e00, e11, e22, e12, cyclic_shear_stress_amplitude, e01)
                    else:
                        O.engines[3].goal = (e00, e11, e22, e12, cyclic_shear_strain_amplitude, e01)
                O.engines[3].progress = 0
                
                current_num_cycle += 0.5
                print(" Next Cycle: " + str(current_num_cycle)) 
                
        
# 載荷が終了した際に呼ばれる関数です。
def finish_simulation():
    O.pause()
    
    # エネルギーの項が全て入っていない場合があるため、もう一度ファイルを読み込んでヘッダーだけ修正します。
    with open(output_file_path, 'r') as f:
        lines = f.readlines()
    
    key_list = "step,s00,s11,s22,s12,s02,s01,e00,e11,e22,e12,e02,e01"
    energy_dict = dict(O.energy.items())
    for temp in list(energy_dict.keys()):
        key_list += "," + temp
    key_list += "\n"
    
    lines[0] = key_list
    
    with open(output_file_path, 'w') as f:
        f.writelines(lines)
    
    
# 図化のための関数です。
def addPlotData():
    global flag_header_done, previous_energy_header
    
    i = O.iter
    s00, s11, s22, s12, s02, s01 = O.engines[3].stress / 1000
    e00, e11, e22, e12, e02, e01 = O.engines[3].strain
    
    # 平均応力 (圧縮が正に変換)
    mean_stress = -(s00 + s11 + s22) / 3
    # ミーゼスの相当応力(偏差応力テンソルの第2次不変量J2の平方根を√3倍したもの)
    deviatoric_stress = ((((s00 - s11) ** 2 + (s11 - s22) ** 2 + (s22 - s00) ** 2) / 6 + s01 ** 2 + s12 ** 2 + s02 ** 2) * 3) ** (1 / 2)
    
    print('progress: {: .5f}'.format(O.engines[3].progress),
          ' p:{:>8.3f}'.format(mean_stress),
          ' q:{:>8.3f}'.format(deviatoric_stress),
          ' s02:{:>8.3f}'.format(s02),
          ' s22:{:>8.2f}'.format(s22),
          ' e02:{:>8.3f}'.format(e02 * 100),  # ひずみは%表記
          ' e22:{:>6.2f}'.format(e22 * 100))
    
    # print(O.engines[3].stressIdeal, O.engines[3].strainRate)
    
    plot.addData(i = i,
                 s00 = s00,
                 s22 = s22,
                 s02 = s02,
                 e00 = e00,
                 e22 = e22,
                 e02 = e02)
    
    # 出力データを準備します
    energy_dict = dict(O.energy.items())
        
    if not flag_header_done:
        flag_header_done = True
        
        # エネルギー項を含めた出力ファイルのヘッダーを用意します。
        key_list = "step,s00,s11,s22,s12,s02,s01,e00,e11,e22,e12,e02,e01"
        energy_dict = dict(O.energy.items())
        for temp in list(energy_dict.keys()):
            key_list += "," + temp
        key_list += "\n"
        with open(output_file_path, 'w') as f:
            f.write(key_list)
    
    energy_values = list(energy_dict.values())
    output_values = [i,s00,s11,s22,s12,s02,s01,e00,e11,e22,e12,e02,e01] + energy_values
    output_values_str = ""
    
    for temp in output_values:
        output_values_str += str(temp) + "," 
        
    output_values_str = output_values_str[:-1] + "\n"
    with open(output_file_path, 'a') as f:
        f.write(output_values_str)

plot.plots = {"i": ("s00", "s22", "s02"),
              "i ": ("e00", "e22", "e02"),
              " e02": ("s02"),
              " s22 ": ("s02")}

plot.plot()
