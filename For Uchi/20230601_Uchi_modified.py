"""
以下は下のプログラムコードを作るにあたって書いた大雑把なメモ書きです。
わからないところがあれば答えられる範囲で答えるので、なんでも聞いてください。
https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Peri3dController

# 基本的な内容
Yadeには境界を制御するためのクラスが用意されているようです。このクラスにはいくつか種類があるのですが、今回はPeri3dControllerというものを使いたいと思います。

"""

from __future__ import print_function
from yade import pack, plot
from yade.gridpfacet import *
import numpy as np

# ターミナルの出力が詰まっていて見づらいための区切り線
print("")
print("-------------------------------------")

frictangle = 35 / 180 * np.pi    # 35度をラジアンに変換
density = 3700.0                 # 密度大きすぎないですか？
young = 3e8
iso_consolidation_stress = 100e3 # 等方圧密時の目標鉛直応力：100kPa=100×10^3Pa
shear_stress_amplitude = 20e3    # 繰り返しせん断時の目標最大せん断振幅：20kPa=
shear_stress_torrelance = 0.5e3  # 繰り返しせん断時のせん断振幅のばらつき：0.5kPa
flag_cyclic_loading = False      # 繰り返し載荷を始めるかどうか(圧密が終了したかどうか)のフラグ
flag_cyclic_forward = True       # 繰り返し載荷の方向を決めるフラグ
target_num_cycle = 10            # 目標の繰り返し回数
current_num_cycle = 0            # 現時点での繰り返し回数

# CohFrictMatクラスはどうやら粘着力を考慮できる材料のクラスみたいですね。ただ
# https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.CohFrictMat
# を見ると、normalCohesionとshearCohesionに正の値を設定しないと、粘着力は考慮しないみたいです。
mat_sp = CohFrictMat(alphaKr=0.5,
                     alphaKtw=0.5,
                     young=young,
                     poisson=0.33,
                     frictionAngle=frictangle,
                     density=density)
O.materials.append(mat_sp)

sp = pack.randomPeriPack(radius=0.10, initSize=Vector3(6, 6, 2), memoizeDb='/tmp/packDb.sqlite')
sp.toSimulation(color=(0, 0, 1))  # pure blue

# 下のコードはどのような意味がある？このコードを有効化すると粒子同士の内部摩擦角が0になってしまうけど...
# O.materials[0].frictionAngle = 0

# 等方圧密用の境界制御のためのクラスです。
Peri3D_iso = Peri3dController(
    # Goalでは応力テンソルかひずみテンソルの対角成分を含む上半分の目標数値を決めます。
    # 順番は00,11,22,12,02,01の順です。
    # まずは最初等方的に圧密をしたいので、(iso_consolidation_stress, iso_consolidation_stress, iso_consolidation_stress, 0, 0, 0)とします。
    
    goal = (-iso_consolidation_stress,
            -iso_consolidation_stress, 
            -iso_consolidation_stress,
            0, 0, 0),
    
    # 説明が少し難しいんですが、順番にします。まず0bというのはpythonでは0bのあとに続く数字は2進数だよという目印です。
    # 問題は111000の部分なんですが6桁の数字のそれぞれが上のテンソルの位置に対応します。
    # 例えば最初の1は、テンソルの00成分に対応します。
    # そして、それぞれの数字は0か1の値を取ります。
    # 0の場合はそのテンソルの目標値がひずみであること、1の場合は応力であることを表します。
    # なので下の例だと、すべての成分の目標値は応力であることを示します。
    stressMask = 0b111111,
    
    # このステップ数が終了するとdoneHookに指定した関数が呼ばれます。
    nSteps = 2000,
    doneHook = "cyclic_shear()",
    
    # 最大のひずみ速度を決めます。大きすぎると負のダイレタンシーが発生した際に応力を下げすぎてしまうので、一旦かなり小さい値にして、繰り返し載荷で発生するひずみ振幅よりも小さく(例えば1/10とか)にすればいいのではないでしょうか？
    maxStrainRate = 0.5
)

# 繰り返し載荷(せん断力を増やす際)の際の境界制御のクラスです。
Peri3D_cyclic_forward = Peri3dController(
    goal = (-iso_consolidation_stress,
            -iso_consolidation_stress, 
            -iso_consolidation_stress,
            0, shear_stress_amplitude, 0),
    stressMask = 0b111111,
    
    nSteps = 2000,
    doneHook = "cyclic_shear()",
    maxStrainRate = 0.05,
    xxPath = ((0, 1), (1, 1)),
    yyPath = ((0, 1), (1, 1)),
    zzPath = ((0, 1), (1, 1))
)

# 繰り返し載荷(せん断力を減らす際際)の際の境界制御のクラスです。
Peri3D_cyclic_backward = Peri3dController(
    goal = (-iso_consolidation_stress,
            -iso_consolidation_stress, 
            -iso_consolidation_stress,
            0, -shear_stress_amplitude, 0),
    stressMask = 0b111111,
    
    nSteps = 2000,
    doneHook = "cyclic_shear()",
    maxStrainRate = 0.05,
    xxPath = ((0, 1), (1, 1)),
    yyPath = ((0, 1), (1, 1)),
    zzPath = ((0, 1), (1, 1))
    )

# エンジンと呼ばれる部分です。
# 詳しくは割愛します。
O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Facet_Aabb(), Bo1_Box_Aabb()]),
    InteractionLoop(
        # interaction loop
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(),
         Ig2_Box_Sphere_ScGeom()],
        [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()],
        [Law2_ScGeom6D_CohFrictPhys_CohesionMoment(
            always_use_moment_law=True), Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    Peri3D_iso,
    NewtonIntegrator(damping=0.2),
    PyRunner(realPeriod=0.5, command="addPlotData()")
]

O.dt = .1 * PWaveTimeStep()

# 繰り返しせん断用の関数です。
def cyclic_shear():
    global flag_cyclic_loading, flag_cyclic_forward
    global iso_consolidation_stress, shear_stress_amplitude, current_num_cycle, target_num_cycle, shear_stress_torrelance
    
    e00, e11, e22, e12, e02, e01 = O.engines[3].strain
    s00, s11, s22, s12, s02, s01 = O.engines[3].stress
    
    if not flag_cyclic_loading:
        print("Consolidation has finished! Now proceeding to cyclic loading")
        flag_cyclic_loading = True
        
        O.engines = O.engines[0:3] + [Peri3D_cyclic_forward] + O.engines[4:]
        O.engines[3].strain = (e00, e11, e22, e12, e02, e01)
        O.engines[3].stressIdeal = (s00, s11, s22, s12, s02, s01)
        O.engines[3].stressRate = (0, 0, 0, 0, 0, 0)

        O.engines[3].progress = 0

        
        
    else:
        if (current_num_cycle > target_num_cycle):
            
            print("Cyclic loading has just finished!")
            finish_simulation()
            
        elif (current_num_cycle <= target_num_cycle) and flag_cyclic_forward:
            
            flag_under_torrelance = (shear_stress_amplitude - s02) < shear_stress_torrelance           
           
            if flag_under_torrelance:
                print("Current Cycle: " + str(current_num_cycle)) 
                
                flag_cyclic_forward = False
                O.engines = O.engines[0:3] + [Peri3D_cyclic_backward] + O.engines[4:]
                O.engines[3].strain = (e00, e11, e22, e12, e02, e01)
                O.engines[3].stressIdeal = (s00, s11, s22, s12, s02, s01)
                O.engines[3].stressRate = (0, 0, 0, 0, 0, 0)
                O.engines[3].progress = 0
                
                
                current_num_cycle += 0.5

                
        elif (current_num_cycle <= target_num_cycle) and not flag_cyclic_forward:
                        
            flag_under_torrelance = (-shear_stress_amplitude - s02) > -shear_stress_torrelance
            
            if flag_under_torrelance: 
                print("Current Cycle: " + str(current_num_cycle))
                
                flag_cyclic_forward = True
                O.engines = O.engines[0:3] + [Peri3D_cyclic_forward] + O.engines[4:]
                O.engines[3].strain = (e00, e11, e22, e12, e02, e01)
                O.engines[3].stressIdeal = (s00, s11, s22, s12, s02, s01)
                O.engines[3].stressRate = (0, 0, 0, 0, 0, 0)
                O.engines[3].progress = 0
                
                current_num_cycle += 0.5
        
# 載荷が終了した際に呼ばれる関数です。
def finish_simulation():
    O.pause()
    
# 図化のための関数です。
def addPlotData():
    i = O.iter
    s00, s11, s22, s12, s02, s01 = O.engines[3].stress / 1000
    e00, e11, e22, e12, e02, e01 = O.engines[3].strain
    gs00, gs11, gs22, gs12, gs02, gs01 = O.engines[3].stressGoal / 1000
    
    # temp = 2 * (O.engines[3].stress - O.engines[3].stressIdeal) - (O.engines[3].stressOld - (O.engines[3].stressIdeal - O.engines[3].stressRate * O.dt))
    print(O.engines[3].stressIdeal, O.engines[3].stressRate)

    print('pro: {: .2f}'.format(O.engines[3].progress),
          ' s00: {: .2f}'.format(s00),
          ' s11: {: .2f}'.format(s11),
          ' s22: {: .2f}'.format(s22),
          ' s12: {: .2f}'.format(s12),
          ' s02: {: .2f}'.format(s02),
          ' s01: {: .2f}'.format(s01))
    
    plot.addData(i = i,
                 s00 = s00,
                 s22 = s22,
                 s02 = s02,
                 e00 = e00,
                 e22 = e22,
                 e02 = e02)

plot.plots = {"i": ("s00", "s22", "s02"),
              "i ": ("e00", "e22", "e02"),
              " e02": ("s02"),
              " s22 ": ("s02")}

plot.plot()
