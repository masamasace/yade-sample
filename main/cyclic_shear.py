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

############## 定数の定義 (単位は全てSI単位系(長さはm, 時間はsec, 質量はkg)) ##############

frictangle = 35 / 180 * np.pi           # (要検討)内部摩擦角、35度をラジアンに変換しています。
density = 2650.0                        # 土粒子単体の密度 (kg/m3)です。
young = 3e8                             # (要検討)粒子単体のヤング率 (kPa)です。
consolidation_stress = 100e3            # 等方圧密時の目標鉛直応力です。：100kPa=100×10^3Pa
cyclic_shear_stress_amplitude = 20e3    # 繰り返しせん断時の目標最大せん断応力振幅です。：20kPa
cyclic_shear_strain_amplitude = 0.05     # 繰り返しせん断時の目標最大せん断ひずみ振幅です。：0.05→5%
shear_stress_torrelance = 0.5e3         # (要修正)繰り返しせん断時のせん断振幅のばらつきです。：0.5kPa
flag_cyclic_loading = False             # 繰り返しせん断を始めるかどうか(圧密が終了したかどうか)のフラグです。
flag_cyclic_forward = True              # 繰り返しせん断の方向を決めるフラグです。
flag_stress_threshold = True            # 繰り返しせん断時の反転を応力で定義するか、ひずみで定義するかを決めます。Trueが応力定義です。
target_num_cycle = 10                   # 目標の繰り返し回数です。
current_num_cycle = 0                   # 現時点での繰り返し回数です。
consolidation_nSteps = 1000             # 圧密時での最低限必要なステップ数です。
cyclic_loading_nSteps = 2000            # 繰り返しせん断時での最低限必要なステップ数です。
consolidation_max_strain_rate = 0.5     # 圧密時での最大ひずみ速度です。
cyclic_loading_max_strain_rate = 0.5  # 圧密時での最大ひずみ速度です。
state_index = 0

############## 接触モデルの定義 ##############
# 2つの粒子の間に働く力関係を計算するためのモデルを定義します。 
# CohFrictMatクラスは粘着力を考慮できる材料のクラスです。
mat_sp = CohFrictMat(alphaKr=0.5,
                     alphaKtw=0.5,
                     young=young,
                     poisson=0.33,
                     frictionAngle=frictangle,
                     density=density)
O.materials.append(mat_sp)

# 粒子をランダムに生成するコードです。粒子径、粒子径のばらつき、生成する範囲、生成した粒子の情報を格納しておく場所、シード値を入力します。
sp = pack.randomPeriPack(radius=0.10, 
                         rRelFuzz=0,
                         initSize=Vector3(6, 6, 2),
                         memoizeDb='/tmp/packDb.sqlite', 
                         seed=-1)

# 上のコードで生成した粒子をOmegaインスタンスに渡します。ここで粒子の色付けを行うことができます。
sp.toSimulation(color=(0, 0, 1))

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
    doneHook = "cyclic_shear()",
    
    # 最大のひずみ速度の値を代入します。
    maxStrainRate = consolidation_max_strain_rate
)

# 繰り返しせん断の反転条件が応力の場合と、ひずみの場合でこの後のインスタンス生成時に代入するパラメータが変わってきます。
# なのでここでif文を使って条件整理をします。
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
    zxPath_temp = ((0.999999, 1),)
else:
    goal_forward_temp = (-consolidation_stress,
                         -consolidation_stress, 
                         -consolidation_stress,
                         0, cyclic_shear_strain_amplitude, 0)
    goal_backward_temp = (-consolidation_stress,
                          -consolidation_stress, 
                          -consolidation_stress,
                          0, cyclic_shear_strain_amplitude, 0)
    
    # ここでstressMaskの説明をもう一回します。今回の場合の対応関係は以下のとおりです。
    # stressMask = 0b    1    0    1    1    1    1
    #                   01   02   12   22   11   00
    # 0の場合はそのテンソルの目標値がひずみであること、1の場合は応力であることを表すので、
    # 今回の場合は02成分(xz成分)のみが目標値がひずみであり、他の5成分は全て応力であることを示します。
    stressMask_temp = 0b101111

# 繰り返し載荷(せん断力を増やす際)の際の境界制御のクラスです。
Peri3D_cyclic_forward = Peri3dController(
    goal = goal_forward_temp,
    stressMask = stressMask_temp,
    nSteps = cyclic_loading_nSteps,
    doneHook = "cyclic_shear()",
    maxStrainRate = cyclic_loading_max_strain_rate,
    xxPath = ((0, 1), (1, 1)),
    yyPath = ((0, 1), (1, 1)),
    zzPath = ((0, 1), (1, 1)),
    zxPath = ((0.999999, 1),)
)

# 繰り返し載荷(せん断力を減らす際際)の際の境界制御のクラスです。
Peri3D_cyclic_backward = Peri3dController(
    goal = goal_backward_temp,
    stressMask = stressMask_temp,
    nSteps = cyclic_loading_nSteps,
    doneHook = "cyclic_shear()",
    maxStrainRate = cyclic_loading_max_strain_rate,
    xxPath = ((0, 1), (1, 1)),
    yyPath = ((0, 1), (1, 1)),
    zzPath = ((0, 1), (1, 1)),
    zxPath = ((0.999999, 1),)
)

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
    Peri3D_iso,
    NewtonIntegrator(damping=0.2),
    PyRunner(realPeriod=0.5, command="addPlotData()")
]

# 解析の際の時間ステップを決めるコードです。
# MEMO: 本来であればすべての2粒子間接触から求まる固有周期の最小値を時間ステップとして採用するべきなんですが、初期状態では接触点の数が0の場合もあるので、近似的に粒子のP波速度の0.1倍を使っています。(Class Referenceと計算式が違うがほとんどのサンプルコードがこれぐらいの値を採用している...)
O.dt = .1 * PWaveTimeStep()

# 繰り返しせん断用の関数です。
# MEMO
def cyclic_shear():
    global flag_cyclic_loading, flag_cyclic_forward
    global consolidation_stress, cyclic_shear_stress_amplitude, current_num_cycle, target_num_cycle, shear_stress_torrelance
    
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
            
            flag_under_torrelance = (cyclic_shear_stress_amplitude - s02) < shear_stress_torrelance           
           
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
                        
            flag_under_torrelance = (-cyclic_shear_stress_amplitude - s02) > -shear_stress_torrelance
            
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
