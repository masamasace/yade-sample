from yade import pack
from yade.gridpfacet import *

dp = 6e-3

periodic = 1
cylinderBottom = 0
cylindersBottom = 1

ground = 0.2*dp
H = 60*dp
L = 20*dp

O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb(
    ), Bo1_GridConnection_Aabb()], allowBiggerThanPeriod=True),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom(),
         Ig2_Sphere_GridConnection_ScGridCoGeom()],
        [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
        [Law2_ScGeom_ViscElPhys_Basic()]
        # [Ip2_FrictMat_FrictMat_FrictPhys()],
        # [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    PyRunner(command='check()', virtPeriod=0.1),
    NewtonIntegrator(damping=0., gravity=(0, 0, -9.81))
]

# Material creation
O.materials.append(ViscElMat(en=0.5, et=1., young=5e6,
                   poisson=0.5, density=2500., frictionAngle=0.4, label='Mat'))
# O.materials.append(FrictMat(young = 5e6, poisson = 0.5, density=2500,frictionAngle=0.4, label='Mat'))

# Periodic Cell or containing box
if periodic:
    O.periodic = True
    O.cell.setBox(L, L, H)
    O.bodies.append(box(center=(L/2., L/2., ground), extents=(100, 100, 0),
                    fixed=True, color=(0., 1., 0.), wire=True, material='Mat'))
else:
    O.bodies.append(box(center=(L/2., L/2., ground), extents=(100, 100, 0), wire=True,
                    fixed=True, color=(0., 1., 0.), material='Mat'))  # Made invisible to see below
    O.bodies.append(box(center=(0, L/2., H/2.), extents=(0, L/2., H/2.), wire=True,
                    fixed=True, color=(0., 1., 0.), material='Mat'))  # Made invisible to see
    O.bodies.append(box(center=(L, L/2., H/2.), extents=(0, L/2.,
                    H/2.), fixed=True, color=(0., 1., 0.), material='Mat'))
    O.bodies.append(box(center=(L/2., 0, H/2.), extents=(L/2., 0, H/2.), wire=True,
                    fixed=True, color=(0., 1., 0.), material='Mat'))  # Made invisible to see
    O.bodies.append(box(center=(L/2., L, H/2.), extents=(L/2.,
                    0, H/2.), fixed=True, color=(0., 1., 0.), material='Mat'))

# Bottom fixed cylinder or fixed particles
if cylinderBottom:
    n = len(O.bodies)
    cylinder(begin=(L/2., 100*L+0.001*dp, ground+dp/2.), end=(L/2., -100*L+0.001*dp, ground+dp/2.), radius=dp /
             2., fixed=True, color=(0, 0, 1), intMaterial='Mat', extMaterial='Mat')  # Made invisible to see inside
else:
    if cylindersBottom:
        ep = 0.001*dp  # a small thing
        nodes = [[L/2., ep, ground+dp/2.], [L/2., 0.33*L+ep, ground+dp/2.],
                 [L/2., 0.33*2*L+ep, ground+dp/2.], [L/2., L+ep, ground+dp/2.]]
        cylinderConnection([[L/2., ep, ground+dp/2.], [L/2., 0.33*L+ep, ground+dp/2.], [
                           L/2., 0.33*2*L+ep, ground+dp/2.], [L/2., L+ep, ground+dp/2.]], dp/2., fixed=True)
    else:
        for n in range(0, int(L/dp)+1):
            O.bodies.append(sphere(center=(L/2., n*dp, ground+dp/2.), radius=dp/2., fixed=True,
                            color=(0, 0, 1), wire=True, material='Mat'))  # Made invisible to see inside

# Particle cloud for gravity deposition
partCloud = pack.SpherePack()
partCloud.makeCloud(minCorner=(4*L/10., 4*L/10., ground+3.*dp),
                    maxCorner=(6*L/10., 6*L/10., 0.9*H), rRelFuzz=0., rMean=dp/2., num=200)
partCloud.toSimulation(material='Mat')

O.dt = 1e-5
O.saveTmp()

# Function to check if the center of a particle is contained inside the cylinder or pseudo cylinder


def check():
    for b in O.bodies:
        if b.dynamic and b.state.pos[2] < ground+dp:
            posCylXZ = Vector3(L/2., 0., ground+dp/2.)
            pRelX = b.state.pos[0] % L - posCylXZ[0]
            pRelZ = b.state.pos[2] - posCylXZ[2]
            # If the particle center is closer than the radius (i.e. if it is completely inside!)
            if sqrt(pow(pRelX, 2) + pow(pRelZ, 2)) < dp/2.:
                delta = dp/2.-sqrt(pow(pRelX, 2) + pow(pRelZ, 2))
                print('particle {0} is inside the cylinder from an amount of {1} diameter'.format(
                    b.id, delta/dp))
                O.pause()
