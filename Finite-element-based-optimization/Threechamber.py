# Environment for three chamber FE model (gym structure with pyansys FE implementation)

from array import array
import csv
import gym
from gym import spaces
from gym.utils import seeding
import math
import matplotlib.pyplot as plt
import numpy as np
from os import path
import pyansys
import pyvista

class ThreechamberEnv(gym.Env):

    def __init__(self):
        s = 0.5
        high = np.array([np.inf] * 2)
        self.action_space = spaces.Box(np.array([-s, -s, -s, -s,-s,-s,-s,-s]),
                                       np.array([s, s, s, s,s,s,s,s]),
                                       dtype=np.float32)
        self.observation_space = spaces.Box(-high, high, dtype=np.float32)
        self.seed()

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]

    def step(self, action):
        displacement, Prev = self.state
        Last = 1

        length = 0.0508 # Dome dimensions are 3" x 2" x 2"
        width = 0.06896

        mapdl = pyansys.launch_mapdl(override=True, additional_switches='-smp',
                                     loglevel='ERROR')

        mapdl.clear('NOSTART')
        mapdl.prep7()

        # Define bounds for chamber vertices (All units SI)
        min_x = 0.01
        min_y = 0.005
        med_y = 0.03
        max_x = 0.025-0.006
        max_y = 0.06339
        depth = 0.05

        # Define ordered pairs for chamber vertices 8-DoF Configuration
        a1x = np.clip(action[0], min_x, max_x - 0.005)
        a1y = np.clip(action[1], min_y, max_y - 0.005)
        #a1y = min_y
        a2x = np.clip(action[2], a1x+0.004, max_x)
        a2y = np.clip(action[3], min_y, max_y - 0.005)
        #a2y = min_y
        a3x = np.clip(action[4], min_x + 0.005, max_x)
        #a3y = max_y
        a3y = np.clip(action[5], a2y + 0.001, max_y)
        a4x = np.clip(action[6], min_x, a3x - 0.004)
        #a4y = max_y
        a4y = np.clip(action[7], a1y + 0.001, max_y)

        listv = [a1x,a1y,a2x,a2y,a3x,a3y,a4x,a4y]
        # Set condition so that slope of L1 & L4, L2 & L3 are not undef. for frontal area calc.
        if a1x == a4x:
            a1x = a1x-0.0001
        if a2x == a3x:
            a3x = a3x-0.0001
        # Conditional Statement for likely geometry overlap
        if a3x <= a1x:
            a3x = a2x+0.001
        if a3x >= a2x:
            a4x = a2x - 0.001 # extra mm added for some mesh problems when segments are too close
        if a4x > a2x:
            a4x = a2x - 0.001
        if a3y < a2y and a3x >= a2x:
            a3y = a2y+ 0.001 # prevents crossing line segments helps connection counter clockwise, often breaks if CW
        if a3y < a1y:
            a3y = a1y + 0.001
        if a4y < a2y:
            a4y = a2y +0.001

        # Define Left chamber vertices by keypoints and connect
        mapdl.k(1, a1x, a1y, 0)
        mapdl.k(2, a2x, a2y, 0)
        mapdl.k(3, a3x, a3y, 0)
        mapdl.k(4, a4x, a4y, 0)
        mapdl.l(1, 2)
        mapdl.l(2, 3)
        mapdl.l(3, 4)
        mapdl.l(4, 1)
        A1 = mapdl.al(1, 2, 3, 4)

        # Mirror Left chamber vertices to form right side
        a5x = length - a1x
        a6x = length - a2x
        a7x = length - a3x
        a8x = length - a4x

        # Define Right chamber vertices by keypoints and connect
        mapdl.k(5, a5x, a1y, 0)
        mapdl.k(6, a6x, a2y, 0)
        mapdl.k(7, a7x, a3y, 0)
        mapdl.k(8, a8x, a4y, 0)
        mapdl.l(5, 6)
        mapdl.l(6, 7)
        mapdl.l(7, 8)
        mapdl.l(8, 5)
        A2 = mapdl.al(5, 6, 7, 8)

        size = 0.00075

        # Define conditions (constant)
        Pressure1 = 68947.6  # Seperate pressure for sides defined for easy study of asymmetric chamber pressures
        Pressure2 = -68947.6

        # Calculate chamber side lengths (symmetric chambers only one side needs to be determined)
        d1 = math.sqrt((a2x - a1x) ** 2 + (a2y - a1y) ** 2)
        d2 = math.sqrt((a3x - a2x) ** 2 + (a3y - a2y) ** 2)
        d3 = math.sqrt((a4x - a3x) ** 2 + (a4y - a3y) ** 2)
        d4 = math.sqrt((a4x - a1x) ** 2 + (a4y - a1y) ** 2)

        # Determine the slopes and angles wrt. horz. for each segment
        Slope1 = (a2y - a1y) / (a2x - a1x)
        Slope2 = (a3y - a2y) / (a3x - a2x)
        Slope3 = (a4y - a3y) / (a4x - a3x)
        Slope4 = (a4y - a1y) / (a4x - a1x)
        th1 = math.atan(Slope1)
        th2 = math.atan(Slope2)
        th3 = math.atan(Slope3)
        th4 = math.atan(Slope4)
        m_area1 = math.ceil(d1 / size) ** -1
        m_area2 = math.ceil(d2 / size) ** -1
        m_area3 = math.ceil(d3 / size) ** -1
        m_area4 = math.ceil(d4 / size) ** -1
        area1 = (depth * d1)
        area2 = (depth * d2)
        area3 = (depth * d3)
        area4 = (depth * d4)

        # Estimate the load for each reigon by components
        Load_est1x = (m_area1 * area1 * Pressure1 * math.cos(th1)) #/ 2
        Load_est1y = (m_area1 * area1 * Pressure1 * math.sin(th1)) #/ 2
        Load_est2x = (m_area2 * area2 * Pressure1 * math.cos(th2)) #/ 2
        Load_est2y = (m_area2 * area2 * Pressure1 * math.sin(th2)) #/ 2
        Load_est3x = (m_area3 * area3 * Pressure1 * math.cos(th3)) #/ 2
        Load_est3y = (m_area3 * area3 * Pressure1 * math.sin(th3)) #/ 2
        Load_est4x = (m_area4 * area4 * Pressure1 * math.cos(th4)) #/ 2
        Load_est4y = (m_area4 * area4 * Pressure1 * math.sin(th4)) #/ 2

        Load_est5x = (m_area1 * area1 * Pressure2 * math.cos(th1)) #/ 2
        Load_est5y = (m_area1 * area1 * Pressure2 * math.sin(th1)) #/ 2
        Load_est6x = (m_area2 * area2 * Pressure2 * math.cos(th2)) #/ 2
        Load_est6y = (m_area2 * area2 * Pressure2 * math.sin(th2)) #/ 2
        Load_est7x = (m_area2 * area3 * Pressure2 * math.cos(th3)) #/ 2
        Load_est7y = (m_area3 * area3 * Pressure2 * math.sin(th3)) #/ 2
        Load_est8x = (m_area4 * area4 * Pressure2 * math.cos(th4)) #/ 2
        Load_est8y = (m_area4 * area4 * Pressure2 * math.sin(th4)) #/ 2

        # Additional Points manual calc. Left chamber
        # Section for segment 1-4
        mid14y = (a4y + a1y) / 2
        mid14x = (a4x + a1x) / 2

        mid14ty = (a4y + mid14y) / 2
        mid14tx = (a4x + mid14x) / 2

        mid14tty = (a4y + mid14ty) / 2
        mid14ttx = (a4x + mid14tx) / 2

        mid14by = (a1y + mid14y) / 2
        mid14bx = (a1x + mid14x) / 2

        mid14bby = (a1y + mid14by) / 2
        mid14bbx = (a1x + mid14bx) / 2

        # Section for segment 1-2
        mid12y = (a2y + a1y) / 2
        mid12x = (a2x + a1x) / 2

        mid12ty = (a2y + mid12y) / 2
        mid12tx = (a2x + mid12x) / 2

        mid12tty = (a2y + mid12ty) / 2
        mid12ttx = (a2x + mid12tx) / 2

        mid12by = (a1y + mid12y) / 2
        mid12bx = (a1x + mid12x) / 2

        mid12bby = (a1y + mid12by) / 2
        mid12bbx = (a1x + mid12bx) / 2

        # Section for segment 2-3
        mid23y = (a3y + a2y) / 2
        mid23x = (a3x + a2x) / 2

        mid23ty = (a3y + mid23y) / 2
        mid23tx = (a3x + mid23x) / 2

        mid23tty = (a3y + mid23ty) / 2
        mid23ttx = (a3x + mid23tx) / 2

        mid23by = (a2y + mid23y) / 2
        mid23bx = (a2x + mid23x) / 2

        mid23bby = (a2y + mid23by) / 2
        mid23bbx = (a2x + mid23bx) / 2

        # section for segment 3-4
        mid34y = (a3y + a4y) / 2
        mid34x = (a3x + a4x) / 2

        mid34ty = (a3y + mid34y) / 2
        mid34tx = (a3x + mid34x) / 2

        mid34tty = (a3y + mid34ty) / 2
        mid34ttx = (a3x + mid34tx) / 2

        mid34by = (a4y + mid34y) / 2
        mid34bx = (a4x + mid34x) / 2

        mid34bby = (a4y + mid34by) / 2
        mid34bbx = (a4x + mid34bx) / 2

        # Additional Points manual calc. Right chamber
        # Section for segment 5-8
        mid58y = mid14y
        mid58x = length - mid14x
        mid58ty = mid14ty
        mid58tx = length - mid14tx
        mid58tty = mid14tty
        mid58ttx = length - mid14ttx
        mid58by = mid14by
        mid58bx = length - mid14bx
        mid58bby = mid14bby
        mid58bbx = length - mid14bbx

        # Section for segment 5-6
        mid56y = mid12y
        mid56x = length - mid12x
        mid56ty = mid12ty
        mid56tx = length - mid12tx
        mid56tty = mid12tty
        mid56ttx = length - mid12ttx
        mid56by = mid12by
        mid56bx = length - mid12bx
        mid56bby = mid12bby
        mid56bbx = length - mid12bbx

        # Section for segment 6-7
        mid67y = mid23y
        mid67x = length - mid23x
        mid67ty = mid23ty
        mid67tx = length - mid23tx
        mid67tty = mid23tty
        mid67ttx = length - mid23ttx
        mid67by = mid23by
        mid67bx = length - mid23bx
        mid67bby = mid23bby
        mid67bbx = length - mid23bbx

        # section for segment 7-8
        mid78y = mid34y
        mid78x = length - mid34x
        mid78ty = mid34ty
        mid78tx = length - mid34tx
        mid78tty = mid34tty
        mid78ttx = length - mid34ttx
        mid78by = mid34by
        mid78bx = length - mid34bx
        mid78bby = mid34bby
        mid78bbx = length - mid34bbx

        # FE/ ANSYS implementation
        # Define elements and material prop.
        mapdl.units('SI')  # SI - International system (m, kg, s, K).

        # Definition of element tyoe and thickness
        mapdl.et(1, "PLANE183", kop3=3)
        mapdl.r(1, depth)

        # Definition of Material Properties
        mapdl.mp('EX', 1, 404E3)  # Youngs modulus (Pa)
        mapdl.mp('DENS', 1, 2330)  # Density (kg/m^3)
        mapdl.mp('NUXY', 1, 0.3988)  # Poisson's Ratio

        # create the Dome area
        rect_anum = mapdl.blc4(width=length, height=width)

        # Remove Chamber areas from dome area
        plate_with_hole_anum = mapdl.asba(rect_anum, "ALL")
        # Specify Element Size and Mesh
        esize = 0.00075
        mapdl.esize(esize)
        mapdl.amesh(plate_with_hole_anum)

        # Constrain bottom of the dome
        mapdl.nsel('S', 'LOC', 'Y', 0)
        mapdl.d('ALL', 'UY')
        mapdl.nsel('R', 'LOC', 'X', 0, length)
        mapdl.d('ALL', 'UX')

        # Loading of chambers using calculated points Horizontal (all)
        # Chamber Vertices Left side
        mapdl.nsel('S', 'LOC', 'Y', a1y)
        mapdl.nsel('R', 'LOC', 'X', a1x)
        mapdl.f('ALL', 'FY', Load_est1x)

        mapdl.nsel('S', 'LOC', 'Y', a2y)
        mapdl.nsel('R', 'LOC', 'X', a2x)
        mapdl.f('ALL', 'FY', Load_est2x)

        mapdl.nsel('S', 'LOC', 'Y', a3y)
        mapdl.nsel('R', 'LOC', 'X', a3x)
        mapdl.f('ALL', 'FY', Load_est3x)

        mapdl.nsel('S', 'LOC', 'Y', a4y)
        mapdl.nsel('R', 'LOC', 'X', a4x)
        mapdl.f('ALL', 'FY', Load_est4x)

        # Chamber Vertices Right side
        mapdl.nsel('S', 'LOC', 'Y', a1y)
        mapdl.nsel('R', 'LOC', 'X', a5x)
        mapdl.f('ALL', 'FY', Load_est5x)

        mapdl.nsel('S', 'LOC', 'Y', a2y)
        mapdl.nsel('R', 'LOC', 'X', a6x)
        mapdl.f('ALL', 'FY', Load_est6x)

        mapdl.nsel('S', 'LOC', 'Y', a3y)
        mapdl.nsel('R', 'LOC', 'X', a7x)
        mapdl.f('ALL', 'FY', Load_est7x)

        mapdl.nsel('S', 'LOC', 'Y', a4y)
        mapdl.nsel('R', 'LOC', 'X', a8x)
        mapdl.f('ALL', 'FY', Load_est8x)

        # Loading of points along chamber edges
        # Line 1-4 (Left Chamber)
        mapdl.nsel('S', 'LOC', 'Y', mid14y)
        mapdl.nsel('R', 'LOC', 'X', mid14x)
        mapdl.f('ALL', 'FY', Load_est4x)

        mapdl.nsel('S', 'LOC', 'Y', mid14ty)
        mapdl.nsel('R', 'LOC', 'X', mid14tx)
        mapdl.f('ALL', 'FY', Load_est4x)

        mapdl.nsel('S', 'LOC', 'Y', mid14tty)
        mapdl.nsel('R', 'LOC', 'X', mid14ttx)
        mapdl.f('ALL', 'FY', Load_est4x)

        mapdl.nsel('S', 'LOC', 'Y', mid14by)
        mapdl.nsel('R', 'LOC', 'X', mid14bx)
        mapdl.f('ALL', 'FY', Load_est4x)

        mapdl.nsel('S', 'LOC', 'Y', mid14bby)
        mapdl.nsel('R', 'LOC', 'X', mid14bbx)
        mapdl.f('ALL', 'FY', Load_est4x)

        # Line 1-2
        mapdl.nsel('S', 'LOC', 'Y', mid12y)
        mapdl.nsel('R', 'LOC', 'X', mid12x)
        mapdl.f('ALL', 'FY', Load_est1x)

        mapdl.nsel('S', 'LOC', 'Y', mid12ty)
        mapdl.nsel('R', 'LOC', 'X', mid12tx)
        mapdl.f('ALL', 'FY', Load_est1x)

        mapdl.nsel('S', 'LOC', 'Y', mid12tty)
        mapdl.nsel('R', 'LOC', 'X', mid12ttx)
        mapdl.f('ALL', 'FY', Load_est1x)

        mapdl.nsel('S', 'LOC', 'Y', mid12by)
        mapdl.nsel('R', 'LOC', 'X', mid12bx)
        mapdl.f('ALL', 'FY', Load_est1x)

        mapdl.nsel('S', 'LOC', 'Y', mid12bby)
        mapdl.nsel('R', 'LOC', 'X', mid12bbx)
        mapdl.f('ALL', 'FY', Load_est1x)

        # Line 2-3
        mapdl.nsel('S', 'LOC', 'Y', mid23y)
        mapdl.nsel('R', 'LOC', 'X', mid23x)
        mapdl.f('ALL', 'FY', Load_est2x)

        mapdl.nsel('S', 'LOC', 'Y', mid23ty)
        mapdl.nsel('R', 'LOC', 'X', mid23tx)
        mapdl.f('ALL', 'FY', Load_est2x)

        mapdl.nsel('S', 'LOC', 'Y', mid23tty)
        mapdl.nsel('R', 'LOC', 'X', mid23ttx)
        mapdl.f('ALL', 'FY', Load_est2x)

        mapdl.nsel('S', 'LOC', 'Y', mid23by)
        mapdl.nsel('R', 'LOC', 'X', mid23bx)
        mapdl.f('ALL', 'FY', Load_est2x)

        mapdl.nsel('S', 'LOC', 'Y', mid23bby)
        mapdl.nsel('R', 'LOC', 'X', mid23bbx)
        mapdl.f('ALL', 'FY', Load_est2x)

        # Line 3-4
        mapdl.nsel('S', 'LOC', 'Y', mid34y)
        mapdl.nsel('R', 'LOC', 'X', mid34x)
        mapdl.f('ALL', 'FY', Load_est3x)

        mapdl.nsel('S', 'LOC', 'Y', mid34ty)
        mapdl.nsel('R', 'LOC', 'X', mid34tx)
        mapdl.f('ALL', 'FY', Load_est3x)

        mapdl.nsel('S', 'LOC', 'Y', mid34tty)
        mapdl.nsel('R', 'LOC', 'X', mid34ttx)
        mapdl.f('ALL', 'FY', Load_est3x)

        mapdl.nsel('S', 'LOC', 'Y', mid34by)
        mapdl.nsel('R', 'LOC', 'X', mid34bx)
        mapdl.f('ALL', 'FY', Load_est3x)

        mapdl.nsel('S', 'LOC', 'Y', mid34bby)
        mapdl.nsel('R', 'LOC', 'X', mid34bbx)
        mapdl.f('ALL', 'FY', Load_est3x)

        # Line 5-8 (Right Chamber)
        mapdl.nsel('S', 'LOC', 'Y', mid58y)
        mapdl.nsel('R', 'LOC', 'X', mid58x)
        mapdl.f('ALL', 'FY', Load_est8x)

        mapdl.nsel('S', 'LOC', 'Y', mid58ty)
        mapdl.nsel('R', 'LOC', 'X', mid58tx)
        mapdl.f('ALL', 'FY', Load_est8x)

        mapdl.nsel('S', 'LOC', 'Y', mid58tty)
        mapdl.nsel('R', 'LOC', 'X', mid58ttx)
        mapdl.f('ALL', 'FY', Load_est8x)

        mapdl.nsel('S', 'LOC', 'Y', mid58by)
        mapdl.nsel('R', 'LOC', 'X', mid58bx)
        mapdl.f('ALL', 'FY', Load_est8x)

        mapdl.nsel('S', 'LOC', 'Y', mid58bby)
        mapdl.nsel('R', 'LOC', 'X', mid58bbx)
        mapdl.f('ALL', 'FY', Load_est8x)

        # Line 5-6
        mapdl.nsel('S', 'LOC', 'Y', mid56y)
        mapdl.nsel('R', 'LOC', 'X', mid56x)
        mapdl.f('ALL', 'FY', Load_est5x)

        mapdl.nsel('S', 'LOC', 'Y', mid56ty)
        mapdl.nsel('R', 'LOC', 'X', mid56tx)
        mapdl.f('ALL', 'FY', Load_est5x)

        mapdl.nsel('S', 'LOC', 'Y', mid56tty)
        mapdl.nsel('R', 'LOC', 'X', mid56ttx)
        mapdl.f('ALL', 'FY', Load_est5x)

        mapdl.nsel('S', 'LOC', 'Y', mid56by)
        mapdl.nsel('R', 'LOC', 'X', mid56bx)
        mapdl.f('ALL', 'FY', Load_est5x)

        mapdl.nsel('S', 'LOC', 'Y', mid56bby)
        mapdl.nsel('R', 'LOC', 'X', mid56bbx)
        mapdl.f('ALL', 'FY', Load_est5x)

        # Line 6-7
        mapdl.nsel('S', 'LOC', 'Y', mid67y)
        mapdl.nsel('R', 'LOC', 'X', mid67x)
        mapdl.f('ALL', 'FY', Load_est6x)

        mapdl.nsel('S', 'LOC', 'Y', mid67ty)
        mapdl.nsel('R', 'LOC', 'X', mid67tx)
        mapdl.f('ALL', 'FY', Load_est6x)

        mapdl.nsel('S', 'LOC', 'Y', mid67tty)
        mapdl.nsel('R', 'LOC', 'X', mid67ttx)
        mapdl.f('ALL', 'FY', Load_est6x)

        mapdl.nsel('S', 'LOC', 'Y', mid67by)
        mapdl.nsel('R', 'LOC', 'X', mid67bx)
        mapdl.f('ALL', 'FY', Load_est6x)

        mapdl.nsel('S', 'LOC', 'Y', mid67bby)
        mapdl.nsel('R', 'LOC', 'X', mid67bbx)
        mapdl.f('ALL', 'FY', Load_est6x)

        # Line 7-8
        mapdl.nsel('S', 'LOC', 'Y', mid78y)
        mapdl.nsel('R', 'LOC', 'X', mid78x)
        mapdl.f('ALL', 'FY', Load_est7x)

        mapdl.nsel('S', 'LOC', 'Y', mid78ty)
        mapdl.nsel('R', 'LOC', 'X', mid78tx)
        mapdl.f('ALL', 'FY', Load_est7x)

        mapdl.nsel('S', 'LOC', 'Y', mid78tty)
        mapdl.nsel('R', 'LOC', 'X', mid78ttx)
        mapdl.f('ALL', 'FY', Load_est7x)

        mapdl.nsel('S', 'LOC', 'Y', mid78by)
        mapdl.nsel('R', 'LOC', 'X', mid78bx)
        mapdl.f('ALL', 'FY', Load_est7x)

        mapdl.nsel('S', 'LOC', 'Y', mid78bby)
        mapdl.nsel('R', 'LOC', 'X', mid78bbx)
        mapdl.f('ALL', 'FY', Load_est7x)

        # Verticle Loading of nodes (Manual)
        # Vertices
        mapdl.nsel('S', 'LOC', 'Y', a1y)
        mapdl.nsel('R', 'LOC', 'X', a1x)
        mapdl.f('ALL', 'FX', Load_est1y)

        mapdl.nsel('S', 'LOC', 'Y', a2y)
        mapdl.nsel('R', 'LOC', 'X', a2x)
        mapdl.f('ALL', 'FX', Load_est2y)

        mapdl.nsel('S', 'LOC', 'Y', a3y)
        mapdl.nsel('R', 'LOC', 'X', a3x)
        mapdl.f('ALL', 'FX', Load_est3y)

        mapdl.nsel('S', 'LOC', 'Y', a4y)
        mapdl.nsel('R', 'LOC', 'X', a4x)
        mapdl.f('ALL', 'FX', Load_est4y)

        # Chamber Vertices Right side
        mapdl.nsel('S', 'LOC', 'Y', a1y)
        mapdl.nsel('R', 'LOC', 'X', a5x)
        mapdl.f('ALL', 'FX', Load_est5y)

        mapdl.nsel('S', 'LOC', 'Y', a2y)
        mapdl.nsel('R', 'LOC', 'X', a6x)
        mapdl.f('ALL', 'FX', Load_est6y)

        mapdl.nsel('S', 'LOC', 'Y', a3y)
        mapdl.nsel('R', 'LOC', 'X', a7x)
        mapdl.f('ALL', 'FX', Load_est7y)

        mapdl.nsel('S', 'LOC', 'Y', a4y)
        mapdl.nsel('R', 'LOC', 'X', a8x)
        mapdl.f('ALL', 'FX', Load_est8y)

        # Loading of points along chamber edges
        # Line 1-4 (Left Chamber)
        mapdl.nsel('S', 'LOC', 'Y', mid14y)
        mapdl.nsel('R', 'LOC', 'X', mid14x)
        mapdl.f('ALL', 'FX', Load_est4y)

        mapdl.nsel('S', 'LOC', 'Y', mid14ty)
        mapdl.nsel('R', 'LOC', 'X', mid14tx)
        mapdl.f('ALL', 'FX', Load_est4y)

        mapdl.nsel('S', 'LOC', 'Y', mid14tty)
        mapdl.nsel('R', 'LOC', 'X', mid14ttx)
        mapdl.f('ALL', 'FX', Load_est4y)

        mapdl.nsel('S', 'LOC', 'Y', mid14by)
        mapdl.nsel('R', 'LOC', 'X', mid14bx)
        mapdl.f('ALL', 'FX', Load_est4y)

        mapdl.nsel('S', 'LOC', 'Y', mid14bby)
        mapdl.nsel('R', 'LOC', 'X', mid14bbx)
        mapdl.f('ALL', 'FX', Load_est4y)

        # Line 1-2
        mapdl.nsel('S', 'LOC', 'Y', mid12y)
        mapdl.nsel('R', 'LOC', 'X', mid12x)
        mapdl.f('ALL', 'FX', Load_est1y)

        mapdl.nsel('S', 'LOC', 'Y', mid12ty)
        mapdl.nsel('R', 'LOC', 'X', mid12tx)
        mapdl.f('ALL', 'FX', Load_est1y)

        mapdl.nsel('S', 'LOC', 'Y', mid12tty)
        mapdl.nsel('R', 'LOC', 'X', mid12ttx)
        mapdl.f('ALL', 'FX', Load_est1y)

        mapdl.nsel('S', 'LOC', 'Y', mid12by)
        mapdl.nsel('R', 'LOC', 'X', mid12bx)
        mapdl.f('ALL', 'FX', Load_est1y)

        mapdl.nsel('S', 'LOC', 'Y', mid12bby)
        mapdl.nsel('R', 'LOC', 'X', mid12bbx)
        mapdl.f('ALL', 'FX', Load_est1y)

        # Line 2-3
        mapdl.nsel('S', 'LOC', 'Y', mid23y)
        mapdl.nsel('R', 'LOC', 'X', mid23x)
        mapdl.f('ALL', 'FX', Load_est2y)

        mapdl.nsel('S', 'LOC', 'Y', mid23ty)
        mapdl.nsel('R', 'LOC', 'X', mid23tx)
        mapdl.f('ALL', 'FX', Load_est2y)

        mapdl.nsel('S', 'LOC', 'Y', mid23tty)
        mapdl.nsel('R', 'LOC', 'X', mid23ttx)
        mapdl.f('ALL', 'FX', Load_est2y)

        mapdl.nsel('S', 'LOC', 'Y', mid23by)
        mapdl.nsel('R', 'LOC', 'X', mid23bx)
        mapdl.f('ALL', 'FX', Load_est2y)

        mapdl.nsel('S', 'LOC', 'Y', mid23bby)
        mapdl.nsel('R', 'LOC', 'X', mid23bbx)
        mapdl.f('ALL', 'FX', Load_est2y)

        # Line 3-4
        mapdl.nsel('S', 'LOC', 'Y', mid34y)
        mapdl.nsel('R', 'LOC', 'X', mid34x)
        mapdl.f('ALL', 'FX', Load_est3y)

        mapdl.nsel('S', 'LOC', 'Y', mid34ty)
        mapdl.nsel('R', 'LOC', 'X', mid34tx)
        mapdl.f('ALL', 'FX', Load_est3y)

        mapdl.nsel('S', 'LOC', 'Y', mid34tty)
        mapdl.nsel('R', 'LOC', 'X', mid34ttx)
        mapdl.f('ALL', 'FX', Load_est3y)

        mapdl.nsel('S', 'LOC', 'Y', mid34by)
        mapdl.nsel('R', 'LOC', 'X', mid34bx)
        mapdl.f('ALL', 'FX', Load_est3y)

        mapdl.nsel('S', 'LOC', 'Y', mid34bby)
        mapdl.nsel('R', 'LOC', 'X', mid34bbx)
        mapdl.f('ALL', 'FX', Load_est3y)

        # Line 5-8 (Right Chamber)
        mapdl.nsel('S', 'LOC', 'Y', mid58y)
        mapdl.nsel('R', 'LOC', 'X', mid58x)
        mapdl.f('ALL', 'FX', Load_est8y)

        mapdl.nsel('S', 'LOC', 'Y', mid58ty)
        mapdl.nsel('R', 'LOC', 'X', mid58tx)
        mapdl.f('ALL', 'FX', Load_est8y)

        mapdl.nsel('S', 'LOC', 'Y', mid58tty)
        mapdl.nsel('R', 'LOC', 'X', mid58ttx)
        mapdl.f('ALL', 'FX', Load_est8y)

        mapdl.nsel('S', 'LOC', 'Y', mid58by)
        mapdl.nsel('R', 'LOC', 'X', mid58bx)
        mapdl.f('ALL', 'FX', Load_est8y)

        mapdl.nsel('S', 'LOC', 'Y', mid58bby)
        mapdl.nsel('R', 'LOC', 'X', mid58bbx)
        mapdl.f('ALL', 'FX', Load_est8y)

        # Line 5-6
        mapdl.nsel('S', 'LOC', 'Y', mid56y)
        mapdl.nsel('R', 'LOC', 'X', mid56x)
        mapdl.f('ALL', 'FX', Load_est5y)

        mapdl.nsel('S', 'LOC', 'Y', mid56ty)
        mapdl.nsel('R', 'LOC', 'X', mid56tx)
        mapdl.f('ALL', 'FX', Load_est5y)

        mapdl.nsel('S', 'LOC', 'Y', mid56tty)
        mapdl.nsel('R', 'LOC', 'X', mid56ttx)
        mapdl.f('ALL', 'FX', Load_est5y)

        mapdl.nsel('S', 'LOC', 'Y', mid56by)
        mapdl.nsel('R', 'LOC', 'X', mid56bx)
        mapdl.f('ALL', 'FX', Load_est5y)

        mapdl.nsel('S', 'LOC', 'Y', mid56bby)
        mapdl.nsel('R', 'LOC', 'X', mid56bbx)
        mapdl.f('ALL', 'FX', Load_est5y)

        # Line 6-7
        mapdl.nsel('S', 'LOC', 'Y', mid67y)
        mapdl.nsel('R', 'LOC', 'X', mid67x)
        mapdl.f('ALL', 'FX', Load_est6y)

        mapdl.nsel('S', 'LOC', 'Y', mid67ty)
        mapdl.nsel('R', 'LOC', 'X', mid67tx)
        mapdl.f('ALL', 'FX', Load_est6y)

        mapdl.nsel('S', 'LOC', 'Y', mid67tty)
        mapdl.nsel('R', 'LOC', 'X', mid67ttx)
        mapdl.f('ALL', 'FX', Load_est6y)

        mapdl.nsel('S', 'LOC', 'Y', mid67by)
        mapdl.nsel('R', 'LOC', 'X', mid67bx)
        mapdl.f('ALL', 'FX', Load_est6y)

        mapdl.nsel('S', 'LOC', 'Y', mid67bby)
        mapdl.nsel('R', 'LOC', 'X', mid67bbx)
        mapdl.f('ALL', 'FX', Load_est6y)

        # Line 7-8
        mapdl.nsel('S', 'LOC', 'Y', mid78y)
        mapdl.nsel('R', 'LOC', 'X', mid78x)
        mapdl.f('ALL', 'FX', Load_est7y)

        mapdl.nsel('S', 'LOC', 'Y', mid78ty)
        mapdl.nsel('R', 'LOC', 'X', mid78tx)
        mapdl.f('ALL', 'FX', Load_est7y)

        mapdl.nsel('S', 'LOC', 'Y', mid78tty)
        mapdl.nsel('R', 'LOC', 'X', mid78ttx)
        mapdl.f('ALL', 'FX', Load_est7y)

        mapdl.nsel('S', 'LOC', 'Y', mid78by)
        mapdl.nsel('R', 'LOC', 'X', mid78bx)
        mapdl.f('ALL', 'FX', Load_est7y)

        mapdl.nsel('S', 'LOC', 'Y', mid78bby)
        mapdl.nsel('R', 'LOC', 'X', mid78bbx)
        mapdl.f('ALL', 'FX', Load_est7y)


        # Solve
        _ = mapdl.allsel()
        mapdl.run('/SOLU')
        mapdl.antype('STATIC')
        mapdl.solve()
        result = mapdl.result

        elem, disp = result.nodal_solution(0)
        dispvalx = disp[:, 0]
        dispvaly = disp[:, -1]

        # Additional Information collected
        max_x_displacement = max(dispvalx)
        max_y_displacement = max(dispvaly)

        len_x = len(dispvalx)
        len_y = len(dispvalx)

        Avgdispx=sum((dispvalx))/len_x
        Avgdispy=sum((dispvaly))/len_y

        Reward = 1000*(Avgdispx -10*Avgdispy)#-10*Avgdispy)

        # Allocate space for writing data to a .csv file
        Full_step_list = [a1x, a1y, a2x, a2y, a3x, a3y, a4x, a4y,max_x_displacement,max_y_displacement,Avgdispx,
                          Avgdispy, Reward]

        Full_step_list = [ '%.4f' % elem for elem in Full_step_list ]
        print(Full_step_list)
        mapdl.exit()

        with open("Output data.csv", "a") as out_file:
            out_String = str(Full_step_list)
            out_file.write(out_String + '\n')

        newdisplacement = Avgdispx
        if newdisplacement < displacement:
            Last=0

        self.state = np.array([newdisplacement, Last])
        return self._get_obs(), Reward, True, {}

    def reset(self):
        high = np.array([1,1])
        self.state = self.np_random.uniform(low=-high, high=high)
        self.last_u = None
        return self._get_obs()

    def _get_obs(self):
        displacement, Prev = self.state
        return np.array([displacement, Prev])

    def _render_(self):
        print("This works")
        return



