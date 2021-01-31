import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
from os import path
import math
from array import array

class ThreechamberEnva(gym.Env):
    metadata = {
        'render.modes': ['human', 'rgb_array'],
        'video.frames_per_second': 30
    }

    def __init__(self):

        s = 0.1
        high = np.array([np.inf] * 2)
        self.action_space = spaces.Box(np.array([-s, -s, -s, -s, -s]), np.array([  s, s, s, s, s]),
                                       dtype=np.float32)
        self.observation_space = spaces.Box(-high, high, dtype=np.float32)
        self.seed()

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]

    def step(self, action):
        displacement1, displacement2 = self.state

        max_x = 0.025  # units are in m
        L1 =   0.06339
        L2 = L1
        L1_p = np.clip(action[0], 0.015, max_x)
        b1 = np.clip(action[1],0.005, 0.01)
        L2_3p = np.clip(action[2], 0.006,0.01 )
        b2 = np.clip(action[3], 0.004, max_x)
        h = np.clip(action[4], 0.04, 0.044)

        x=0
        x1 = 0
        x2 = L1/2
        E = 885523.3145

        c1 = b1/L1_p
        c2 = b2/L2_3p

        if c1 >0.98:
            L1_p = b1 +0.001

        if c1 >0.05:
            L1_p =  b1 +0.001

        b0 = b1
        h2 = h
        h1 = h

        L2_2p = ((L2 * L2) + ((b2 - L2_3p) * (b2 - L2_3p))) ** (0.5)
        L1_2p = ((L1 *L1) + ((L1_p - b1) * (L1_p - b1 )))** (0.5)
        d1 = (0.0508 - ((2 * b0) + b2))
        d11 = (0.0508 - ((2 * L1_p) + L2_3p))
        cos1 = (L1 / L1_2p)
        A11 = h1 * L1_2p
        A22 = h2*L2_2p

        P1 = -(103421)  # Units are in Pa
        P3 = (103421)
        P2 = (-P1 + P3)

        A1 = (h * A11)
        q1 = ((P1 * (A1 * cos1)) / L1)
        alpha = ((L1_p - b0) / L1)
        E = 885523.3145
        R = (q1 / ((alpha ** 4) * E * h))
        F1 = q1 * L1

        C2 = (3 * R * b0)
        C4 = ((3 / 4) * R * (b0 * b0))
        C3 = (C4 / (L1_p * L1_p)) - (C2 / L1_p) - ((1.5) * R * ((math.log(L1_p)) + 1))
        C1 = -((C2*(math.log(L1_p))) + (C3*L1_p) + (C4 / L1_p) + (1.5 * R * (math.log(L1_p)) * L1_p))

        Y1 = (C1+(C2 * (math.log((alpha * x) + b0))) + (C3 * ((alpha * x) + b0)) + (C4 / ((alpha * x) + b0)) + (
                    1.5 * R * (math.log((alpha * x) + b0)) * ((alpha * x) + b0)))
        Y11 = (C1 + (C2 * (math.log((alpha * x2) + b0))) + (C3 * ((alpha * x2) + b0)) + (C4 / ((alpha * x2) + b0)) + (
                    1.5 * R * (math.log((alpha * x2) + b0)) * ((alpha * x2) + b0)))

        cos2 = (L2 / L2_2p)
        A2 = (h2 * A22)
        q2 = ((P2 * (A2 * cos2)) / L2)
        alpha2 = ((L2_3p - b2) / L2)
        R2 = (q2 / ((alpha2 ** 4) * E * h2))
        F2 = q2 * L2

        C22 = ((6 / 4) * R2 * b2)
        C42 = ((3 / 8) * R2 * (b2 * b2))
        C32 = (C42 / (L2_3p  * L2_3p)) - (C22 / L2_3p) - ((3 / 4) * R2 * ((math.log(L2_3p)) + 1))
        C12 = -((C22 * (math.log(L2_3p))) + (C32 * L2_3p) + (C42 / L2_3p) + ((3 / 4) * R2 * (math.log(L2_3p)) * L2_3p))

        Y2 = (C12 + (C22 * (math.log((alpha2*x) + b2))) + (C32 * ((alpha2 * x) + b2)) + (C42 / ((alpha2 * x) + b2)) + (
                    (3 / 4) * R2 * (math.log((alpha2 * x) + b2)) * ((alpha2 * x) + b2)))
        Y22 = (C12 + (C22 * (math.log((alpha2 * x2) + b2))) + (C32 * ((alpha2 * x2) + b2)) + (C42 / ((alpha2 * x2) + b2)) + (
                    (3 / 4) * R2*(math.log((alpha2 * x2) + b2)) * ((alpha2 * x2) + b2)))

        A3 = (h * A11)
        q3 = ((P3 * (A3 * cos1)) / L1)
        alpha3 = ((L1_p - b0) / L1)
        R3 = (q3 / ((alpha ** 4) * E * h))
        F3 = q3 * L1
        L3_p = L1_p
        b3 = b1
        C23 = (3 * R3 * b3)
        C43 = ((3 / 4) * R3 * (b3 * b3))
        C33 = (C43 / (L3_p * L3_p)) - (C23 / L3_p) - ((1.5) * R3 * ((math.log(L3_p)) + 1))
        C13 = -((C23 * (math.log(L3_p))) + (C33 * L3_p) + (C43 / L3_p) + (1.5 * R3 * (math.log(L3_p)) * L3_p))

        Y3 = -(C13 + (C23 * (math.log((alpha3 * x) + b3))) + (C33 * ((alpha3 * x) + b3)) + (C43 / ((alpha3 * x) + b3)) + (
                    1.5 * R3 * (math.log((alpha3 * x) + b0)) * ((alpha3 * x) + b3)))
        Y33 = -(C13 + (C23 * (math.log((alpha3 * x2) + b3))) + (C33 * ((alpha3 * x2) + b3)) + (
                    C43 / ((alpha3 * x2) + b3)) + (1.5 * R3 * (math.log((alpha3 * x2) + b0)) * ((alpha3 * x2) + b3)))

        K1 = abs(F1 / Y1);
        K2 = abs(F2 / Y2);
        K3 = abs(F3 / Y3);
        K11 = abs(F1 / Y11);
        K22 = abs(F2 / Y22);
        K33 = abs(F3 / Y33);
        Ks = (6.721*2)/(0.0508-b1-b2-b3)  * 2

        miu = 1;
        miu1 = 1;
        miu2 = 1;
        miu3 = 1;

        Yo2 = ((Y2) + ((Ks * Y1) / (K2 * (1 + (Ks / K1)))) - ((Ks * K3 * Y3) / (K2 * (K3 + Ks)))) / (
                    1 - ((Ks * Ks) / (K2 * (K3 + Ks))) - ((Ks * Ks) / (K2 * (K1 + Ks))) + ((2 * Ks) / K2));
        Yo22 = ((Y22) + ((Ks * Y11) / (K22 * (1 + (Ks / K11)))) - ((Ks * K33 * Y33) / (K22 * (K33 + Ks)))) / (
                    1 - ((Ks * Ks) / (K2 * (K33 + Ks))) - ((Ks * Ks) / (K22 * (K11 + Ks))) + ((2 * Ks) / K22));

        delta_L = ((19994.8 * 0.0508 * 0.0508 * L1) / ((2 * b0 + L2_3p) * E * h2))** 2
        Yo1 = -((((Ks / K1) * (Yo2)) + (Y1)) / (1 + (Ks / K1)))
        Yo11 = ((((Ks / K11) * (Yo22)) + (Y11)) / (1 + (Ks / K11)))
        Yo3 = -((((Ks / K3) * (Yo2)) - (Y3)) / (1 + (Ks / K3)))
        Yo33 = ((((Ks / K33) * (Yo22)) - (Y33)) / (1 + (Ks / K33)))
        print(Yo11,Yo22,Yo33)
        newdisplacement1 = 100*(Yo1 )-(Yo22+Yo33)/2
        + miu1 * 100 * ((max(0, (((((((Ks / K1) * (Yo2)) + (Y1)) / (1 + (Ks / K1))) - b0) + (Yo2 + b2 / 2)) - d1)))**2)
        + miu1 * 10 * ((max(0, -((2 * b0) + L2_3p - 0.03)))** 2)
        + miu3 * 100000000 * ((max(0, (0.0508 - ((2 * abs(b0)) + abs(b2) + abs(d1)))))** 2)
        + miu3 * 10000000 * ((max(0, -(0.0508 - ((2 * abs(b0)) + abs(b2) + abs(d1)))))** 2)
        + miu3 * 100000 * ((max(0, (0.0508 - ((2 * abs(L1_p)) + abs(L2_3p) + abs(d11)))))** 2)
        + miu * 100000 * ((max(0, (-(0.0508 - ((2 * abs(L1_p)) + abs(L2_3p) + abs(d11))))))** 2)
        + miu2 * 1000000000 * ((max(0, (-d1)))** 2)
        +miu * ((max(0, (-((b0 / L1_p) - 0.1))))** 2)
        + miu * 1000000 * ((max(0, ((b0 / L1_p) - 0.9)))** 2)
        +miu * ((max(0, (-((L2_3p / b2) - 0.1))))** 2) +miu * 1000000 *((max(0, ((L2_3p / b2) - 0.9)))** 2)

        if Yo11 >0.05 or Yo11<-0.05 or Yo22>0.05 or Yo22 <-0.05 or Yo33 >0.05 or Yo33 <-0.05:# or newdisplacement1 >1:
            Reward = ""
        newdisplacement2 = 0
        if newdisplacement1 > displacement1:
            newdisplacement2 =1

        total_x1 = b1+L2_2p
        total_x2 = L2_3p+L1_p
        Reward = newdisplacement1/100

        #  Update and display data
        Full_step_list = [L1, L1_p,L1_2p, L2,L2_2p, L2_3p,b1,b2,h,Yo1, Yo2, Yo3,Reward]
        Full_step_list = [ '%.4f' % elem for elem in Full_step_list ]
        print(Full_step_list)
        with open("Output_data_analytical.csv", "a") as out_file:
            out_String = str(Full_step_list)
            out_file.write(out_String + '\n')

        self.state = np.array([newdisplacement1, newdisplacement2])
        return self._get_obs(), Reward, True, {}

    def reset(self):
        high = np.array([1, 1])
        self.state = self.np_random.uniform(low=-high, high=high)
        self.last_u = None
        return self._get_obs()

    def _get_obs(self):
        displacement1, displacement2 = self.state
        return np.array([displacement1, displacement2])

    def _render_(self):
        print("Complete")
        return

