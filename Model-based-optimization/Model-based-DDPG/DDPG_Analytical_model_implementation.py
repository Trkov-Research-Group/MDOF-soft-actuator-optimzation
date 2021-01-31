#  Nicholas Pagliocca
#  DDPG implementation using the fork from OpenAI baselines called Stable baselines

import gym
import numpy as np
import pandas as pd
from stable_baselines.ddpg.policies import MlpPolicy
from stable_baselines.common.noise import OrnsteinUhlenbeckActionNoise
from stable_baselines import DDPG

with open("Output_data_analytical.csv", "w") as out_file:
    out_String = str( "L1,L1_p,L1_2p,L2,L1_2p,L2_3p,b1,L2_pb2,h,Y01,Y02,Y03,Reward")
    out_file.write(out_String + '\n')

i = 0

while True:
    i = i+1
    df = pd.read_csv('Output_data_analytical.csv')

    max_disp = df['Reward'].max()

    if i > 2:

        last_val = df['Reward'].iloc[-1]
        n_last = last_val
        print(n_last)
        if last_val >= max_disp:
            i = 5
    print(i)
    print(max_disp)

    env = gym.make('Threechambera-v0')

    # the noise objects for DDPG
    n_actions = env.action_space.shape[-1]
    param_noise = None
    action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(n_actions), sigma=float(0.5) * np.ones(n_actions))

    model = DDPG(MlpPolicy, env, verbose=1, param_noise=param_noise, action_noise=action_noise, batch_size=8)
    model.learn(total_timesteps=500)
    model.save("test_analytical")

    #del model # remove to demonstrate saving and loading

    model = DDPG.load("test_analytical")
    count = 0
    obs = env.reset()
    while True:
        action, _states = model.predict(obs)
        obs, rewards, dones, info = env.step(action)
        count = count +1
        break
    if i == 5:
        break



