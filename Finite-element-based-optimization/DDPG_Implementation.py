#  Stable Baselines implementation of FE model based optimization of MDOF Soft Actuator

import gym

import numpy as np
import pandas as pd

from stable_baselines.ddpg.policies import MlpPolicy
from stable_baselines.common.noise import OrnsteinUhlenbeckActionNoise
from stable_baselines import DDPG


# Write data to csv file for post processing
with open("Output data.csv", "w") as out_file:
    out_String = str("a1x,a1y,a2x,a2y,a3x,a3y,a4x,a4y,max_x_disp,max_y_disp,Avgdispx,Avgdispy,Reward")
    out_file.write(out_String + '\n')

i = 0  # Main loop index

# Main Loop
while True:
    i = i+1
    df = pd.read_csv('Output data.csv')

    if i > 2:
        # Generate 1 full run of data on the second pull data for max and last entry
        last_val = df['Reward'].iloc[-1]
        max_disp = df['Reward'].max()

        # Main loop termination criteria
        if last_val >= max_disp:
            i = 10

    # Define Environment
    env = gym.make('Threechamber-v0')  # Placed in GYM folder in Lib

    # DDPG implementation
    n_actions = env.action_space.shape[-1]
    param_noise = None
    action_noise = OrnsteinUhlenbeckActionNoise(mean=np.zeros(n_actions), sigma=float(0.5) * np.ones(n_actions))
    model = DDPG(MlpPolicy, env, verbose=1, param_noise=param_noise, action_noise=action_noise, batch_size=16)
    model.learn(total_timesteps=500)
    model.save("test_optimization")

    del model # sep. consideration for each starting geometry
    model = DDPG.load("test_optimization")
    count = 0  # RL section count
    obs = env.reset()

    while True:
        action, _states = model.predict(obs)
        obs, rewards, dones, info = env.step(action)
        count = count +1
        break

    if i == 5:
        break
