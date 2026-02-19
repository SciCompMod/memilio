import memilio.simulation.jolly as jl
import os
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabc
import random

pyabc.settings.set_figure_params('pyabc')

simulation_filepath = "/home/kilian/Documents/projects/jolly/common_data/farms.csv"
observation_filepath = "/home/kilian/Documents/projects/jolly/first_phase/observations.csv"
tmax = 32
dt = 1

sei_duck = 1.7452
sei_chicken = 2.7603
ei_duck = 20.5745
ei_chicken = 0.6410
id_duck = 0.1011
id_chicken = 0.6412
death_duck = 0.00013
death_chicken = 0.00025

mean, sd = 10000, 2000
sigma = np.sqrt(np.log(1 + (sd / mean) ** 2))
mu = np.log(mean) - 0.5 * sigma ** 2


def model(parameter):
    seed = random.randint(0, 1000000000)
    result = jl.simulate(
        farm_file=simulation_filepath,
        tmax=tmax,
        dt=dt,
        suspicion_threshold=parameter["suspicion_threshold"],
        sensitivity=parameter["sensitivity"],
        h0=1.0,
        r0=parameter["r0"],
        alpha=parameter["alpha"],
        A0_SEI=sei_duck,
        A0_EI=ei_duck,
        A0_ID=id_duck,
        A0_DeathRate=death_duck,
        A1_SEI=sei_duck,
        A1_EI=ei_duck,
        A1_ID=id_duck,
        A1_DeathRate=death_duck,
        A2_SEI=sei_chicken,
        A2_EI=ei_chicken,
        A2_ID=id_chicken,
        A2_DeathRate=death_chicken,
        A3_SEI=sei_chicken,
        A3_EI=ei_chicken,
        A3_ID=id_chicken,
        A3_DeathRate=death_chicken,
        A4_SEI=sei_chicken,
        A4_EI=ei_chicken,
        A4_ID=id_chicken,
        A4_DeathRate=death_chicken,
        foi_inner_factor0=parameter["foi_inner_factor0"],
        foi_outer_factor0=parameter["foi_outer_factor0"],
        foi_inner_factor1=parameter["foi_inner_factor1"],
        foi_outer_factor1=parameter["foi_outer_factor1"],
        foi_inner_factor2=parameter["foi_inner_factor2"],
        foi_outer_factor2=parameter["foi_outer_factor2"],
        foi_inner_factor3=parameter["foi_inner_factor3"],
        foi_outer_factor3=parameter["foi_outer_factor3"],
        foi_inner_factor4=parameter["foi_inner_factor4"],
        foi_outer_factor4=parameter["foi_outer_factor4"],
        first_infection_day=parameter["first_infection_day"],
        second_infection_day=parameter["second_infection_day"],
        third_infection_day=parameter["third_infection_day"],
        seed=seed)
    return {"data": np.array(result), "seed": seed}


prior = pyabc.Distribution(
    suspicion_threshold=pyabc.RV("uniform", 0.1, 0.4),
    sensitivity=pyabc.RV("uniform", 0.95, 0.05),
    # h0=pyabc.RV("uniform", 0, 100),
    r0=pyabc.RV("lognorm", s=sigma, scale=np.exp(mu)),
    alpha=pyabc.RV("uniform", 0, 40),
    # A0_SEI = pyabc.RV("uniform", 0,1),
    # A0_EI = pyabc.RV("uniform", 0,1),
    # A0_ID = pyabc.RV("uniform", 0,1),
    # A0_DeathRate = pyabc.RV("uniform", 0,1),
    # A1_SEI = pyabc.RV("uniform", 0,1),
    # A1_EI = pyabc.RV("uniform", 0,1),
    # A1_ID = pyabc.RV("uniform", 0,1),
    # A1_DeathRate = pyabc.RV("uniform", 0,1),
    # A2_SEI = pyabc.RV("uniform", 0,1),
    # A2_EI = pyabc.RV("uniform", 0,1),
    # A2_ID = pyabc.RV("uniform", 0,1),
    # A2_DeathRate = pyabc.RV("uniform", 0,1),
    # A3_SEI = pyabc.RV("uniform", 0,1),
    # A3_EI = pyabc.RV("uniform", 0,1),
    # A3_ID = pyabc.RV("uniform", 0,1),
    # A3_DeathRate = pyabc.RV("uniform", 0,1),
    # A4_SEI = pyabc.RV("uniform", 0,1),
    # A4_EI = pyabc.RV("uniform", 0,1),
    # A4_ID = pyabc.RV("uniform", 0,1),
    # A4_DeathRate = pyabc.RV("uniform", 0,1),
    foi_inner_factor0=pyabc.RV("uniform", 0, 100),
    foi_outer_factor0=pyabc.RV("uniform", 0, 1),
    foi_inner_factor1=pyabc.RV("uniform", 0, 100),
    foi_outer_factor1=pyabc.RV("uniform", 0, 1),
    foi_inner_factor2=pyabc.RV("uniform", 0, 100),
    foi_outer_factor2=pyabc.RV("uniform", 0, 1),
    foi_inner_factor3=pyabc.RV("uniform", 0, 100),
    foi_outer_factor3=pyabc.RV("uniform", 0, 1),
    foi_inner_factor4=pyabc.RV("uniform", 0, 100),
    foi_outer_factor4=pyabc.RV("uniform", 0, 1),
    first_infection_day=pyabc.RV("uniform", 0, 12),
    second_infection_day=pyabc.RV("uniform", 0, 14),
    third_infection_day=pyabc.RV("uniform", 0, 14)
)

observation_data = pd.read_csv(observation_filepath)
obs_data = dict(data=observation_data["date_confirmed"].values)


def distance_not_confirmed(data, observation):
    obs = data['data']
    obs_idx = obs == -1
    sim = observation['data']
    sim_idx = sim == -1
    return np.sum(np.abs(sim_idx.astype(int) - obs_idx.astype(int))) / len(obs)


def distance_confirmed(data, observation):
    obs = data['data']
    obs_idx = obs != -1
    sim = observation['data']
    sim_idx = sim != -1
    joint_idx = obs_idx & sim_idx
    return np.sum(np.abs(obs[joint_idx] - sim[joint_idx])) / (np.sum(obs_idx) * tmax)


distance = pyabc.AdaptiveAggregatedDistance(
    [distance_not_confirmed, distance_confirmed], adaptive=False)

abc = pyabc.ABCSMC(model, prior, distance, population_size=100,
                   transitions=pyabc.transition.LocalTransition(20, 5))
db_path = "sqlite:///" + os.path.join(tempfile.gettempdir(), "tmp3.db")
abc.new(db_path, obs_data)
history = abc.run(max_nr_populations=6, minimum_epsilon=0.1)

# print(model(0.1))
